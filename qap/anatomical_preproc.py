base_test_dir = "/tdata/QAP/qc_test"


def anatomical_reorient_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_scan

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.maths as fsl

    from nipype.interfaces.afni import preprocess

    from workflow_utils import check_input_resources
    
    
    check_input_resources(resource_pool, "anatomical_scan")


    anat_deoblique = pe.Node(interface=preprocess.Refit(),
                                name='anat_deoblique')

    anat_deoblique.inputs.in_file = resource_pool["anatomical_scan"]
    anat_deoblique.inputs.deoblique = True


    anat_reorient = pe.Node(interface=preprocess.Resample(),
                            name='anat_reorient')

    anat_reorient.inputs.orientation = 'RPI'
    anat_reorient.inputs.outputtype = 'NIFTI_GZ'


    workflow.connect(anat_deoblique, 'out_file', anat_reorient, 'in_file')

   
    resource_pool["anatomical_reorient"] = (anat_reorient, 'out_file')


    return workflow, resource_pool



def run_anatomical_reorient(anatomical_scan, run=True):

    # stand-alone runner for anatomical reorient workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='anatomical_reorient_workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "anatomical_reorient")
    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1


    resource_pool["anatomical_scan"] = anatomical_scan
    
    workflow, resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_anatomical_reorient')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["anatomical_reorient"]

    workflow.connect(node, out_file, ds, 'anatomical_reorient')


    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "anatomical_reorient",\
                                         "*"))[0]

        return outpath

    else:

        return workflow, workflow.base_dir



def anatomical_skullstrip_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_reorient

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util

    from nipype.interfaces.afni import preprocess


    if "anatomical_reorient" not in resource_pool.keys():
        
        from anatomical_preproc import anatomical_reorient_workflow

        workflow, resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config)


    #check_input_resources(resource_pool, "anatomical_reorient")


    anat_skullstrip = pe.Node(interface=preprocess.SkullStrip(),
                              name='anat_skullstrip')
    
    anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'
    

    anat_skullstrip_orig_vol = pe.Node(interface=preprocess.Calc(),
                                       name='anat_skullstrip_orig_vol')

    anat_skullstrip_orig_vol.inputs.expr = 'a*step(b)'
    anat_skullstrip_orig_vol.inputs.outputtype = 'NIFTI_GZ'


    if len(resource_pool["anatomical_reorient"]) == 2:
        node, out_file = resource_pool["anatomical_reorient"]
        workflow.connect(node, out_file, anat_skullstrip, 'in_file')
    else:
        anat_skullstrip.inputs.in_file = \
            resource_pool["anatomical_reorient"]


    if len(resource_pool["anatomical_reorient"]) == 2:
        node, out_file = resource_pool["anatomical_reorient"]
        workflow.connect(node, out_file,
                             anat_skullstrip_orig_vol, 'in_file_a')
    else:
        anat_skullstrip_orig_vol.inputs.in_file_a = \
            resource_pool["anatomical_reorient"]


    workflow.connect(anat_skullstrip, 'out_file',
                        anat_skullstrip_orig_vol, 'in_file_b')


    resource_pool["anatomical_brain"] = (anat_skullstrip_orig_vol, 'out_file')


    return workflow, resource_pool



def run_anatomical_skullstrip(anatomical_reorient, run=True):

    # stand-alone runner for anatomical skullstrip workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='anatomical_skullstrip_workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "anatomical_skullstrip")
    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1


    resource_pool["anatomical_reorient"] = anatomical_reorient
    
    workflow, resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_anatomical_skullstrip')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["anatomical_brain"]

    workflow.connect(node, out_file, ds, 'anatomical_brain')

    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "anatomical_brain", \
                                         "*"))[0]

        return outpath

    else:

        return workflow, workflow.base_dir



def flirt_anatomical_linear_registration(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_brain

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util

    import nipype.interfaces.fsl as fsl

    from workflow_utils import check_input_resources, \
                               check_config_settings

    from nipype.interfaces.fsl.base import Info
    
    if "template_brain_for_anat" not in config:
        config["template_brain_for_anat"] = Info.standard_image("MNI152_T1_2mm_brain.nii.gz")
    check_config_settings(config, "template_brain_for_anat")


    if "anatomical_brain" not in resource_pool.keys():

        from anatomical_preproc import anatomical_skullstrip_workflow

        workflow, resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config)


    #check_input_resources(resource_pool, "anatomical_brain")

    calc_flirt_warp = pe.Node(interface=fsl.FLIRT(), name='calc_flirt_warp')

    calc_flirt_warp.inputs.cost = 'corratio'


    if len(resource_pool["anatomical_brain"]) == 2:
        node, out_file = resource_pool["anatomical_brain"]
        workflow.connect(node, out_file, calc_flirt_warp, 'in_file')
    else:
        calc_flirt_warp.inputs.in_file = resource_pool["anatomical_brain"]


    calc_flirt_warp.inputs.reference = config["template_brain_for_anat"]


    resource_pool["flirt_affine_xfm"] = (calc_flirt_warp, 'out_matrix_file')

    resource_pool["flirt_linear_warped_image"] = (calc_flirt_warp, 'out_file')


    return workflow, resource_pool



def run_flirt_anatomical_linear_registration(anatomical_brain, \
                                                 template_brain, run=True):

    # stand-alone runner for FSL FLIRT anatomical linear registration workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='flirt_anatomical_linear_registration_' \
                                'workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "flirt_anatomical_linear_" \
                                    "registration")
    workflow.base_dir = workflow_dir


    num_cores_per_subject = 1


    resource_pool = {}
    config = {}
    
    
    resource_pool["anatomical_brain"] = anatomical_brain
    config["template_brain_for_anat"] = template_brain
    
    workflow, resource_pool = \
        flirt_anatomical_linear_registration(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_flirt_anatomical_linear_' \
                                      'registration')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["flirt_linear_warped_image"]

    workflow.connect(node, out_file, ds, 'flirt_linear_warped_image')

    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})


        outpath = glob.glob(os.path.join(workflow_dir, "flirt_linear_" \
                                         "warped_image", "*"))[0]

        return outpath

    else:

        return workflow, workflow.base_dir



def ants_anatomical_linear_registration(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_brain

    # linear ANTS registration takes roughly 2.5 minutes per subject running
    # on one core of an Intel Core i7-4800MQ CPU @ 2.70GHz

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    

    import nipype.interfaces.utility as util

    from anatomical_preproc_utils import ants_lin_reg, \
                                         separate_warps_list

    from workflow_utils import check_input_resources, \
                               check_config_settings
    from nipype.interfaces.fsl.base import Info
    
    if "template_brain_for_anat" not in config:
        config["template_brain_for_anat"] = Info.standard_image("MNI152_T1_2mm_brain.nii.gz")
    check_config_settings(config, "template_brain_for_anat")


    if "anatomical_brain" not in resource_pool.keys():

        from anatomical_preproc import anatomical_skullstrip_workflow

        workflow, resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config)


    #check_input_resources(resource_pool, "anatomical_brain")


    calc_ants_warp = pe.Node(interface=util.Function(
                                 input_names=['anatomical_brain',
                                              'reference_brain'],
                                 output_names=['warp_list',
                                               'warped_image'],
                                 function=ants_lin_reg),
                                 name='calc_ants_linear_warp')


    select_forward_initial = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=separate_warps_list), name='select_forward_initial')

    select_forward_initial.inputs.selection = "Initial"


    select_forward_rigid = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=separate_warps_list), name='select_forward_rigid')

    select_forward_rigid.inputs.selection = "Rigid"


    select_forward_affine = pe.Node(util.Function(input_names=['warp_list',
            'selection'], output_names=['selected_warp'],
            function=separate_warps_list), name='select_forward_affine')

    select_forward_affine.inputs.selection = "Affine"


    if len(resource_pool["anatomical_brain"]) == 2:
        node, out_file = resource_pool["anatomical_brain"]
        workflow.connect(node, out_file, calc_ants_warp, 'anatomical_brain')
    else:
       calc_ants_warp.inputs.anatomical_brain = \
            resource_pool["anatomical_brain"]


    calc_ants_warp.inputs.reference_brain = config["template_brain_for_anat"]


    workflow.connect(calc_ants_warp, 'warp_list',
                         select_forward_initial, 'warp_list')

    workflow.connect(calc_ants_warp, 'warp_list',
                         select_forward_rigid, 'warp_list')

    workflow.connect(calc_ants_warp, 'warp_list',
                         select_forward_affine, 'warp_list')


    resource_pool["ants_initial_xfm"] = \
        (select_forward_initial, 'selected_warp')

    resource_pool["ants_rigid_xfm"] = (select_forward_rigid, 'selected_warp')

    resource_pool["ants_affine_xfm"] = \
        (select_forward_affine, 'selected_warp')

    resource_pool["ants_linear_warped_image"] = \
        (calc_ants_warp, 'warped_image')


    return workflow, resource_pool
    
    
    
def run_ants_anatomical_linear_registration(anatomical_brain, \
                                                template_brain, num_cores=1, \
                                                run=True):

    # stand-alone runner for anatomical skullstrip workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='ants_anatomical_linear_registration_' \
                                'workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "ants_anatomical_linear_" \
                                    "registration")
    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1
    
    
    os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] = str(num_cores)


    resource_pool["anatomical_brain"] = anatomical_brain
    config["template_brain_for_anat"] = template_brain
    
    workflow, resource_pool = \
          ants_anatomical_linear_registration(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_ants_anatomical_linear_' \
                                      'registration')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["ants_linear_warped_image"]

    workflow.connect(node, out_file, ds, 'ants_linear_warped_image')

    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "ants_linear_warped_" \
                                         "image", "*"))[0]

        return outpath

    else:

        return workflow



def segmentation_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_brain

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.ants as ants
    import nipype.interfaces.utility as util

    from anatomical_preproc_utils import pick_seg_type

    from workflow_utils import check_input_resources, \
                               check_config_settings


    if "anatomical_brain" not in resource_pool.keys():

        from anatomical_preproc import anatomical_skullstrip_workflow

        workflow, resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config)


    segment = pe.Node(interface=fsl.FAST(), name='segmentation')

    segment.inputs.img_type = 1
    segment.inputs.segments = True
    segment.inputs.probability_maps = True
    segment.inputs.out_basename = 'segment'


    if len(resource_pool["anatomical_brain"]) == 2:
        node, out_file = resource_pool["anatomical_brain"]
        workflow.connect(node, out_file, segment, 'in_files')
    else:
        segment.inputs.in_files = resource_pool["anatomical_brain"]


    # process masks

    seg_types = ["gm", "wm", "csf"]

    for seg in seg_types:

        pick_seg = pe.Node(interface=util.Function(
                           input_names=['probability_maps',
                                        'seg_type'],
                           output_names=['filename'],
                           function=pick_seg_type),
                           name='pick_%s' % seg)

        pick_seg.inputs.seg_type = seg
        

        workflow.connect(segment, 'tissue_class_files',
                             pick_seg, 'probability_maps')

        
        resource_pool["anatomical_%s_mask" % seg] = (pick_seg, 'filename')


    return workflow, resource_pool



def run_segmentation_workflow(anatomical_brain, run=True):

    # stand-alone runner for segmentation workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='segmentation_workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "segmentation")
    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1


    resource_pool["anatomical_brain"] = anatomical_brain  
    
    workflow, resource_pool = \
            segmentation_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_segmentation')
    ds.inputs.base_directory = workflow_dir
    
    
    seg_types = ["gm", "wm", "csf"]

    for seg in seg_types:
    
        node, out_file = resource_pool["anatomical_%s_mask" % seg]

        workflow.connect(node, out_file, ds, 'anatomical_%s_mask' % seg)


    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "anatomical_*_mask", \
                                         "*"))

        return outpath

    else:

        return workflow, workflow.base_dir

