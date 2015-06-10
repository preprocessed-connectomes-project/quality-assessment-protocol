
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



def run_anatomical_reorient(anatomical_scan):

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


    workflow.run(plugin='MultiProc', plugin_args= \
                     {'n_procs': num_cores_per_subject})


    outpath = glob.glob(os.path.join(workflow_dir, "anatomical_reorient", \
                                     "*"))[0]


    return outpath



def test_run_anatomical_reorient():

    import os
    import commands


    if "anatomical_reorient" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_scan = os.path.join(base_test_dir, \
                    "raw_data/0050002/session_1/anat_1/mprage.nii.gz")

    ref_out = os.path.join(base_test_dir, \
                  "ABIDE/0050002_session_1/anatomical_reorient/" \
                  "mprage_RPI.nii.gz")

    output = run_anatomical_reorient(anat_scan)

    cmd = "3ddot -demean %s %s" % (output, ref_out)

    correlation = commands.getoutput(cmd)

    os.system("rm -R anatomical_reorient")


    assert correlation > 0.98



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



def run_anatomical_skullstrip(anatomical_reorient):

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


    workflow.run(plugin='MultiProc', plugin_args= \
                     {'n_procs': num_cores_per_subject})


    outpath = glob.glob(os.path.join(workflow_dir, "anatomical_brain", \
                                     "*"))[0]


    return outpath



def test_run_anatomical_skullstrip():

    import os
    import commands


    if "anatomical_skullstrip" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_reorient = os.path.join(base_test_dir, \
                        "ABIDE/0050002_session_1/anatomical_reorient/" \
                        "mprage_RPI.nii.gz")

    ref_out = os.path.join(base_test_dir, \
                  "ABIDE/0050002_session_1/anatomical_brain/" \
                  "mprage_RPI_calc.nii.gz")

    output = run_anatomical_skullstrip(anat_reorient)

    cmd = "3ddot -demean %s %s" % (output, ref_out)

    correlation = commands.getoutput(cmd)

    os.system("rm -R anatomical_brain")


    assert correlation > 0.98



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

    #workflow.connect(calc_ants_warp, 'warp_list',
    #                     select_forward_warp, 'warp_list')

    #workflow.connect(calc_ants_warp, 'warp_list',
    #                     select_inverse_warp, 'warp_list')


    resource_pool["ants_initial_xfm"] = \
        (select_forward_initial, 'selected_warp')

    resource_pool["ants_rigid_xfm"] = (select_forward_rigid, 'selected_warp')

    resource_pool["ants_affine_xfm"] = \
        (select_forward_affine, 'selected_warp')

    resource_pool["ants_linear_warped_image"] = \
        (calc_ants_warp, 'warped_image')


    return workflow, resource_pool



def segmentation_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_brain
    #     ants_initial_xfm
    #     ants_rigid_xfm
    #     ants_affine_xfm

    # configuration should have:
    #     prior_csf
    #     prior_gm
    #     prior_wm
    #     csf_threshold
    #     gm_threshold
    #     wm_threshold


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


    if ("ants_affine_xfm" not in resource_pool.keys()) or \
           ("ants_rigid_xfm" not in resource_pool.keys()) or \
               ("ants_initial_xfm" not in resource_pool.keys()):

        from anatomical_preproc import ants_anatomical_linear_registration

        workflow, resource_pool = \
            ants_anatomical_linear_registration(workflow, resource_pool,
                                                config)


    if "anatomical_brain" not in resource_pool.keys():

        from anatomical_preproc import anatomical_skullstrip_workflow

        workflow, resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config)


    #check_input_resources(resource_pool, "ants_initial_xfm")
    #check_input_resources(resource_pool, "ants_rigid_xfm")
    #check_input_resources(resource_pool, "ants_affine_xfm")
    #check_input_resources(resource_pool, "anatomical_brain")

    check_config_settings(config, "prior_csf")
    check_config_settings(config, "prior_gm")
    check_config_settings(config, "prior_wm")


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

    def form_threshold_string(threshold):
        return '-thr %f -bin ' % (threshold)

    seg_types = ["gm", "wm", "csf"]

    for seg in seg_types:

        collect_linear_transforms = pe.Node(util.Merge(3), 
                                            name='%s_collect_xfm' % seg)


        tissueprior_mni_to_t1 = pe.Node(interface=ants.ApplyTransforms(),
                                        name='%s_prior_mni_to_t1' % seg)

        tissueprior_mni_to_t1.inputs.invert_transform_flags = \
            [True, True, True]
        tissueprior_mni_to_t1.inputs.interpolation = 'NearestNeighbor'


        pick_seg = pe.Node(interface=util.Function(
                           input_names=['probability_maps',
                                        'seg_type'],
                           output_names=['filename'],
                           function=pick_seg_type),
                           name='pick_%s' % seg)


        overlap_segmentmap_with_prior = pe.Node(interface= \
                                                fsl.MultiImageMaths(),
                                       name='overlap_%s_map_with_prior' % seg)

        overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '


        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(),
                                           name='binarize_threshold_%s' % seg)


        segment_mask = pe.Node(interface=fsl.MultiImageMaths(),
                               name='%s_mask' % seg)

        segment_mask.inputs.op_string = '-mas %s '


        if len(resource_pool["anatomical_brain"]) == 2:
            node, out_file = resource_pool["anatomical_brain"]
            workflow.connect(node, out_file,
                                 tissueprior_mni_to_t1, 'reference_image')
        else:
            tissueprior_mni_to_t1.inputs.reference_image = \
                resource_pool["anatomical_brain"]


        if len(resource_pool["ants_initial_xfm"]) == 2:
            node, out_file = resource_pool["ants_initial_xfm"]
            workflow.connect(node, out_file,
                                 collect_linear_transforms, 'in1')
        else:
            collect_linear_transforms.inputs.in1 = \
                resource_pool["ants_initial_xfm"]


        if len(resource_pool["ants_rigid_xfm"]) == 2:
            node, out_file = resource_pool["ants_rigid_xfm"]
            workflow.connect(node, out_file,
                                 collect_linear_transforms, 'in2')
        else:
            collect_linear_transforms.inputs.in2 = \
                resource_pool["ants_rigid_xfm"]


        if len(resource_pool["ants_affine_xfm"]) == 2:
            node, out_file = resource_pool["ants_affine_xfm"]
            workflow.connect(node, out_file,
                                 collect_linear_transforms, 'in3')
        else:
            collect_linear_transforms.inputs.in3 = \
                resource_pool["ants_affine_xfm"]


        tissueprior_mni_to_t1.inputs.input_image = config["prior_%s" % seg]

        binarize_threshold_segmentmap.inputs.op_string = \
            form_threshold_string(config["%s_threshold" % seg])

        pick_seg.inputs.seg_type = seg
        


        workflow.connect(segment, 'probability_maps',
                             pick_seg, 'probability_maps')


        workflow.connect(pick_seg, 'filename',
                             overlap_segmentmap_with_prior, 'in_file')


        workflow.connect(collect_linear_transforms, 'out',
                             tissueprior_mni_to_t1, 'transforms')


        # overlapping
        workflow.connect(tissueprior_mni_to_t1, 'output_image',
                             overlap_segmentmap_with_prior, 'operand_files')


        # binarize
        workflow.connect(overlap_segmentmap_with_prior, 'out_file',
                             binarize_threshold_segmentmap, 'in_file')


        # create segment mask
        workflow.connect(binarize_threshold_segmentmap, 'out_file',
                             segment_mask, 'in_file')

        workflow.connect(tissueprior_mni_to_t1, 'output_image',
                             segment_mask, 'operand_files')


        resource_pool["anatomical_%s_mask" % seg] = (segment_mask, 'out_file')


    return workflow, resource_pool



def run_segmentation_workflow(anatomical_brain):

    pass






