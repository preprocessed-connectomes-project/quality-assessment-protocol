base_test_dir = "/tdata/QAP/qc_test"


def anatomical_reorient_workflow(workflow, resource_pool, config, name="_"):

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
                                name='anat_deoblique%s' % name)

    anat_deoblique.inputs.in_file = resource_pool["anatomical_scan"]
    anat_deoblique.inputs.deoblique = True


    anat_reorient = pe.Node(interface=preprocess.Resample(),
                            name='anat_reorient%s' % name)

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
        try:
            workflow.run(plugin='ResourceMultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "anatomical_reorient",\
                                             "*"))[0]
            return outpath
        except:
            workflow.run(plugin='MultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "anatomical_reorient",\
                                             "*"))[0]
            return outpath

    else:

        return workflow, workflow.base_dir



def anatomical_skullstrip_workflow(workflow, resource_pool, config, name="_"):

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
            anatomical_reorient_workflow(workflow, resource_pool, config, name)


    #check_input_resources(resource_pool, "anatomical_reorient")


    anat_skullstrip = pe.Node(interface=preprocess.SkullStrip(),
                              name='anat_skullstrip%s' % name)
    
    anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'
    

    anat_skullstrip_orig_vol = pe.Node(interface=preprocess.Calc(),
                                       name='anat_skullstrip_orig_vol%s' % name)

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
        try:
            workflow.run(plugin='ResourceMultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "anatomical_brain", \
                                             "*"))[0]
            return outpath
        except:
            workflow.run(plugin='MultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "anatomical_brain", \
                                             "*"))[0]
            return outpath

    else:

        return workflow, workflow.base_dir



def afni_anatomical_linear_registration(workflow, resource_pool, \
    config, name="_"):

    # resource pool should have:
    #     anatomical_brain

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.afni as afni

    from workflow_utils import check_input_resources, \
                               check_config_settings

    check_config_settings(config, "template_skull_for_anat")


    if "anatomical_reorient" not in resource_pool.keys():

        from anatomical_preproc import anatomical_reorient_workflow

        workflow, resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config, name)

    calc_3dallineate_warp = pe.Node(interface=afni.Allineate(),
                                    name='calc_3dAllineate_warp%s' % name)
    calc_3dallineate_warp.inputs.outputtype = "NIFTI_GZ"

    if len(resource_pool["anatomical_reorient"]) == 2:
        node, out_file = resource_pool["anatomical_reorient"]
        workflow.connect(node, out_file, calc_3dallineate_warp, 'in_file')
    else:
        calc_3dallineate_warp.inputs.in_file = \
            resource_pool["anatomical_reorient"]


    calc_3dallineate_warp.inputs.reference = config["template_skull_for_anat"]

    calc_3dallineate_warp.inputs.out_file = "allineate_warped_brain.nii.gz"
    calc_3dallineate_warp.inputs.out_matrix = "3dallineate_warp"


    resource_pool["allineate_linear_xfm"] = \
        (calc_3dallineate_warp, 'matrix')

    resource_pool["afni_linear_warped_image"] = \
        (calc_3dallineate_warp, 'out_file')


    return workflow, resource_pool



def run_afni_anatomical_linear_registration(reference_skull,
                                            anatomical_skull=None,
                                            anatomical_scan=None,
                                            run=True):

    # stand-alone runner for anatomical skullstrip workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    from workflow_utils import raise_smart_exception

    workflow = pe.Workflow(name='3dallineate_workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "3dallineate")
    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    if anatomical_skull:
        resource_pool["anatomical_reorient"] = anatomical_skull
    elif anatomical_scan:
        resource_pool["anatomical_scan"] = anatomical_scan
    else:
        err = "No anatomical_reorient or anatomical_scan provided!"
        raise_smart_exception(locals(),err)

    config["template_skull_for_anat"] = reference_skull
    
    workflow, resource_pool = \
            afni_anatomical_linear_registration(workflow, resource_pool, \
                                                config)


    ds = pe.Node(nio.DataSink(), name='datasink_3dallineate')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["afni_linear_warped_image"]
    workflow.connect(node, out_file, ds, 'afni_linear_warped_image')

    node, out_file = resource_pool["allineate_linear_xfm"]
    workflow.connect(node, out_file, ds, 'allineate_linear_xfm')

    if run == True:
        try:
            workflow.run(plugin='ResourceMultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, \
                                             "afni_linear_warped_image", \
                                             "*"))[0]
            return outpath
        except:
            workflow.run(plugin='MultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, \
                                             "afni_linear_warped_image", \
                                             "*"))[0]
            return outpath

    else:

        return workflow, workflow.base_dir



def afni_segmentation_workflow(workflow, resource_pool, config, name="_"):

    # resource pool should have:
    #     anatomical_brain

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    from nipype.interfaces.afni import preprocess

    from workflow_utils import check_input_resources, \
                               check_config_settings


    if "anatomical_brain" not in resource_pool.keys():

        from anatomical_preproc import anatomical_skullstrip_workflow

        workflow, resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config, name)


    segment = pe.Node(interface=preprocess.Seg(), name='segmentation%s' % name)

    segment.inputs.mask = 'AUTO'

    if len(resource_pool["anatomical_brain"]) == 2:
        node, out_file = resource_pool["anatomical_brain"]
        workflow.connect(node, out_file, segment, 'in_file')
    else:
        segment.inputs.in_file = resource_pool["anatomical_brain"]

    # output processing
    AFNItoNIFTI = pe.Node(interface=preprocess.AFNItoNIFTI(),
                          name="segment_AFNItoNIFTI%s" % name)

    AFNItoNIFTI.inputs.out_file = "classes.nii.gz"
    

    workflow.connect(segment, 'out_file', AFNItoNIFTI, 'in_file')

    # break out each of the three tissue types into
    # three separate NIFTI files
    extract_CSF = pe.Node(interface=preprocess.Calc(),
                          name='extract_CSF_mask%s' % name)
    extract_CSF.inputs.expr = "within(a,1,1)"
    extract_CSF.inputs.out_file = "anatomical_csf_mask.nii.gz"

    extract_GM = pe.Node(interface=preprocess.Calc(),
                          name='extract_GM_mask%s' % name)
    extract_GM.inputs.expr = "within(a,2,2)"
    extract_GM.inputs.out_file = "anatomical_gm_mask.nii.gz"

    extract_WM = pe.Node(interface=preprocess.Calc(),
                          name='extract_WM_mask%s' % name)
    extract_WM.inputs.expr = "within(a,3,3)"
    extract_WM.inputs.out_file = "anatomical_wm_mask.nii.gz"

    workflow.connect(AFNItoNIFTI, 'out_file', extract_CSF, 'in_file_a')
    workflow.connect(AFNItoNIFTI, 'out_file', extract_GM, 'in_file_a')
    workflow.connect(AFNItoNIFTI, 'out_file', extract_WM, 'in_file_a')


    resource_pool["anatomical_csf_mask"] = (extract_CSF, 'out_file')
    resource_pool["anatomical_gm_mask"] = (extract_GM, 'out_file')
    resource_pool["anatomical_wm_mask"] = (extract_WM, 'out_file')


    return workflow, resource_pool



def run_afni_segmentation_workflow(anatomical_brain, run=True):

    # stand-alone runner for segmentation workflow

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='afni_segmentation_workflow')

    current_dir = os.getcwd()

    workflow_dir = os.path.join(current_dir, "afni_segmentation")
    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1


    resource_pool["anatomical_brain"] = anatomical_brain  
    
    workflow, resource_pool = \
            afni_segmentation_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_afni_segmentation')
    ds.inputs.base_directory = workflow_dir
    
    
    seg_types = ["gm", "wm", "csf"]

    for seg in seg_types:
    
        node, out_file = resource_pool["anatomical_%s_mask" % seg]

        workflow.connect(node, out_file, ds, 'anatomical_%s_mask' % seg)


    if run == True:
        try:
            workflow.run(plugin='ResourceMultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "anatomical_*_mask", \
                                             "*"))            
        except:
            workflow.run(plugin='MultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "anatomical_*_mask", \
                                             "*"))

        return outpath

    else:

        return workflow, workflow.base_dir

