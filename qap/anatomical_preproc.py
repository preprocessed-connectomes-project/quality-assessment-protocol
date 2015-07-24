
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


    resource_pool["ants_initial_xfm"] = \
        (select_forward_initial, 'selected_warp')

    resource_pool["ants_rigid_xfm"] = (select_forward_rigid, 'selected_warp')

    resource_pool["ants_affine_xfm"] = \
        (select_forward_affine, 'selected_warp')

    resource_pool["ants_linear_warped_image"] = \
        (calc_ants_warp, 'warped_image')


    return workflow, resource_pool
    
    
    
def run_ants_anatomical_linear_registration(anatomical_brain, \
                                                template_brain, num_cores=1):

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


    workflow.run(plugin='MultiProc', plugin_args= \
                     {'n_procs': num_cores_per_subject})


    outpath = glob.glob(os.path.join(workflow_dir, "ants_linear_warped_" \
                                     "image", "*"))[0]


    return outpath



def segmentation_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     anatomical_brain

    # configuration should have:
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
        

        if config["process_segmentation_maps"] == True:

            workflow.connect(segment, 'probability_maps', \
                                 pick_seg, 'probability_maps')

            resource_pool["%s_probability_map" % seg] = (pick_seg, 'filename')

        else:

            workflow.connect(segment, 'tissue_class_files',
                                 pick_seg, 'probability_maps')

            resource_pool["anatomical_%s_mask" % seg] = (pick_seg, 'filename')


    return workflow, resource_pool



def run_segmentation_workflow(anatomical_brain, csf_threshold, gm_threshold, \
                              wm_threshold):

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
    
    config["csf_threshold"] = csf_threshold
    config["gm_threshold"] = gm_threshold
    config["wm_threshold"] = wm_threshold
    
    
    workflow, resource_pool = \
            segmentation_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_segmentation')
    ds.inputs.base_directory = workflow_dir
    
    
    seg_types = ["gm", "wm", "csf"]

    for seg in seg_types:
    
        node, out_file = resource_pool["anatomical_%s_mask" % seg]

        workflow.connect(node, out_file, ds, 'anatomical_%s_mask' % seg)


    workflow.run(plugin='MultiProc', plugin_args= \
                     {'n_procs': num_cores_per_subject})


    outpath = glob.glob(os.path.join(workflow_dir, "anatomical_*_mask", \
                                     "*"))


    return outpath



def process_seg_maps_workflow(workflow, resource_pool, config):

    # resource pool should have:
    #     csf_probability_map
    #     gm_probability_map
    #     wm_probability_map

    # configuration should have:
    #     csf_threshold
    #     gm_threshold
    #     wm_threshold
    #     csf_prior
    #     gm_prior
    #     wm_prior


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


    # for the "binarize_threshold_segmentmap" node
    def form_threshold_string(threshold):
        return '-thr %f -bin ' % (threshold)


    if ("csf_probability_map" not in resource_pool.keys()) or \
        ("gm_probability_map" not in resource_pool.keys()) or \
            ("wm_probability_map" not in resource_pool.keys()):

        from anatomical_preproc import segmentation_workflow

        # this should already be set if this workflow is running, 
        # but let's be safe for now
        config["process_segmentation_maps"] = True

        workflow, resource_pool = \
            segmentation_workflow(workflow, resource_pool, config)


    if ("csf_prior_standard_to_subject_warped" not in resource_pool.keys()) or \
        ("gm_prior_standard_to_subject_warped" not in resource_pool.keys()) or \
            ("wm_prior_standard_to_subject_warped" not in resource_pool.keys()):

        seg_types = ["csf", "gm", "wm"]

        for seg in seg_types:

            ''' for now '''
            config["registration_option"] = "ANTS"

            if config["registration_option"] == "ANTS":

                if ("ants_initial_xfm" not in resource_pool.keys()) or \
                    ("ants_rigid_xfm" not in resource_pool.keys()) or \
                        ("ants_affine_xfm" not in resource_pool.keys()):

                    from anatomical_preproc import \
                        ants_anatomical_linear_registration

                    workflow, resource_pool = \
                        ants_anatomical_linear_registration(workflow, \
                            resource_pool, config)


                # invert the subject-to-standard warps, then apply these
                # inverted warps to the segmentation priors, bringing them
                # from standard-space to subject-space
                from anatomical_preproc import apply_ants_transforms_workflow

                wf_config = {}
                wf_config["name"] = "%s_prior_standard_to_subject" % seg
                wf_config["input_image"] = config["%s_prior" % seg]
                wf_config["ref_image"] = resource_pool["anatomical_brain"]
                wf_config["interpolation"] = "NearestNeighbor"
                wf_config["invert"] = True

                wf_config["warps_list"] = [resource_pool["ants_initial_xfm"],\
                                           resource_pool["ants_rigid_xfm"], \
                                           resource_pool["ants_affine_xfm"]]

                # returns:
                # resource_pool["{seg_type}_prior_standard_to_subject_warped"]

                workflow, resource_pool = \
                    apply_ants_transforms_workflow(workflow, resource_pool, \
                                                       wf_config)


    for seg in seg_types:

        # nodes

        overlap_segmentmap_with_prior = pe.Node(interface=fsl.MultiImageMaths(), name='overlap_%s_map_with_prior' % seg)
        overlap_segmentmap_with_prior.inputs.op_string = '-mas %s '

        binarize_threshold_segmentmap = pe.Node(interface=fsl.ImageMaths(), name='binarize_threshold_%s' % seg)

        segment_mask = pe.Node(interface=fsl.MultiImageMaths(), name='%s_mask' % seg)
        segment_mask.inputs.op_string = '-mas %s '


        # overlapping

        if len(resource_pool["%s_probability_map" % seg]) == 2:
            node, out_file = resource_pool["%s_probability_map" % seg]
            workflow.connect(node, out_file, overlap_segmentmap_with_prior, 'in_file')
        else:
            overlap_segmentmap_with_prior.inputs.in_file = resource_pool["%s_probability_map" % seg]


        if len(resource_pool["%s_prior_standard_to_subject_warped" % seg]) == 2:
            node, out_file = resource_pool["%s_prior_standard_to_subject_warped" % seg]
            workflow.connect(node, out_file, overlap_segmentmap_with_prior, 'operand_files')
        else:
            overlap_segmentmap_with_prior.inputs.operand_files = resource_pool["%s_prior_standard_to_subject_warped" % seg]


        # binarize
        workflow.connect(overlap_segmentmap_with_prior, 'out_file', binarize_threshold_segmentmap, 'in_file')
        
        binarize_threshold_segmentmap.inputs.op_string = form_threshold_string(config["%s_threshold" % seg])


        # create segment mask
        workflow.connect(binarize_threshold_segmentmap, 'out_file', segment_mask, 'in_file')

        if len(resource_pool["%s_prior_standard_to_subject_warped" % seg]) == 2:
            node, out_file = resource_pool["%s_prior_standard_to_subject_warped" % seg]
            workflow.connect(node, out_file, segment_mask, 'operand_files')
        else:
            segment_mask.inputs.operand_files = resource_pool["%s_prior_standard_to_subject_warped" % seg]


        resource_pool["anatomical_%s_mask" % seg] = (segment_mask, 'out_file')


    return workflow, resource_pool



def apply_ants_transforms_workflow(workflow, resource_pool, wf_config):

    # input image
    # reference
    # warps
    # invert? yes or no
    # interpolation?

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.fsl as fsl
    import nipype.interfaces.ants as ants
    import nipype.interfaces.utility as util

    from workflow_utils import check_input_resources, \
                               check_config_settings

    apply_name = wf_config["name"]

    input_image = wf_config["input_image"]
    ref_image = wf_config["ref_image"]

    list_of_warps = wf_config["warps_list"]

    invert = wf_config["invert"]
    interpolation = wf_config["interpolation"]

    num_warps = len(list_of_warps)


    # nodes

    collect_linear_transforms = pe.Node(util.Merge(num_warps), \
                                        name='%s_collect_linear_transforms' \
                                        % apply_name)

    apply_ants_xfm = pe.Node(interface=ants.ApplyTransforms(), \
                             name='%s_apply_ants_xfm' % apply_name)
    apply_ants_xfm.inputs.interpolation = interpolation #'NearestNeighbor'

    if invert == True:
        # assuming we're inverting all warps provided
        apply_ants_xfm.inputs.invert_transform_flags = [True, True, True]

    '''
    if invert == True:

        # assuming we're inverting all warps provided
        invert_flags = []
        for warp in list_of_warps:
            invert_flags.append(True)

        apply_ants_xfm.inputs.invert_transform_flags = invert_flags
    '''

    # connections

    if len(input_image) == 2:
        node, out_file = input_image
        workflow.connect(node, out_file, apply_ants_xfm, 'input_image')
    else:
        apply_ants_xfm.inputs.input_image = input_image


    if len(ref_image) == 2:
        node, out_file = ref_image
        workflow.connect(node, out_file, apply_ants_xfm, 'reference_image')
    else:
        apply_ants_xfm.inputs.reference_image = ref_image


    if len(list_of_warps[0]) == 2:
        node, out_file = list_of_warps[0]
        workflow.connect(node, out_file, collect_linear_transforms, 'in1')
    else:
        collect_linear_transforms.inputs.in1 = list_of_warps[0]


    if len(list_of_warps[1]) == 2:
        node, out_file = list_of_warps[1]
        workflow.connect(node, out_file, collect_linear_transforms, 'in2')
    else:
        collect_linear_transforms.inputs.in2 = list_of_warps[1]


    if len(list_of_warps[2]) == 2:
        node, out_file = list_of_warps[2]
        workflow.connect(node, out_file, collect_linear_transforms, 'in3')
    else:
        collect_linear_transforms.inputs.in3 = list_of_warps[2]


    '''
    warpcount = 1

    for warp in list_of_warps:

        if len(warp) == 2:
            node, out_file = warp
            workflow.connect(node, out_file, collect_linear_transforms,
                                 'in%d' % warpcount)
        else:
            datasource = nio.DataGrabber(name = "%s_grabber" \
                                             % warp.split("/")[-1])
            datasource.inputs.base_directory = os.getcwd()
            datasource.inputs.template = warp
            workflow.connect(datasource, 'outfiles', 
                                collect_linear_transforms, 'in%d' % warpcount)
            # NEED A WAY TO HANDLE THIS!
            #collect_linear_transforms.inputs.inX = warp

        warpcount += 1
    '''

    workflow.connect(collect_linear_transforms, 'out',
                         apply_ants_xfm, 'transforms')


    resource_pool["%s_warped" % apply_name] = \
        (apply_ants_xfm, 'output_image')


    return workflow, resource_pool


