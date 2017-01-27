base_test_dir = "/tdata/QAP/qc_test"


def anatomical_reorient_workflow(workflow, resource_pool, config, name="_"):
    """Build a Nipype workflow to deoblique and reorient an anatomical scan
    from a NIFTI file.

    - This is a seminal workflow that can only take an input directly from
      disk (i.e. no Nipype workflow connections/pointers, and this is where
      the pipeline will actually begin). For the sake of building the
      pipeine in reverse, if this workflow is called when there is no input
      file available, this function will return the unmodified workflow and
      resource pool directly back.
    - In conjunction with the other workflow-building functions, if this
      function returns the workflow and resource pool unmodified, each
      function up will do the same until it reaches the top level, allowing
      the pipeline builder to continue "searching" for a base-level input
      without crashing at this one.

    Expected Resources in Resource Pool
      - anatomical_scan: The raw anatomical scan in a NIFTI image.

    New Resources Added to Resource Pool
      - anatomical_reorient: The deobliqued, reoriented anatomical scan.

    Workflow Steps
      1. AFNI's 3drefit to deoblique the anatomical scan.
      2. AFNI's 3dresample to reorient the deobliqued anatomical scan to RPI.

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.
    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    import nipype.pipeline.engine as pe
    from nipype.interfaces.afni import preprocess

    if "anatomical_scan" not in resource_pool.keys():
        return workflow, resource_pool

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


def run_anatomical_reorient(anatomical_scan, out_dir=None, run=True):
    """Run the 'anatomical_reorient_workflow' function to execute the modular
    workflow with the provided inputs.

    :type anatomical_scan: str
    :param anatomical_scan: The filepath to the raw anatomical image in a
                            NIFTI file.
    :type out_dir: str
    :param out_dir: (default: None) The output directory to write the results
                    to; if left as None, will write to the current directory.
    :type run: bool
    :param run: (default: True) Will run the workflow; if set to False, will
                connect the Nipype workflow and return the workflow object
                instead.
    :rtype: str
    :return: (if run=True) The filepath of the generated anatomical_reorient
             file.
    :rtype: Nipype workflow object
    :return: (if run=False) The connected Nipype workflow object.
    :rtype: str
    :return: (if run=False) The base directory of the workflow if it were to
             be run.
    """

    import os
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "anatomical_reorient"

    workflow = pe.Workflow(name='anatomical_reorient_workflow')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
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


def anatomical_skullstrip_workflow(workflow, resource_pool, config, name="_"):
    """Build a Nipype workflow to skullstrip an anatomical image using AFNI's
    3dSkullStrip.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - anatomical_reorient: The deobliqued, reoriented anatomical scan.

    New Resources Added to Resource Pool
      - anatomical_brain: The skull-stripped anatomical image (brain only).

    Workflow Steps
      1. AFNI 3dSkullStrip to create a binary mask selecting only the brain.
      2. AFNI 3dcalc to multiply the anatomical image with this mask.

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.
    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    import copy
    import nipype.pipeline.engine as pe

    from nipype.interfaces.afni import preprocess

    if "anatomical_reorient" not in resource_pool.keys():

        from anatomical_preproc import anatomical_reorient_workflow
        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            anatomical_reorient_workflow(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

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


def run_anatomical_skullstrip(anatomical_reorient, out_dir=None, run=True):
    """Run the 'anatomical_skullstrip_workflow' function to execute the 
    modular workflow with the provided inputs.

    :type anatomical_reorient: str
    :param anatomical_reorient: The filepath of the deobliqued, reoriented
                                anatomical image NIFTI file.
    :type out_dir: str
    :param out_dir: (default: None) The output directory to write the results
                    to; if left as None, will write to the current directory.
    :type run: bool
    :param run: (default: True) Will run the workflow; if set to False, will
                connect the Nipype workflow and return the workflow object
                instead.
    :rtype: str
    :return: (if run=True) The filepath of the generated anatomical_reorient
             file.
    :rtype: Nipype workflow object
    :return: (if run=False) The connected Nipype workflow object.
    :rtype: str
    :return: (if run=False) The base directory of the workflow if it were to
             be run.
    """

    import os
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "anatomical_brain"

    workflow = pe.Workflow(name='anatomical_skullstrip_workflow')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
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


def afni_anatomical_linear_registration(workflow, resource_pool, \
    config, name="_"):
    """Build Nipype workflow to calculate the linear registration (participant
    to template) of an anatomical image using AFNI's 3dAllineate.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Settings in Configuration
      - skull_on_registration: (optional- default: True) Whether or not to
                               accept anatomical_reorient or anatomical_brain
                               as the input for registration.
      - template_head_for_anat: (for skull-on registration) The reference
                                template of the whole head.
      - template_brain_for_anat: (for skull-off registration) The reference
                                 template of the brain without skull.

    Expected Resources in Resource Pool
      - anatomical_reorient: The deobliqued, reoriented anatomical scan.
        OR
      - anatomical_brain: The skull-stripped anatomical image (brain only).

    New Resources Added to Resource Pool
      - afni_linear_warped_image: The anatomical image transformed to the
                                  template (using linear warps).
      - allineate_linear_xfm: The text file containing the linear warp matrix
                              produced by AFNI's 3dAllineate.

    Workflow Steps
      1. AFNI's 3dAllineate to calculate the linear registration.

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.
    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    import copy
    import nipype.pipeline.engine as pe
    import nipype.interfaces.afni as afni

    if "skull_on_registration" not in config.keys():
        config["skull_on_registration"] = True

    calc_allineate_warp = pe.Node(interface=afni.Allineate(),
                                    name='calc_3dAllineate_warp%s' % name)
    calc_allineate_warp.inputs.outputtype = "NIFTI_GZ"

    if config["skull_on_registration"]:

        if "anatomical_reorient" not in resource_pool.keys():

            from anatomical_preproc import anatomical_reorient_workflow
            old_rp = copy.copy(resource_pool)
            workflow, new_resource_pool = \
                anatomical_reorient_workflow(workflow, resource_pool,
                                             config, name)

            if resource_pool == old_rp:
                return workflow, resource_pool

        if len(resource_pool["anatomical_reorient"]) == 2:
            node, out_file = resource_pool["anatomical_reorient"]
            workflow.connect(node, out_file, calc_allineate_warp, 'in_file')
        else:
            calc_allineate_warp.inputs.in_file = \
                resource_pool["anatomical_reorient"]

        calc_allineate_warp.inputs.reference = \
            config["template_head_for_anat"]

        calc_allineate_warp.inputs.out_file = "allineate_warped_head.nii.gz"

    else:

        if "anatomical_brain" not in resource_pool.keys():

            from anatomical_preproc import anatomical_skullstrip_workflow
            old_rp = copy.copy(resource_pool)
            workflow, new_resource_pool = \
                anatomical_skullstrip_workflow(workflow, resource_pool, \
                                               config, name)

            if resource_pool == old_rp:
                return workflow, resource_pool

        if len(resource_pool["anatomical_brain"]) == 2:
            node, out_file = resource_pool["anatomical_brain"]
            workflow.connect(node, out_file, calc_allineate_warp, 'in_file')
        else:
            calc_allineate_warp.inputs.in_file = \
                resource_pool["anatomical_brain"]

        calc_allineate_warp.inputs.reference = \
            config["template_brain_for_anat"]

        calc_allineate_warp.inputs.out_file = \
            "allineate_warped_brain.nii.gz"

    calc_allineate_warp.inputs.out_matrix = "3dallineate_warp"

    resource_pool["allineate_linear_xfm"] = \
        (calc_allineate_warp, 'matrix')

    resource_pool["afni_linear_warped_image"] = \
        (calc_allineate_warp, 'out_file')

    return workflow, resource_pool


def run_afni_anatomical_linear_registration(input_image, reference_image,
                                            skull_on=True, out_dir=None,
                                            run=True):
    """Run the 'afni_anatomical_linear_registration' function to execute the 
    modular workflow with the provided inputs.

    :type input_image: str
    :param input_image: Filepath to either the deobliqued, reoriented
                        anatomical image with skull, or the skullstripped
                        brain.
    :type reference_image: str
    :param reference_image: Filepath to the template image to calculate the
                            registration towards.
    :type skull_on: bool
    :param skull_on: (default: True) Whether to run the registration using the
                      whole head or just the brain (use this to mark which
                     type of input_image you are providing).
    :type out_dir: str
    :param out_dir: (default: None) The output directory to write the results
                    to; if left as None, will write to the current directory.
    :type run: bool
    :param run: (default: True) Will run the workflow; if set to False, will
                connect the Nipype workflow and return the workflow object
                instead.
    :rtype: str
    :return: (if run=True) The filepath of the generated anatomical_reorient
             file.
    :rtype: Nipype workflow object
    :return: (if run=False) The connected Nipype workflow object.
    :rtype: str
    :return: (if run=False) The base directory of the workflow if it were to
             be run.
    """

    import os
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "afni_linear_warped_image"

    workflow = pe.Workflow(name='3dallineate_workflow')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    config["skull_on_registration"] = skull_on

    if skull_on:
        resource_pool["anatomical_reorient"] = input_image
        config["template_head_for_anat"] = reference_image
    else:
        resource_pool["anatomical_brain"] = input_image
        config["template_brain_for_anat"] = reference_image

    workflow, resource_pool = \
            afni_anatomical_linear_registration(workflow, resource_pool,
                                                config)

    ds = pe.Node(nio.DataSink(), name='datasink_3dallineate')
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool["afni_linear_warped_image"]
    workflow.connect(node, out_file, ds, 'afni_linear_warped_image')

    node, out_file = resource_pool["allineate_linear_xfm"]
    workflow.connect(node, out_file, ds, 'allineate_linear_xfm')

    if run:
        workflow.run(plugin='MultiProc', plugin_args=
                     {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir,
                                         "afni_linear_warped_image",
                                         "*"))[0]
        return outpath
    else:
        return workflow, workflow.base_dir



def afni_segmentation_workflow(workflow, resource_pool, config, name="_"):
    """Build a Nipype workflow to generate anatomical tissue segmentation maps
    using AFNI's 3dSeg.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      anatomical_brain: The skull-stripped anatomical image (brain only).

    New Resources Added to Resource Pool
      anatomical_csf_mask: The binary mask mapping the CSF voxels.
      anatomical_gm_mask: The binary mask mapping the gray matter voxels.
      anatomical_wm_mask: The binary mask mapping the white matter voxels.

    Workflow Steps
      1. AFNI 3dSeg to run tissue segmentation on the anatomical brain.
      2. AFNI 3dAFNItoNIFTI to convert the AFNI-format 3dSeg output into a
         NIFTI file (as of Oct 2016 3dSeg cannot be configured to write to
         NIFTI).
      3. AFNI 3dcalc to separate the three masks within the output file into
         three separate images.

    :type workflow: Nipype workflow object
    :param workflow: A Nipype workflow object which can already contain other
                     connected nodes; this function will insert the following
                     workflow into this one provided.
    :type resource_pool: dict
    :param resource_pool: A dictionary defining input files and pointers to
                          Nipype node outputs / workflow connections; the keys
                          are the resource names.
    :type config: dict
    :param config: A dictionary defining the configuration settings for the
                   workflow, such as directory paths or toggled options.
    :type name: str
    :param name: (default: "_") A string to append to the end of each node
                 name.
    :rtype: Nipype workflow object
    :return: The Nipype workflow originally provided, but with this function's
              sub-workflow connected into it.
    :rtype: dict
    :return: The resource pool originally provided, but updated (if
             applicable) with the newest outputs and connections.
    """

    import copy
    import nipype.pipeline.engine as pe
    from nipype.interfaces.afni import preprocess

    if "anatomical_brain" not in resource_pool.keys():

        from anatomical_preproc import anatomical_skullstrip_workflow
        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            anatomical_skullstrip_workflow(workflow, resource_pool, config,
                                           name)

        if resource_pool == old_rp:
            return workflow, resource_pool

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


def run_afni_segmentation(anatomical_brain, out_dir=None, run=True):
    """Run the 'afni_segmentation_workflow' function to execute the modular
    workflow with the provided inputs.

    :type anatomical_brain: str
    :param anatomical_brain: Filepath to the skull-stripped brain image in a
                             NIFTI file.
    :type out_dir: str
    :param out_dir: (default: None) The output directory to write the results
                    to; if left as None, will write to the current directory.
    :type run: bool
    :param run: (default: True) Will run the workflow; if set to False, will
                connect the Nipype workflow and return the workflow object
                instead.
    :rtype: str
    :return: (if run=True) The filepath of the generated anatomical_reorient
             file.
    :rtype: Nipype workflow object
    :return: (if run=False) The connected Nipype workflow object.
    :rtype: str
    :return: (if run=False) The base directory of the workflow if it were to
             be run.
    """

    import os
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "anatomical_afni_segmentation_masks"

    workflow = pe.Workflow(name='afni_segmentation_workflow')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
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
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, "anatomical_*_mask", \
                                         "*"))
        return outpath
    else:
        return workflow, workflow.base_dir
