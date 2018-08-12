import nipype.pipeline.engine as pe
import nipype.interfaces.utility as util
import nipype.interfaces.afni as afni
import qap.cloud_utils as cloud_utils

if hasattr(afni.preprocess, 'Calc'):
    afni_utils = afni.preprocess
else:
    afni_utils = afni.utils
afni_preprocess = afni.preprocess


def get_idx(in_files, stop_idx=None, start_idx=None):
    """Validate and return the first and the last volume of the functional 
    timeseries (if selected).

    Runs in a nipype node, so it must contain all of its includes.

    :type in_files: str
    :param in_files: The filepath to the NIFTI image of the functional
                     timeseries.
    :type start_idx: int
    :param start_idx: (default: None) First volume to be considered, specified
                       by user in the configuration file.
    :type stop_idx: int or string
    :param stop_idx: (default: None) Last volume to be considered, specified
                     by user in the configuration file. If 'End' or None the last
                     volume of the file (nvols - 1) will be used.
    :rtype: int
    :return: Value of first volume to consider for the functional run.
    :rtype: int
    :return: Value of last volume to consider for the functional run.
    """
    import nibabel

    if isinstance(stop_idx, str):
        if stop_idx.lower() == "end":
            stop_idx = None
        else:
            stop_idx = int(stop_idx)

    nvols = nibabel.load(in_files).shape[3]

    if (start_idx is None) or (start_idx < 0) or (start_idx > (nvols - 1)):
        startidx = 0
    else:
        startidx = start_idx

    if (stop_idx is None) or (stop_idx > (nvols - 1)):
        stopidx = nvols - 1
    else:
        stopidx = stop_idx

    return stopidx, startidx


def func_preproc_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to deoblique and reorient a functional
    scan from a NIFTI file.

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
      - functional_scan: The raw functional 4D timeseries in a NIFTI file.

    New Resources Added to Resource Pool
      - func_reorient: The deobliqued, reoriented functional timeseries.

    Workflow Steps
      1. get_idx function node (if a start_idx and/or stop_idx is set in the
         configuration) to generate the volume range to keep in the timeseries
      2. AFNI 3dcalc to drop volumes not included in the range (if a start_idx
         and/or stop_idx has been set in the configuration only)
      3. AFNI 3drefit to deoblique the file
      4. AFNI 3dresample to reorient the file to RPI

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

    # if this workflows output is already in the resource pool, we have nothing to do, so move on
    if "func_reorient" in resource_pool:
        return workflow, resource_pool

    # since most of the necessary files can be made form the functional_scan, raising a failure here
    # should trickle down through the entire workflows.
    if "functional_scan" not in resource_pool.keys():
        raise ValueError("Error! No functional scan was found in the resource pool and do not know how to make one.")
    elif "s3://" in resource_pool["functional_scan"]:
        resource_pool["functional_scan"] = \
            cloud_utils.download_single_s3_path(resource_pool["functional_scan"], config)

    if "start_idx" not in config.keys():
        config["start_idx"] = 0

    if "stop_idx" not in config.keys():
        config["stop_idx"] = None

    drop_trs = False
    if (config["start_idx"] != 0) and (config["stop_idx"] is not None):
        drop_trs = True

    func_get_idx = pe.Node(util.Function(input_names=['in_files',
                                                      'stop_idx',
                                                      'start_idx'],
                                         output_names=['stopidx',
                                                       'startidx'],
                                         function=get_idx),
                           name='func_get_idx{0}'.format(name))

    func_get_idx.inputs.in_files = resource_pool["functional_scan"]
    func_get_idx.inputs.start_idx = config["start_idx"]
    func_get_idx.inputs.stop_idx = config["stop_idx"]

    func_deoblique = pe.Node(interface=afni_utils.Refit(),
                             name='func_deoblique{0}'.format(name))

    func_deoblique.inputs.deoblique = True

    if drop_trs:
        func_drop_trs = pe.Node(interface=afni_utils.Calc(),
                                name='func_drop_trs{0}'.format(name))

        func_drop_trs.inputs.in_file_a = resource_pool["functional_scan"]
        func_drop_trs.inputs.expr = 'a'
        func_drop_trs.inputs.outputtype = 'NIFTI_GZ'

        workflow.connect(func_get_idx, 'startidx',
                         func_drop_trs, 'start_idx')

        workflow.connect(func_get_idx, 'stopidx',
                         func_drop_trs, 'stop_idx')

        workflow.connect(func_drop_trs, 'out_file',
                         func_deoblique, 'in_file')
    else:
        func_deoblique.inputs.in_file = resource_pool["functional_scan"]

    func_reorient = pe.Node(interface=afni_utils.Resample(),
                            name='func_reorient{0}'.format(name))

    func_reorient.inputs.orientation = 'RPI'
    func_reorient.inputs.outputtype = 'NIFTI_GZ'

    workflow.connect(func_deoblique, 'out_file',
                     func_reorient, 'in_file')

    resource_pool["func_reorient"] = (func_reorient, 'out_file')

    return workflow, resource_pool


def func_motion_correct_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to calculate the motion correction
    parameters of a functional timeseries using AFNI's 3dvolreg.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - func_reorient: The deobliqued, reoriented functional timeseries.

    New Resources Added to Resource Pool:
      - func_motion_correct: The motion-corrected functional timeseries.
      - coordinate_transformation: The matrix transformation from AFNI's
                                   3dvolreg (--1Dmatrix_save option).

    Workflow Steps
      1. AFNI 3dcalc to extract the first volume of the functional timeseries
         for the basefile for 3dvolreg.
      2. AFNI 3dvolreg to calculate the motion correction parameters.

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

    # if this workflows output is already in the resource pool, we have nothing to do, so move on
    if "func_motion_correct" in resource_pool and "func_coordinate_transformation" in resource_pool and \
            "func_motion_estimates" in resource_pool:
        return workflow, resource_pool

    if "func_reorient" not in resource_pool.keys():
        # this should raise an exception that will kill the workflow builder
        # if func_preproc_workflow could not be installed
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)
    
    # get the first volume of the time series
    get_func_volume = pe.Node(interface=afni_utils.Calc(),
                              name='get_func_volume{0}'.format(name))

    get_func_volume.inputs.expr = 'a'
    get_func_volume.inputs.single_idx = 0
    get_func_volume.inputs.outputtype = 'NIFTI_GZ'

    if isinstance(resource_pool["func_reorient"], tuple):
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, get_func_volume, 'in_file_a')
    else:
        get_func_volume.inputs.in_file_a = resource_pool["func_reorient"]
        
    # calculate motion parameters
    func_motion_correct = pe.Node(interface=afni_preprocess.Volreg(),
                                  name='func_motion_correct{0}'.format(name))

    func_motion_correct.inputs.args = '-Fourier -twopass'
    func_motion_correct.inputs.zpad = 4
    func_motion_correct.inputs.outputtype = 'NIFTI_GZ'
    
    if isinstance(resource_pool["func_reorient"], tuple):
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_motion_correct, 'in_file')
    else:
        func_motion_correct.inputs.in_file = resource_pool["func_reorient"]

    workflow.connect(get_func_volume, 'out_file',
                     func_motion_correct, 'basefile')

    resource_pool["func_motion_correct"] = (func_motion_correct, 'out_file')
    resource_pool["func_coordinate_transformation"] = \
        (func_motion_correct, 'oned_matrix_save')
    resource_pool["func_motion_estimates"] = \
        (func_motion_correct, 'oned_file')

    return workflow, resource_pool


def functional_brain_mask_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to generate a functional brain mask
    using AFNI's 3dAutomask.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - func_reorient: The deobliqued, reoriented functional timeseries.

    New Resources Added to Resource Pool
      - functional_brain_mask: The binary brain mask of the functional time
                               series.

    Workflow Steps
      1. AFNI's 3dAutomask to generate the mask.

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
    # if this workflows output is already in the resource pool, we have nothing to do, so move on
    if "func_brain_mask" in resource_pool:
        return workflow, resource_pool

    if "func_reorient" not in resource_pool.keys():
        # this should raise an exception that will kill the workflow builder
        # if func_preproc_workflow could not be installed
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)
  
    func_get_brain_mask = pe.Node(interface=afni_preprocess.Automask(),
                                  name='func_get_brain_mask{0}'.format(name))
    func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    if isinstance(resource_pool["func_reorient"], tuple):
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_get_brain_mask, 'in_file')
    else:
        func_get_brain_mask.inputs.in_file = \
            resource_pool["func_reorient"]

    resource_pool["func_brain_mask"] = (func_get_brain_mask, 'out_file')

    return workflow, resource_pool


def invert_functional_brain_mask_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to generate a background mask of a
    functional scan (the inversion of the functional brain mask) using AFNI's
    3dCalc.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool
      - functional_brain_mask: The binary brain mask of the functional time
                               series.

    New Resources Added to Resource Pool
      - inverted_functional_brain_mask: The inversion of the functional brain
                                        mask, a binary brain mask of the
                                        background of the functional time
                                        series.

    Workflow Steps:
      1. AFNI's 3dcalc to invert the functional brain mask

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
    # if this workflows output is already in the resource pool, we have nothing to do, so move on
    if "func_inverted_brain_mask" in resource_pool:
        return workflow, resource_pool

    if "func_brain_mask" not in resource_pool.keys():
        # this should raise an exception that will kill the workflow builder
        # if func_preproc_workflow could not be installed
        workflow, resource_pool = \
            functional_brain_mask_workflow(workflow, resource_pool, config, name)

    # 3dcalc to invert the binary functional brain mask
    invert_mask = pe.Node(interface=afni_utils.Calc(),
                          name='invert_mask{0}'.format(name))

    invert_mask.inputs.expr = "iszero(a)"
    invert_mask.inputs.outputtype = "NIFTI_GZ"

    # functional_brain_mask -> 3dcalc        
    if isinstance(resource_pool["func_brain_mask"], tuple):
        node, out_file = resource_pool["func_brain_mask"]
        workflow.connect(node, out_file, invert_mask, 'in_file_a')
    else:
        invert_mask.inputs.in_file_a = resource_pool["func_brain_mask"]

    resource_pool["func_inverted_brain_mask"] = (invert_mask, 'out_file')
    return workflow, resource_pool


def mean_functional_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to generate a one-volume image from a
    functional timeseries comprising of the mean of its timepoint values,
    using AFNI's 3dTstat.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.
    - For QAP: This workflow will NOT remove background noise from the
      image, to maintain as accurate of a quality metric as possible.

    Expected Resources in Resource Pool
      - func_reorient: The deobliqued, reoriented functional timeseries.

    New Resources Added to Resource Pool
      - mean_functional: The one-volume image of the averaged timeseries.

    Workflow Steps:
      1. AFNI 3dTstat to calculate the mean of the functional timeseries.

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
    # if this workflows output is already in the resource pool, we have nothing to do, so move on
    if "func_mean" in resource_pool:
        return workflow, resource_pool

    if "func_reorient" not in resource_pool.keys():
        # this should raise an exception that will kill the workflow builder
        # if func_preproc_workflow could not be installed
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)

    func_mean_tstat = pe.Node(interface=afni_utils.TStat(),
                              name='func_mean_tstat{0}'.format(name))

    func_mean_tstat.inputs.options = '-mean'
    func_mean_tstat.inputs.outputtype = 'NIFTI_GZ'

    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_mean_tstat, 'in_file')
    else:
        func_mean_tstat.inputs.in_file = \
            resource_pool["func_reorient"]

    resource_pool["func_mean"] = (func_mean_tstat, 'out_file')

    return workflow, resource_pool


def tstd_functional_workflow(workflow, resource_pool, config, name="_"):
    """Build and run a Nipype workflow to generate a one-volume image from a
    functional timeseries comprising of the tstd of its time point values,
    using AFNI's 3dTstat.

    - If any resources/outputs required by this workflow are not in the
      resource pool, this workflow will call pre-requisite workflow builder
      functions to further populate the pipeline with workflows which will
      calculate/generate these necessary pre-requisites.
    - For QAP: This workflow will NOT remove background noise from the
      image, to maintain as accurate of a quality metric as possible.

    Expected Resources in Resource Pool
      - func_reorient: The deobliqued, reoriented functional timeseries.

    New Resources Added to Resource Pool
      - mean_functional: The one-volume image of the averaged timeseries.

    Workflow Steps:
      1. AFNI 3dTstat to calculate the tstd of the functional timeseries.

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
    # if this workflows output is already in the resource pool, we have nothing to do, so move on
    if "func_tstd" in resource_pool:
        return workflow, resource_pool

    if "func_reorient" not in resource_pool.keys():
        # this should raise an exception that will kill the workflow builder
        # if func_preproc_workflow could not be installed
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)

    func_tstd_tstat = pe.Node(interface=afni_utils.TStat(),
                              name='func_tstd_tstat{0}'.format(name))

    func_tstd_tstat.inputs.options = '-stdev'
    func_tstd_tstat.inputs.outputtype = 'NIFTI_GZ'

    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_tstd_tstat, 'in_file')
    else:
        func_tstd_tstat.inputs.in_file = \
            resource_pool["func_reorient"]

    resource_pool["func_tstd"] = (func_tstd_tstat, 'out_file')

    return workflow, resource_pool
