

def get_idx(in_files, stop_idx=None, start_idx=None):
    """Validate and return the first and the last volume of the functional 
    timeseries (if selected).

    :type in_files: str
    :param in_files: The filepath to the NIFTI image of the functional
                     timeseries.
    :type start_idx: int
    :param start_idx: (default: None) First volume to be considered, specified
                       by user in the configuration file.
    :type stop_idx: int
    :param stop_idx: (default: None) Last volume to be considered, specified
                     by user in the configuration file.
    :rtype: int
    :return: Value of first volume to consider for the functional run.
    :rtype: int
    :return: Value of last volume to consider for the functional run.
    """

    from nibabel import load

    nvols = load(in_files).shape[3]

    if (start_idx == None) or (start_idx < 0) or (start_idx > (nvols - 1)):
        startidx = 0
    else:
        startidx = start_idx

    if (stop_idx == None) or (stop_idx > (nvols - 1)):
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

    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util
    from nipype.interfaces.afni import preprocess

    if "functional_scan" not in resource_pool.keys():
        return workflow, resource_pool

    if "start_idx" not in config.keys():
        config["start_idx"] = 0

    if "stop_idx" not in config.keys():
        config["stop_idx"] = None

    drop_trs = False
    if (config["start_idx"] != 0) and (config["stop_idx"] != None):
        drop_trs = True

    func_get_idx = pe.Node(util.Function(input_names=['in_files', 
                                                      'stop_idx', 
                                                      'start_idx'],
                                         output_names=['stopidx', 
                                                       'startidx'],
                                         function=get_idx),
                                         name='func_get_idx%s' % name)

    func_get_idx.inputs.in_files = resource_pool["functional_scan"]
    func_get_idx.inputs.start_idx = config["start_idx"]
    func_get_idx.inputs.stop_idx = config["stop_idx"]
    
    if drop_trs:
        func_drop_trs = pe.Node(interface=preprocess.Calc(),
                                name='func_drop_trs%s' % name)

        func_drop_trs.inputs.in_file_a = resource_pool["functional_scan"]
        func_drop_trs.inputs.expr = 'a'
        func_drop_trs.inputs.outputtype = 'NIFTI_GZ'

        workflow.connect(func_get_idx, 'startidx',
                         func_drop_trs, 'start_idx')

        workflow.connect(func_get_idx, 'stopidx',
                         func_drop_trs, 'stop_idx')
    

    func_deoblique = pe.Node(interface=preprocess.Refit(),
                            name='func_deoblique%s' % name)
    func_deoblique.inputs.deoblique = True
    
    if drop_trs:
        workflow.connect(func_drop_trs, 'out_file',
                         func_deoblique, 'in_file')
    else:
        func_deoblique.inputs.in_file = resource_pool["functional_scan"]

    func_reorient = pe.Node(interface=preprocess.Resample(),
                               name='func_reorient%s' % name)
    func_reorient.inputs.orientation = 'RPI'
    func_reorient.inputs.outputtype = 'NIFTI_GZ'

    workflow.connect(func_deoblique, 'out_file',
                    func_reorient, 'in_file')

    resource_pool["func_reorient"] = (func_reorient, 'out_file')

    return workflow, resource_pool


def run_func_preproc(functional_scan, start_idx=None, stop_idx=None, \
                         out_dir=None, run=True):
    """Run the 'func_preproc_workflow' function to execute the modular
    workflow with the provided inputs.

    :type functional_scan: str
    :param functional_scan: Filepath to the raw functional timeseries image in
                             a NIFTI file.
    :type start_idx: int
    :param start_idx: (default: None) The timeseries timepoint/volume to start
                       with - i.e. will only include this timepoint and the
                      ones after it - setting this to None will include all
                      timepoints from the beginning.
    :type stop_idx: int
    :param stop_idx: (default: None) The timeseries timepoint/volume to end
                     with - i.e. will only include this timepoint and the
                     ones before it - setting this to None will include all
                     timepoints up until the end of the timeseries.
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

    output = "func_preproc"

    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["functional_scan"] = functional_scan

    if start_idx:
        config["start_idx"] = start_idx
    if stop_idx:
        config["stop_idx"] = stop_idx
    
    workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_func_motion_correct')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["func_reorient"]

    workflow.connect(node, out_file, ds, 'func_reorient')

    if run:
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, "func_reorient",\
                                         "*"))[0]
        return outpath
        
    else:
        return workflow, workflow.base_dir
    

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

    import copy
    import nipype.pipeline.engine as pe
    from nipype.interfaces.afni import preprocess

    if "func_reorient" not in resource_pool.keys():
        from functional_preproc import func_preproc_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool
    
    # get the first volume of the time series
    get_func_volume = pe.Node(interface=preprocess.Calc(),
                              name='get_func_volume%s' % name)
         
    get_func_volume.inputs.expr = 'a'
    get_func_volume.inputs.single_idx = 0
    get_func_volume.inputs.outputtype = 'NIFTI_GZ'

    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, get_func_volume, 'in_file_a')
    else:
        get_func_volume.inputs.in_file_a = resource_pool["func_reorient"]
        
    # calculate motion parameters
    func_motion_correct = pe.Node(interface=preprocess.Volreg(),
                             name='func_motion_correct%s' % name)

    func_motion_correct.inputs.args = '-Fourier -twopass'
    func_motion_correct.inputs.zpad = 4
    func_motion_correct.inputs.outputtype = 'NIFTI_GZ'
    
    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_motion_correct, 'in_file')
    else:
        func_motion_correct.inputs.in_file = resource_pool["func_reorient"]

    workflow.connect(get_func_volume, 'out_file',
                     func_motion_correct, 'basefile')

    resource_pool["func_motion_correct"] = (func_motion_correct, 'out_file')
    resource_pool["coordinate_transformation"] = \
        (func_motion_correct, 'oned_matrix_save')

    return workflow, resource_pool


def run_func_motion_correct(func_reorient, out_dir=None, run=True):
    """Run the 'func_motion_correct_workflow' function to execute the modular
    workflow with the provided inputs.

    :type func_reorient: str
    :param func_reorient: Filepath to the deobliqued, reoriented functional
                          timeseries.
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

    output = "func_motion_correct"

    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["func_reorient"] = func_reorient

    
    workflow, resource_pool = \
            func_motion_correct_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_func_motion_correct')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["func_motion_correct"]

    workflow.connect(node, out_file, ds, 'func_motion_correct')
    
    
    ds = pe.Node(nio.DataSink(), name='datasink_coordinate_transformation')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["coordinate_transformation"]

    workflow.connect(node, out_file, ds, 'coordinate_transformation')

    if run:
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, "func_motion_correct",\
                                         "*"))[0]
        return outpath      
    else:
        return workflow, workflow.base_dir


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

    import copy
    import nipype.pipeline.engine as pe
    from nipype.interfaces.afni import preprocess

    if "func_reorient" not in resource_pool.keys():

        from functional_preproc import func_preproc_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool
  
    func_get_brain_mask = pe.Node(interface=preprocess.Automask(),
                                  name='func_get_brain_mask%s' % name)
    func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_get_brain_mask, 'in_file')
    else:
        func_get_brain_mask.inputs.in_file = \
            resource_pool["func_reorient"]

    resource_pool["functional_brain_mask"] = (func_get_brain_mask, 'out_file')

    return workflow, resource_pool


def run_functional_brain_mask(func_reorient, out_dir=None, run=True):
    """Run the 'functional_brain_mask_workflow' function to execute the 
    modular workflow with the provided inputs.

    :type func_reorient: str
    :param func_reorient: Filepath to the deobliqued, reoriented functional
                          timeseries.
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

    output = "functional_brain_mask"

    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)

    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["func_reorient"] = func_reorient
    
    workflow, resource_pool = \
            functional_brain_mask_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, "functional_brain" \
                                "_mask", "*"))[0]
        return outpath      
    else:
        return workflow, workflow.base_dir


def invert_functional_brain_mask_workflow(workflow, resource_pool, config,
    name="_"):
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

    import copy
    import nipype.pipeline.engine as pe
    from nipype.interfaces.afni import preprocess

    if "functional_brain_mask" not in resource_pool.keys():

        from functional_preproc import functional_brain_mask_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            functional_brain_mask_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool
  
    # 3dcalc to invert the binary functional brain mask
    invert_mask = pe.Node(interface=preprocess.Calc(), 
                          name='invert_mask%s' % name)

    invert_mask.inputs.expr = "iszero(a)"
    invert_mask.inputs.outputtype = "NIFTI_GZ"

    # functional_brain_mask -> 3dcalc        
    if len(resource_pool["functional_brain_mask"]) == 2:
        node, out_file = resource_pool["functional_brain_mask"]
        workflow.connect(node, out_file, invert_mask, 'in_file_a')
    else:
        invert_mask.inputs.in_file_a = resource_pool["functional_brain_mask"]

    resource_pool["inverted_functional_brain_mask"] = (invert_mask, 'out_file')

    return workflow, resource_pool


def run_invert_functional_brain_mask(functional_brain_mask, out_dir=None,
    run=True):
    """Run the 'invert_functional_brain_mask_workflow' function to execute the
    modular workflow with the provided inputs.

    :type functional_brain_mask: str
    :param functional_brain_mask: Filepath to the binary functional brain
                                  mask.
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

    output = "inverted_functional_brain_mask"

    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["functional_brain_mask"] = functional_brain_mask
    
    workflow, resource_pool = \
            invert_functional_brain_mask_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, "inverted_functional_"\
                                "brain_mask", "*"))[0]
        return outpath      
    else:
        return workflow, workflow.base_dir


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

    import copy
    import nipype.pipeline.engine as pe
    from nipype.interfaces.afni import preprocess

    if "func_reorient" not in resource_pool.keys():

        from functional_preproc import func_preproc_workflow
        old_rp = copy.copy(resource_pool)
        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)
        if resource_pool == old_rp:
            return workflow, resource_pool
   
    func_mean_tstat = pe.Node(interface=preprocess.TStat(),
                           name='func_mean_tstat%s' % name)

    func_mean_tstat.inputs.options = '-mean'
    func_mean_tstat.inputs.outputtype = 'NIFTI_GZ'

    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_mean_tstat, 'in_file')
    else:
        func_mean_tstat.inputs.in_file = \
            resource_pool["func_reorient"]

    resource_pool["mean_functional"] = (func_mean_tstat, 'out_file')

    return workflow, resource_pool

 
def run_mean_functional(func_reorient, out_dir=None, run=True):
    """Run the 'mean_functional_workflow' function to execute the modular
    workflow with the provided inputs.

    - This workflow will NOT remove background noise.

    :type func_reorient: str
    :param func_reorient: Filepath to the deobliqued, reoriented functional
                          timeseries.
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

    output = "mean_functional"
    workflow = pe.Workflow(name='%s_workflow' % output)

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["func_reorient"] = func_reorient

    workflow, resource_pool = \
        mean_functional_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run:
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, "mean_functional",
                                         "*"))[0] 
        return outpath
    else:
        return workflow, workflow.base_dir
