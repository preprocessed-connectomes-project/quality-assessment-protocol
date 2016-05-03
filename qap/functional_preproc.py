

def get_idx(in_files, stop_idx=None, start_idx=None):

    """
    Method to get the first and the last volume for
    the functional run. It verifies the user specified
    first and last volume. If the values are not valid, it 
    calculates and returns the very first and the last slice 
    
    Parameters
    ----------
    in_file : string (nifti file)
       Path to input functional run
        
    stop_idx : int
        Last volume to be considered, specified by user
        in the configuration file 
    
    stop_idx : int
        First volume to be considered, specified by user 
        in the configuration file 
    
    Returns
    -------
    stop_idx :  int
        Value of first volume to consider for the functional run 
        
    start_idx : int 
        Value of last volume to consider for the functional run
        
    """

    #stopidx = None
    #startidx = None
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

    # resource pool should have:
    #     functional_scan

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.maths as fsl

    from nipype.interfaces.afni import preprocess

    from workflow_utils import check_input_resources, \
                               check_config_settings


    check_input_resources(resource_pool, "functional_scan")

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



def run_func_preproc(functional_scan, start_idx, stop_idx, run=True):

    # stand-alone runner for functional preproc workflow

    import os
    import sys
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='func_preproc_workflow')

    current_dir = os.getcwd()
    workflow_dir = os.path.join(current_dir, "func_preproc")

    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["functional_scan"] = functional_scan

    config["start_idx"] = start_idx
    config["stop_idx"] = stop_idx

    
    workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_func_motion_correct')
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool["func_reorient"]

    workflow.connect(node, out_file, ds, 'func_reorient')


    if run == True:
        try:
            workflow.run(plugin='ResourceMultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "func_reorient", \
                                "*"))[0]
        except:
            workflow.run(plugin='MultiProc', plugin_args= \
                             {'n_procs': num_cores_per_subject})
            outpath = glob.glob(os.path.join(workflow_dir, "func_reorient",\
                                             "*"))[0]

        return outpath
        
    else:
    
        return workflow, workflow.base_dir
    


def func_motion_correct_workflow(workflow, resource_pool, config, name="_"):

    # resource pool should have:
    #     func_reorient

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.maths as fsl

    from nipype.interfaces.afni import preprocess

    from workflow_utils import check_input_resources, \
                               check_config_settings


    if "func_reorient" not in resource_pool.keys():

        from functional_preproc import func_preproc_workflow

        workflow, resource_pool = \
            func_preproc_workflow(workflow, resource_pool, config, name)

    
    func_get_mean_RPI = pe.Node(interface=preprocess.TStat(),
                            name='func_get_mean_RPI%s' % name)
    func_get_mean_RPI.inputs.options = '-mean'
    func_get_mean_RPI.inputs.outputtype = 'NIFTI_GZ'
    
    if len(resource_pool["func_reorient"]) == 2:
        node, out_file = resource_pool["func_reorient"]
        workflow.connect(node, out_file, func_get_mean_RPI, 'in_file')
    else:
        func_get_mean_RPI.inputs.in_file = resource_pool["func_reorient"]


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

    #workflow.connect(func_get_mean_RPI, 'out_file',
    workflow.connect(get_func_volume, 'out_file',
                     func_motion_correct, 'basefile')


    resource_pool["func_motion_correct"] = (func_motion_correct, 'out_file')
    resource_pool["coordinate_transformation"] = \
        (func_motion_correct, 'oned_matrix_save')


    return workflow, resource_pool



def run_func_motion_correct(func_reorient, run=True):

    # stand-alone runner for functional motion correct workflow

    import os
    import sys
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    workflow = pe.Workflow(name='func_motion_correct_workflow')

    current_dir = os.getcwd()
    workflow_dir = os.path.join(current_dir, "func_motion_correct")

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


    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "func_motion_correct",\
                                         "*"))[0]

        return outpath
        
    else:
    
        return workflow, workflow.base_dir



def functional_brain_mask_workflow(workflow, resource_pool, config, name="_"):

    # resource pool should have:
    #     func_motion_correct

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl as fsl
    
    from nipype.interfaces.afni import preprocess


    #check_input_resources(resource_pool, "func_motion_correct")

    if "use_bet" not in config.keys():
        config["use_bet"] = False


    if "func_motion_correct" not in resource_pool.keys():

        from functional_preproc import func_motion_correct_workflow

        workflow, resource_pool = \
            func_motion_correct_workflow(workflow, resource_pool, config, name)

  
    if config["use_bet"] == False:

        func_get_brain_mask = pe.Node(interface=preprocess.Automask(),
                                      name='func_get_brain_mask%s' % name)

        func_get_brain_mask.inputs.outputtype = 'NIFTI_GZ'

    else:

        func_get_brain_mask = pe.Node(interface=fsl.BET(),
                                      name='func_get_brain_mask_BET%s' % name)

        func_get_brain_mask.inputs.mask = True
        func_get_brain_mask.inputs.functional = True

        erode_one_voxel = pe.Node(interface=fsl.ErodeImage(),
                                  name='erode_one_voxel%s' % name)

        erode_one_voxel.inputs.kernel_shape = 'box'
        erode_one_voxel.inputs.kernel_size = 1.0
       

    #if isinstance(tuple, resource_pool["func_motion_correct"]):
        
    if len(resource_pool["func_motion_correct"]) == 2:
        node, out_file = resource_pool["func_motion_correct"]
        workflow.connect(node, out_file, func_get_brain_mask, 'in_file')
    else:
        func_get_brain_mask.inputs.in_file = \
            resource_pool["func_motion_correct"]


    if config["use_bet"] == False:

        resource_pool["functional_brain_mask"] = (func_get_brain_mask, \
                                                     'out_file')

    else:

        workflow.connect(func_get_brain_mask, 'mask_file',
                             erode_one_voxel, 'in_file')

        resource_pool["functional_brain_mask"] = (erode_one_voxel, 'out_file')


    return workflow, resource_pool



def run_functional_brain_mask(func_motion_correct, use_bet=False, run=True):

    # stand-alone runner for functional brain mask workflow

    import os
    import sys
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "functional_brain_mask"

    workflow = pe.Workflow(name='%s_workflow' % output)

    current_dir = os.getcwd()
    workflow_dir = os.path.join(current_dir, output)

    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["func_motion_correct"] = func_motion_correct

    config["use_bet"] = use_bet
    
    workflow, resource_pool = \
            functional_brain_mask_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "functional_brain" \
                                "_mask", "*"))[0]

        return outpath
        
    else:
    
        return workflow, workflow.base_dir



def mean_functional_workflow(workflow, resource_pool, config, name="_"):

    # resource pool should have:
    #     func_motion_correct

    ''' this version does NOT remove background noise '''

    import os
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.maths as fsl
    
    from nipype.interfaces.afni import preprocess

    from workflow_utils import check_input_resources
        

    #check_input_resources(resource_pool, "func_motion_correct")
    #check_input_resources(resource_pool, "functional_brain_mask")


    if "func_motion_correct" not in resource_pool.keys():

        from functional_preproc import func_motion_correct_workflow

        workflow, resource_pool = \
            func_motion_correct_workflow(workflow, resource_pool, config, name)
            
   
    func_mean_skullstrip = pe.Node(interface=preprocess.TStat(),
                           name='func_mean_skullstrip%s' % name)

    func_mean_skullstrip.inputs.options = '-mean'
    func_mean_skullstrip.inputs.outputtype = 'NIFTI_GZ'


    if len(resource_pool["func_motion_correct"]) == 2:
        node, out_file = resource_pool["func_motion_correct"]
        workflow.connect(node, out_file, func_mean_skullstrip, 'in_file')#func_edge_detect, 'in_file_a')
    else:
        func_mean_skullstrip.inputs.in_file = \
            resource_pool["func_motion_correct"]
 

    resource_pool["mean_functional"] = (func_mean_skullstrip, 'out_file')


    return workflow, resource_pool


 
def run_mean_functional(func_motion_correct, run=True):

    # stand-alone runner for mean functional workflow

    ''' this version does NOT remove background noise '''

    import os
    import sys
    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "mean_functional"

    workflow = pe.Workflow(name='%s_workflow' % output)

    current_dir = os.getcwd()
    workflow_dir = os.path.join(current_dir, output)

    workflow.base_dir = workflow_dir


    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["func_motion_correct"] = func_motion_correct

    
    workflow, resource_pool = \
            mean_functional_workflow(workflow, resource_pool, config)


    ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
    ds.inputs.base_directory = workflow_dir
    
    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run == True:

        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})

        outpath = glob.glob(os.path.join(workflow_dir, "mean_functional", \
                                         "*"))[0] 

        return outpath
        
    else:
    
        return workflow, workflow.base_dir
        

