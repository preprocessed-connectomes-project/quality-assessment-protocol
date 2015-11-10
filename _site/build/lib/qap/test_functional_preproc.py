

test_sub_dir = "test_data/1019436/session_1"


def test_get_idx_whole_timeseries():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, "End", 0)
    
    
    assert idx_tuple == (151,0)
    
    
    
def test_get_idx_partial_timeseries():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, 100, 20)
    
    
    assert idx_tuple == (100,20)
    
    
    
def test_get_idx_partial_timeseries_overshoot():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, 250, 20)
    
    
    assert idx_tuple == (151,20)
    
                                 
                                    
def test_workflow_func_motion_correct_no_slice_time():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct
    from qap.workflow_utils import build_test_case
   
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "func_motion_correct", \
                                    "graph_func_motion_correct.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "func_motion_correct", \
                                     "wf_inputs.txt"))    
                                   
    # build the workflow
    wf, base_dir = run_func_motion_correct(func_scan, 0, "End", False, False)

    
    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILE_A_HERE", 1)
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILES_HERE")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def test_workflow_func_motion_correct_slice_time():

    import os

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct
    from qap.workflow_utils import build_test_case

   
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "func_motion_correct_slice_time", \
                                    "graph_func_motion_correct" \
                                    "_slice_time.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "func_motion_correct_slice_time", \
                                     "wf_inputs.txt"))    
                                   
    # build the workflow
    wf, base_dir = run_func_motion_correct(func_scan, 0, "End", True, False)
         
    
    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILE_A_HERE",\
                                                1)
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILES_HERE")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err
 
    

def test_workflow_functional_brain_mask_3dautomask():

    import os

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask
    from qap.workflow_utils import build_test_case

   
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "functional_brain_mask_3dautomask", \
                                    "graph_functional_brain_mask" \
                                    "_3dautomask.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "functional_brain_mask_3dautomask", \
                                    "wf_inputs.txt"))
                   

    # build the workflow
    wf, base_dir = run_functional_brain_mask(func_motion, use_bet=False, \
                                             run=False)

    
    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_motion, "IN_FILE_HERE")
            

    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err
 
    
    
def test_workflow_functional_brain_mask_BET():

    import os
 
    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask
    from qap.workflow_utils import build_test_case
   
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                                                        
    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "functional_brain_mask_BET", \
                                    "graph_functional_brain_mask_BET.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "functional_brain_mask_BET", \
                                     "wf_inputs.txt"))
                                   
    # build the workflow
    wf, base_dir = run_functional_brain_mask(func_motion, use_bet=True, \
                                             run=False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_motion, "IN_FILE_HERE")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def test_workflow_mean_functional():

    import os

    import pkg_resources as p

    from qap.functional_preproc import run_mean_functional
    from qap.workflow_utils import build_test_case

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "mean_functional", \
                                    "graph_mean_functional.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "mean_functional", \
                                     "wf_inputs.txt"))
                                   
    # build the workflow
    wf, base_dir = run_mean_functional(func_motion, False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_motion, "IN_FILE_HERE")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def run_all_tests_functional_preproc():

    test_get_idx_whole_timeseries()
    test_get_idx_partial_timeseries()
    test_get_idx_partial_timeseries_overshoot()
    test_workflow_func_motion_correct_no_slice_time()
    test_workflow_func_motion_correct_slice_time()
    test_workflow_functional_brain_mask_3dautomask()
    test_workflow_functional_brain_mask_BET()
    test_workflow_mean_functional()   
    
    
    