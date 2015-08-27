

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


    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "graph.dot")
    wf.write_graph(dotfilename=out_graph, simple_form=False)
    
    
    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())
          
        
    # get the reference workflow input string
    with open(ref_inputs,"r") as f:
        ref_inputs_string = f.read()

    ref_inputs_string = ref_inputs_string.replace("\n","")
           
    
    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILE_A_HERE", 1)
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILES_HERE")

            

    # clear generated graph.dot
    try:
        os.system("rm %s" % os.path.join(os.getcwd(), "graph.dot"))
    except:
        pass


    # test the case
    flag = 0
    
    err = ""

    if ref_graph_lines == out_graph_lines:
        flag += 1
    else:
        err = err + "\n" + str(ref_graph_lines)
        err = err + "\n\n" + str(out_graph_lines)
        err = err + "\n\nDependency graph TEST FAIL\n\n"
        
    if ref_inputs_string == wf_inputs_string:
        flag += 1
    else:
        err = err + "\n" + ref_inputs_string
        err = err + "\n\n" + wf_inputs_string
        err = err + "\n\nWorkflow inputs TEST FAIL\n\n"

        
    assert flag == 2, err



def test_run_func_motion_correct_no_slice_time():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct


    if "func_motion_correct" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "func_motion_correct", \
                                  "rest_calc_tshift_resample_volreg.nii.gz"))
                                    
    # run the workflow
    output = run_func_motion_correct(func_scan, 0, "End", False)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R func_motion_correct")

    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()



def test_workflow_func_motion_correct_slice_time():

    import os

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct

   
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


    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "graph.dot")
    wf.write_graph(dotfilename=out_graph, simple_form=False)
    
    
    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())
          
        
    # get the reference workflow input string
    with open(ref_inputs,"r") as f:
        ref_inputs_string = f.read()

    ref_inputs_string = ref_inputs_string.replace("\n","")
           
    
    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILE_A_HERE",\
                                                1)
    wf_inputs_string = wf_inputs_string.replace(func_scan, "IN_FILES_HERE")

            

    # clear generated graph.dot
    try:
        os.system("rm %s" % os.path.join(os.getcwd(), "graph.dot"))
    except:
        pass


    # test the case
    flag = 0
    
    err = ""

    if ref_graph_lines == out_graph_lines:
        flag += 1
    else:
        err = err + "\n" + str(ref_graph_lines)
        err = err + "\n\n" + str(out_graph_lines)
        err = err + "\n\nDependency graph TEST FAIL\n\n"
        
    if ref_inputs_string == wf_inputs_string:
        flag += 1
    else:
        err = err + "\n" + ref_inputs_string
        err = err + "\n\n" + wf_inputs_string
        err = err + "\n\nWorkflow inputs TEST FAIL\n\n"

        
    assert flag == 2, err



def test_run_func_motion_correct_slice_time():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct


    if "func_motion_correct" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_scan", \
                                    "rest.nii.gz"))
    
    ''' NEED A SLICE TIME CORRECTED VERSION OF THIS!!!! NOT COMPLETE '''                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "func_motion_correct", \
                                  "rest_calc_tshift_resample_volreg.nii.gz"))
                                    
    # run the workflow
    output = run_func_motion_correct(func_scan, 0, "End", True)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R func_motion_correct")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()
    
    

def test_workflow_functional_brain_mask_3dautomask():

    import os

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask

   
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


    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "graph.dot")
    wf.write_graph(dotfilename=out_graph, simple_form=False)
    
    
    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())
          
        
    # get the reference workflow input string
    with open(ref_inputs,"r") as f:
        ref_inputs_string = f.read()

    ref_inputs_string = ref_inputs_string.replace("\n","")
           
    
    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(func_motion, "IN_FILE_HERE")

            

    # clear generated graph.dot
    try:
        os.system("rm %s" % os.path.join(os.getcwd(), "graph.dot"))
    except:
        pass


    # test the case
    flag = 0
    
    err = ""

    if ref_graph_lines == out_graph_lines:
        flag += 1
    else:
        err = err + "\n" + str(ref_graph_lines)
        err = err + "\n\n" + str(out_graph_lines)
        err = err + "\n\nDependency graph TEST FAIL\n\n"
        
    if ref_inputs_string == wf_inputs_string:
        flag += 1
    else:
        err = err + "\n" + ref_inputs_string
        err = err + "\n\n" + wf_inputs_string
        err = err + "\n\nWorkflow inputs TEST FAIL\n\n"

        
    assert flag == 2, err



def test_run_functional_brain_mask_3dautomask():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask


    if "functional_brain_mask" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "functional_brain_mask", \
                                  "rest_calc_tshift_resample_volreg_" \
                                  "mask.nii.gz"))
                                    
    # run the workflow
    output = run_functional_brain_mask(func_motion, False)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R functional_brain_mask")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()
    
    
    
def test_run_functional_brain_mask_BET():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask


    if "functional_brain_mask" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "functional_brain_mask", \
                                  "rest_calc_tshift_resample_volreg_" \
                                  "mask_BET.nii.gz"))
                                    
    # run the workflow
    output = run_functional_brain_mask(func_motion, True)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R functional_brain_mask")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()
    
    
    
def test_run_mean_functional():

    import os
    import nibabel as nb

    import pkg_resources as p

    from qap.functional_preproc import run_mean_functional


    if "mean_functional" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)

    
    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "mean_functional", \
                                  "rest_calc_tshift_resample_volreg_" \
                                  "tstat.nii.gz"))
                                    
    # run the workflow
    output = run_mean_functional(func_motion)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R mean_functional")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()



def run_all_tests_functional_preproc():

    test_get_idx_whole_timeseries()
    test_get_idx_partial_timeseries()
    test_get_idx_partial_timeseries_overshoot()
    test_workflow_func_motion_correct_no_slice_time()
    test_workflow_func_motion_correct_slice_time()
    test_workflow_functional_brain_mask_3dautomask()
    #test_run_functional_brain_mask_BET()
    #test_run_mean_functional()   
    
    
    