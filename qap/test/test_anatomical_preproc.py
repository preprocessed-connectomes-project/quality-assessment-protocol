

test_sub_dir = "test_data/1019436/session_1"


def test_anatomical_reorient():

    ''' unit test for the anatomical reorient workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_reorient


    anat_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_scan", \
                                    "mprage.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_reorient", \
                                    "graph.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_reorient", \
                                     "wf_inputs.txt"))

    # build the workflow and return it
    wf = run_anatomical_reorient(anat_scan, False)

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "anatomical_reorient_test", \
                                 "graph.dot")
    
    wf.write_graph(dotfilename=out_graph)
    
    # write the workflow inputs of the workflow we are testing
    out_wf_inputs = os.path.join(os.getcwd(), "anatomical_reorient_test", \
                                     "out_wf_inputs.txt")
    
    with open(out_wf_inputs,"wt") as f:
        print >>f, wf.inputs


    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = f.readlines()

    with open(out_graph,"r") as f:
        out_graph_lines = f.readlines()
              
        
    # get the reference workflow inputs
    with open(ref_inputs,"r") as f:
        ref_inputs_lines = f.readlines()
        
    ref_inputs_string = ""
    
    for line in ref_inputs_lines:
        ref_inputs_string = ref_inputs_string + line
        
    base_dir = os.path.join(os.getcwd(), "anatomical_reorient")
        
    ref_inputs_string = ref_inputs_string.replace("IN_FILE_HERE", anat_scan)
    ref_inputs_string = ref_inputs_string.replace("BASE_DIR_HERE", base_dir)
    
    
    # get the workflow inputs of the workflow being tested
    with open(out_wf_inputs,"r") as f:
        out_wf_inputs_lines = f.readlines()
        
    out_wf_inputs_string = ""
    
    for line in out_wf_inputs_lines:
        out_wf_inputs_string = out_wf_inputs_string + line
        

    # clear temporary working files
    try:
        os.system("rm -R %s" % os.path.join(os.getcwd(), \
                      "anatomical_reorient_test"))
    except:
        pass


    # test the case
    flag = 0
    
    if ref_graph_lines == out_graph_lines:
        flag += 1
        
    if ref_inputs_string == out_wf_inputs_string:
        flag += 1

        
    assert flag == 2



def test_run_anatomical_reorient():

    import os
    import commands

    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_reorient


    if "anatomical_reorient" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_scan", \
                                    "mprage.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "anatomical_reorient", \
                                  "mprage_resample.nii.gz"))

    # run the workflow
    output = run_anatomical_reorient(anat_scan)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R anatomical_reorient")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()
    
    
    
def test_anatomical_skullstrip():

    ''' unit test for the anatomical skullstrip workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_skullstrip


    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_brain", \
                                    "graph.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "wf_inputs.txt"))

    # build the workflow and return it
    wf = run_anatomical_skullstrip(anat_reorient, False)

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "anatomical_skullstrip_test", \
                                 "graph.dot")
    
    wf.write_graph(dotfilename=out_graph)
    
    # write the workflow inputs of the workflow we are testing
    out_wf_inputs = os.path.join(os.getcwd(), "anatomical_skullstrip_test", \
                                     "out_wf_inputs.txt")
    
    with open(out_wf_inputs,"wt") as f:
        print >>f, wf.inputs


    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = f.readlines()

    with open(out_graph,"r") as f:
        out_graph_lines = f.readlines()
              
        
    # get the reference workflow inputs
    with open(ref_inputs,"r") as f:
        ref_inputs_lines = f.readlines()
        
    ref_inputs_string = ""
    
    for line in ref_inputs_lines:
        ref_inputs_string = ref_inputs_string + line
        
    base_dir = os.path.join(os.getcwd(), "anatomical_skullstrip")
        
    ref_inputs_string = ref_inputs_string.replace("IN_FILE_HERE", anat_reorient)
    ref_inputs_string = ref_inputs_string.replace("IN_FILE_A_HERE", anat_reorient)
    ref_inputs_string = ref_inputs_string.replace("BASE_DIR_HERE", base_dir)
    
    
    # get the workflow inputs of the workflow being tested
    with open(out_wf_inputs,"r") as f:
        out_wf_inputs_lines = f.readlines()
        
    out_wf_inputs_string = ""
    
    for line in out_wf_inputs_lines:
        out_wf_inputs_string = out_wf_inputs_string + line
        

    # clear temporary working files
    try:
        os.system("rm -R %s" % os.path.join(os.getcwd(), \
                      "anatomical_skullstrip_test"))
    except:
        pass


    # test the case
    flag = 0
    
    if ref_graph_lines == out_graph_lines:
        flag += 1
        
    if ref_inputs_string == out_wf_inputs_string:
        flag += 1
              
        
    assert flag == 2



def test_flirt_anatomical_linear_registration():

    ''' unit test for the anatomical reorient workflow BUILDER '''

    import os
    import pkg_resources as p

    from qap.anatomical_preproc import run_flirt_anatomical_linear_registration

    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    template_brain = p.resource_filename("qap", os.path.join("test_data", \
                                         "MNI152_T1_2mm_brain.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "flirt_linear_warped_image", \
                                    "graph.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "flirt_linear_warped_image", \
                                     "wf_inputs.txt"))

    # build the workflow and return it
    wf = run_flirt_anatomical_linear_registration(anat_brain, template_brain,\
                                                      False)

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "flirt_linear_reg_test", \
                                 "graph.dot")
    
    wf.write_graph(dotfilename=out_graph)
    
    # write the workflow inputs of the workflow we are testing
    out_wf_inputs = os.path.join(os.getcwd(), "flirt_linear_reg_test", \
                                     "out_wf_inputs.txt")
    
    with open(out_wf_inputs,"wt") as f:
        print >>f, wf.inputs


    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = f.readlines()

    with open(out_graph,"r") as f:
        out_graph_lines = f.readlines()
              
        
    # get the reference workflow inputs
    with open(ref_inputs,"r") as f:
        ref_inputs_lines = f.readlines()
        
    ref_inputs_string = ""
    
    for line in ref_inputs_lines:
        ref_inputs_string = ref_inputs_string + line
        
    base_dir = os.path.join(os.getcwd(), \
                                "flirt_anatomical_linear_registration")
        
    ref_inputs_string = ref_inputs_string.replace("IN_FILE_HERE", anat_brain)
    ref_inputs_string = ref_inputs_string.replace("REFERENCE_HERE", template_brain)
    ref_inputs_string = ref_inputs_string.replace("BASE_DIR_HERE", base_dir)
    
    
    # get the workflow inputs of the workflow being tested
    with open(out_wf_inputs,"r") as f:
        out_wf_inputs_lines = f.readlines()
        
    out_wf_inputs_string = ""
    
    for line in out_wf_inputs_lines:
        out_wf_inputs_string = out_wf_inputs_string + line
        

    # clear temporary working files
    try:
        os.system("rm -R %s" % os.path.join(os.getcwd(), \
                      "flirt_linear_reg_test"))
    except:
        pass


    # test the case
    flag = 0
    
    if ref_graph_lines == out_graph_lines:
        flag += 1
        
    if ref_inputs_string == out_wf_inputs_string:
        flag += 1
              
        
    assert flag == 2



def test_run_flirt_anatomical_linear_registration():

    import os
    import pkg_resources as p

    from qap.anatomical_preproc import run_flirt_anatomical_linear_registration

    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    template_brain = p.resource_filename("qap", os.path.join("test_data", \
                                         "MNI152_T1_2mm_brain.nii.gz"))

    ref_wf_input_file = p.resource_filename("qap", os.path.join(test_sub_dir,\
                                            "anat_1", \
                                            "flirt_affine_xfm", \
                                            "wf_input_string.txt"))

    wf = run_flirt_anatomical_linear_registration(anat_brain, \
                                                      template_brain, 0)


    calc_flirt_warp = wf.get_node("calc_flirt_warp")
    
    inputs = []
    ref_inputs = []
    
    inputs.append(calc_flirt_warp.inputs.in_file)
    inputs.append(calc_flirt_warp.inputs.reference)
    inputs.append(calc_flirt_warp.inputs.cost)
    
    ref_inputs.append(anat_brain)
    ref_inputs.append(template_brain)
    ref_inputs.append("corratio")
    
    flag = 0
    
    for test,ref in zip(inputs,ref_inputs):
        if test == ref:
            flag += 1


    assert flag == 3



def test_run_ants_anatomical_linear_registration():

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_ants_anatomical_linear_registration


    if "ants_anatomical_linear_registration" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    template_brain = p.resource_filename("qap", os.path.join("test_data", \
                                         "MNI152_T1_2mm_brain.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "ants_linear_warped_image", \
                                  "transform_Warped.nii.gz"))

    # run the workflow
    output = run_ants_anatomical_linear_registration(anat_brain, \
                                                         template_brain)

    # make the correlation
    ref_out_data = nb.load(ref_out).get_data()  
    output_data = nb.load(output).get_data()
    
    os.system("rm -R ants_anatomical_linear_registration")


    # create a vector of True and False values
    bool_vector = ref_out_data == output_data

    assert bool_vector.all()
    
    
    
def test_segmentation_workflow():

    ''' unit test for the segmentation workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_segmentation_workflow

    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_csf_mask", \
                                    "graph.dot"))

    ref_inputs = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_csf_mask", \
                                     "wf_inputs.txt"))
                                  


    # build the workflow and return it
    wf = run_segmentation_workflow(anat_brain, False)

    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "segmentation_test", \
                                 "graph.dot")
    
    wf.write_graph(dotfilename=out_graph)
    
    # write the workflow inputs of the workflow we are testing
    out_wf_inputs = os.path.join(os.getcwd(), "segmentation_test", \
                                     "out_wf_inputs.txt")
    
    with open(out_wf_inputs,"wt") as f:
        print >>f, wf.inputs


    # load both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = f.readlines()

    with open(out_graph,"r") as f:
        out_graph_lines = f.readlines()
              
        
    # get the reference workflow inputs
    with open(ref_inputs,"r") as f:
        ref_inputs_lines = f.readlines()
        
    ref_inputs_string = ""
    
    for line in ref_inputs_lines:
        ref_inputs_string = ref_inputs_string + line
        
    base_dir = os.path.join(os.getcwd(), "segmentation")
        
    ref_inputs_string = ref_inputs_string.replace("IN_FILES_HERE", anat_brain)
    ref_inputs_string = ref_inputs_string.replace("BASE_DIR_HERE", base_dir)
    
    
    # get the workflow inputs of the workflow being tested
    with open(out_wf_inputs,"r") as f:
        out_wf_inputs_lines = f.readlines()
        
    out_wf_inputs_string = ""
    
    for line in out_wf_inputs_lines:
        out_wf_inputs_string = out_wf_inputs_string + line
        

    # clear temporary working files
    try:
        os.system("rm -R %s" % os.path.join(os.getcwd(), \
                      "segmentation_test"))
    except:
        pass


    # test the case
    flag = 0
    
    if ref_graph_lines == out_graph_lines:
        flag += 1
        
    if ref_inputs_string == out_wf_inputs_string:
        flag += 1
              
        
    assert flag == 2



def test_run_segmentation_workflow():

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_segmentation_workflow


    if "segmentation" in os.listdir(os.getcwd()):

        err = "\n[!] The output folder for this workflow already exists.\n"

        raise Exception(err)


    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    ref_csf_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "anat_1", \
                                      "anatomical_csf_mask", \
                                      "segment_seg_0.nii.gz"))
                                  
    ref_gm_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_gm_mask", \
                                     "segment_seg_1.nii.gz"))
                                  
    ref_wm_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_wm_mask", \
                                     "segment_seg_2.nii.gz"))

    # run the workflow
    output = run_segmentation_workflow(anat_brain, 0.98, 0.7, 0.98)

    ref_list = [ref_csf_out, ref_gm_out, ref_wm_out]

    correlation_count = 0

    # calculate the correlation
    for out, ref in zip(output, ref_list):
    
        # make the correlation
        ref_out_data = nb.load(ref).get_data()  
        output_data = nb.load(out).get_data()
    
        # create a vector of True and False values
        bool_vector = ref_out_data == output_data

        if bool_vector.all():
            correlation_count += 1


    os.system("rm -R segmentation")


    assert correlation_count == 3



def run_unit_tests_anatomical_preproc():

    test_anatomical_reorient()
    test_anatomical_skullstrip()
    test_flirt_anatomical_linear_registration()
    test_segmentation_workflow()



