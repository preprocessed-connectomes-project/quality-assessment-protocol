

test_sub_dir = "test_data/1019436/session_1"


def test_workflow_anatomical_reorient():

    ''' unit test for the anatomical reorient workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_reorient
    from qap.workflow_utils import build_test_case


    anat_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_scan", \
                                    "mprage.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "anatomical_reorient", \
                                    "graph_anatomical_reorient.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "anatomical_reorient", \
                                     "wf_inputs.txt"))


    # build the workflow and return it
    wf, base_dir = run_anatomical_reorient(anat_scan, False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "BASE_DIRECTORY_HERE")
    wf_inputs_string = wf_inputs_string.replace(anat_scan, "IN_FILE_HERE")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err

  
    
def test_workflow_anatomical_skullstrip():

    ''' unit test for the anatomical skullstrip workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_anatomical_skullstrip
    from qap.workflow_utils import build_test_case


    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "anatomical_skullstrip", \
                                    "graph_anatomical_skullstrip.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "anatomical_skullstrip", \
                                     "wf_inputs.txt"))


    # build the workflow and return it
    wf, base_dir = run_anatomical_skullstrip(anat_reorient, False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here")
    wf_inputs_string = wf_inputs_string.replace(anat_reorient, "in_file_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_reorient, "in_file_a_here")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def test_workflow_flirt_anatomical_linear_registration():

    ''' unit test for the anatomical reorient workflow BUILDER '''

    import os
    import pkg_resources as p

    from qap.anatomical_preproc import run_flirt_anatomical_linear_registration
    from qap.workflow_utils import build_test_case

    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    template_brain = p.resource_filename("qap", os.path.join("test_data", \
                                         "MNI152_T1_2mm_brain.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "flirt_anatomical_linear_registration", \
                                    "graph_flirt_anatomical_linear" \
                                    "_registration.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "flirt_anatomical_linear_registration", \
                                     "wf_inputs.txt"))

    # build the workflow and return it
    wf, base_dir = run_flirt_anatomical_linear_registration(anat_brain, \
                                                            template_brain, \
                                                            False)

    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here")
    wf_inputs_string = wf_inputs_string.replace(anat_brain, "in_file_here")
    wf_inputs_string = wf_inputs_string.replace(template_brain, \
                                                    "reference_here")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err

       
    
def test_workflow_segmentation():

    ''' unit test for the segmentation workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.anatomical_preproc import run_segmentation_workflow
    from qap.workflow_utils import build_test_case

    anat_brain = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "anat_1", \
                                     "anatomical_brain", \
                                     "mprage_resample_calc.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "segmentation", \
                                    "graph_segmentation.dot"))

    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "segmentation", \
                                     "wf_inputs.txt"))


    # build the workflow and return it
    wf, base_dir = run_segmentation_workflow(anat_brain, False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here")

    list_input = "['" + anat_brain + "']"

    wf_inputs_string = wf_inputs_string.replace(list_input, "in_files_here")


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def run_all_tests_anatomical_preproc():

    test_workflow_anatomical_reorient()
    test_workflow_anatomical_skullstrip()
    test_workflow_flirt_anatomical_linear_registration()
    test_workflow_segmentation()


