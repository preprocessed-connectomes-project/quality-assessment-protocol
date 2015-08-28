test_sub_dir = "test_data/1019436/session_1"


def test_workflow_qap_mask():

    ''' unit test for the QAP head mask workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.qap_workflows import run_qap_mask
    from qap.workflow_utils import build_test_case


    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))

    flirt_affine_xfm = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                           "anat_1", \
                                           "flirt_affine_xfm", \
                                           "flirt_affine_xfm.mat"))

    template_skull = p.resource_filename("qap", os.path.join("test_data", \
                                         "MNI152_T1_2mm.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "qap_head_mask", \
                                    "graph_qap_head_mask.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "qap_head_mask", \
                                     "wf_inputs.txt"))


    # build the workflow and return it
    wf, base_dir = run_qap_mask(anat_reorient, flirt_affine_xfm, \
                                    template_skull, False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_reorient, \
                           "input_skull_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_reorient, \
                           "infile_here", 1)
    wf_inputs_string = wf_inputs_string.replace(template_skull, \
                           "standard_here", 1)
    wf_inputs_string = wf_inputs_string.replace(flirt_affine_xfm, \
                           "transform_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_reorient, \
                           "in_file_here", 1)


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err


