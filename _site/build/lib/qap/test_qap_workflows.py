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



def test_workflow_qap_anatomical_spatial():

    ''' unit test for the QAP anatomical spatial workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.qap_workflows import run_single_qap_anatomical_spatial
    from qap.workflow_utils import build_test_case


    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))

    qap_head_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "qap_head_mask", \
                                        "mprage_resample_thresh_maths_maths" \
                                        "_maths.nii.gz"))

    anat_csf_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_csf_mask", \
                                        "segment_seg_0.nii.gz"))

    anat_gm_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                       "anat_1", \
                                       "anatomical_gm_mask", \
                                       "segment_seg_1.nii.gz"))

    anat_wm_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                       "anat_1", \
                                       "anatomical_wm_mask", \
                                       "segment_seg_2.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "qap_anatomical_spatial", \
                                    "graph_qap_anatomical_spatial.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "qap_anatomical_spatial", \
                                     "wf_inputs.txt"))


    subject_id = "1019436"
    session_id = "session_1"
    scan_id = "anat_1"
    site_name = "site_1"


    # build the workflow and return it
    wf, base_dir = run_single_qap_anatomical_spatial(anat_reorient, \
                                                     qap_head_mask, \
                                                     anat_csf_mask, \
                                                     anat_gm_mask, \
                                                     anat_wm_mask, \
                                                     subject_id, session_id, \
                                                     scan_id, site_name, \
                                                     False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_csf_mask, \
                           "anatomical_csf_mask_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_gm_mask, \
                           "anatomical_gm_mask_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_reorient, \
                           "anatomical_reorient_here", 1)
    wf_inputs_string = wf_inputs_string.replace(anat_wm_mask, \
                           "anatomical_wm_mask_here", 1)
    wf_inputs_string = wf_inputs_string.replace(qap_head_mask, \
                           "head_mask_path_here", 1)


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def test_workflow_qap_functional_spatial():

    ''' unit test for the QAP functional spatial workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.qap_workflows import run_single_qap_functional_spatial
    from qap.workflow_utils import build_test_case


    mean_func = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "mean_functional", \
                                    "rest_calc_tshift_resample_volreg" \
                                    "_tstat.nii.gz"))

    func_brain_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "rest_1", \
                                          "functional_brain_mask", \
                                          "rest_calc_tshift_resample_volreg" \
                                          "_mask.nii.gz"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "qap_functional_spatial", \
                                    "graph_qap_functional_spatial.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "qap_functional_spatial", \
                                     "wf_inputs.txt"))


    subject_id = "1019436"
    session_id = "session_1"
    scan_id = "rest_1"
    site_name = "site_1"

    ghost_direction = "x"


    # build the workflow and return it
    wf, base_dir = run_single_qap_functional_spatial(mean_func, \
                                                     func_brain_mask, \
                                                     subject_id, session_id, \
                                                     scan_id, site_name, \
                                                     ghost_direction, \
                                                     False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here", 1)
    wf_inputs_string = wf_inputs_string.replace(func_brain_mask, \
                           "func_brain_mask_here", 1)
    wf_inputs_string = wf_inputs_string.replace(mean_func, \
                           "mean_epi_here", 1)


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def test_workflow_qap_functional_temporal():

    ''' unit test for the QAP functional temporal workflow BUILDER '''

    import os
    import commands
    
    import pkg_resources as p

    from qap.qap_workflows import run_single_qap_functional_temporal
    from qap.workflow_utils import build_test_case


    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_volreg" \
                                      ".nii.gz"))

    func_brain_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "rest_1", \
                                          "functional_brain_mask", \
                                          "rest_calc_tshift_resample_volreg" \
                                          "_mask.nii.gz"))

    coord_xfm = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "coordinate_transformation", \
                                    "rest_calc_tshift_resample.aff12.1D"))

    ref_graph = p.resource_filename("qap", os.path.join("test_data", \
                                    "workflow_reference", \
                                    "qap_functional_temporal", \
                                    "graph_qap_functional_temporal.dot"))
                                    
    ref_inputs = p.resource_filename("qap", os.path.join("test_data", \
                                     "workflow_reference", \
                                     "qap_functional_temporal", \
                                     "wf_inputs.txt"))


    subject_id = "1019436"
    session_id = "session_1"
    scan_id = "rest_1"
    site_name = "site_1"


    # build the workflow and return it
    wf, base_dir = run_single_qap_functional_temporal(func_motion, \
                                                      func_brain_mask, \
                                                      subject_id, session_id,\
                                                      scan_id, site_name, \
                                        coordinate_transformation=coord_xfm, \
                                                      run=False)


    # get the workflow inputs of the workflow being tested
    wf_inputs_string = str(wf.inputs).replace("\n","")
    
    wf_inputs_string = wf_inputs_string.replace(base_dir, \
                           "base_directory_here", 1)
    wf_inputs_string = wf_inputs_string.replace(coord_xfm, \
                           "coord_xfm_matrix_here", 1)
    wf_inputs_string = wf_inputs_string.replace(func_brain_mask, \
                           "func_brain_mask_here", 1)
    wf_inputs_string = wf_inputs_string.replace(func_motion, \
                           "func_motion_correct_here", 1)


    flag, err = build_test_case(wf, ref_inputs, ref_graph, wf_inputs_string)

        
    assert flag == 2, err



def run_all_tests_qap_workflows():

    test_workflow_qap_mask()
    test_workflow_qap_anatomical_spatial()
    test_workflow_qap_functional_spatial()
    test_workflow_qap_functional_temporal()