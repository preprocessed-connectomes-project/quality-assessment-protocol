
test_sub_dir = "test_data/1019436/session_1"


def test_fd_jenkinson():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import fd_jenkinson

    coord_xfm = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "coordinate_transformation", \
                                    "rest_calc_tshift_resample.aff12.1D"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "fd_jenkinson", \
                                  "FD_J.1D"))                                    

    meanfd = fd_jenkinson(coord_xfm)
    
    # do da check
    with open(meanfd,"r") as f:
        test_fd_lines = f.readlines()
        
    with open(ref_out,"r") as f:
        ref_fd_lines = f.readlines()
        
        
    os.system("rm FD_J.1D")
        
    
    assert test_fd_lines == ref_fd_lines
    
    

def test_summarize_fd():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import summarize_fd

    coord_xfm = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "coordinate_transformation", \
                                    "rest_calc_tshift_resample.aff12.1D"))
                   
    out_tuple = summarize_fd(coord_xfm)

    assert out_tuple == (0.050015007171052638, 0.0, 0.0)
    
    
    
def test_summarize_fd_threshold_01():   

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import summarize_fd

    coord_xfm = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "coordinate_transformation", \
                                    "rest_calc_tshift_resample.aff12.1D"))
         
    out_tuple = summarize_fd(coord_xfm, threshold=0.1)

    assert out_tuple == (0.050015007171052638, 14.0, 9.15032679738562)
    
    
    
def test_outlier_timepoints():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import outlier_timepoints

    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                  
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_brain_mask", \
                                    "rest_calc_tshift_resample_volreg" \
                                    "_mask.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "outlier_timepoints", \
                                  "outlier_timepoints_ref_out.p"))
                                    
    out_list = outlier_timepoints(func_motion, func_mask)

    with open(ref_out, "r") as f:
        ref_list = pickle.load(f)
        
    
    assert out_list == ref_list
    
    
    
def test_quality_timepoints():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import quality_timepoints

    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "quality_timepoints", \
                                  "quality_timepoints_output.p"))
                                    
    out_list = quality_timepoints(func_motion)

    with open(ref_out, "r") as f:
        ref_list = pickle.load(f)
        
    
    assert out_list == ref_list



def test_quality_timepoints_no_automask():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import quality_timepoints

    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "quality_timepoints", \
                                  "quality_timepoints_nomask_output.p"))
                                    
    out_list = quality_timepoints(func_motion, False)

    with open(ref_out, "r") as f:
        ref_list = pickle.load(f)
        
    
    assert out_list == ref_list



def test_global_correlation():

    import os
    import pkg_resources as p

    from qap.temporal_qc import global_correlation

    func_motion = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "rest_1", \
                                      "func_motion_correct", \
                                      "rest_calc_tshift_resample_" \
                                      "volreg.nii.gz"))
                                  
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_brain_mask", \
                                    "rest_calc_tshift_resample_volreg" \
                                    "_mask.nii.gz"))

    gcor = global_correlation(func_motion, func_mask)


    assert gcor == 0.0090767564485253263



def run_all_tests_temporal_qc():

    test_fd_jenkinson()
    test_summarize_fd()
    test_summarize_fd_threshold_01()
    test_outlier_timepoints()
    test_quality_timepoints()
    test_quality_timepoints_no_automask()
    test_global_correlation()
    
    
    
