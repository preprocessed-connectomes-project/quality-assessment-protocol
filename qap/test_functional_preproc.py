

test_sub_dir = "test_data"


def test_get_idx_whole_timeseries():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, "End", 0)
    
    
    assert idx_tuple == (123,0)
    
    
    
def test_get_idx_partial_timeseries():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, 100, 20)
    
    
    assert idx_tuple == (100,20)
    
    
    
def test_get_idx_partial_timeseries_overshoot():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))     
    
    idx_tuple = get_idx(func_scan, 250, 20)
    
    
    assert idx_tuple == (123,20)



def test_run_func_preproc():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.functional_preproc import run_func_preproc

    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "func_reorient.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_func_preproc")
    output = run_func_preproc(func_scan, out_dir=out_dir)

    ref_out_sform = nb.load(ref_output).get_header().get_sform()
    out_sform = nb.load(output).get_header().get_sform()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_sform, out_sform)



def test_run_func_motion_correct():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.functional_preproc import run_func_motion_correct
   
    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "func_motion_correct.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_func_preproc")
    output = run_func_motion_correct(func_reorient, out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def test_run_functional_brain_mask_3dAutoMask():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask
   
    func_motion_correct = p.resource_filename("qap", \
                                              os.path.join(test_sub_dir, \
                                              "func_motion_correct.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "functional_brain_mask_3dAutoMask" \
                                     ".nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_func_preproc")
    output = run_functional_brain_mask(func_motion_correct, out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def test_run_functional_brain_mask_BET():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.functional_preproc import run_functional_brain_mask
   
    func_motion_correct = p.resource_filename("qap", \
                                              os.path.join(test_sub_dir, \
                                              "func_motion_correct.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "functional_brain_mask_BET.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_func_preproc")
    output = run_functional_brain_mask(func_motion_correct, use_bet=True, \
                                           out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)
 


def test_run_mean_functional():

    import os
    import numpy as np
    import nibabel as nb
    import shutil

    import pkg_resources as p

    from qap.functional_preproc import run_mean_functional
   
    func_motion_correct = p.resource_filename("qap", \
                                              os.path.join(test_sub_dir, \
                                              "func_motion_correct.nii.gz"))

    ref_output = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "mean_functional.nii.gz"))

    out_dir = os.path.join(os.getcwd(), "unit_tests_func_preproc")
    output = run_mean_functional(func_motion_correct, out_dir=out_dir)

    ref_out_data = nb.load(ref_output).get_data()
    out_data = nb.load(output).get_data()

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    np.testing.assert_array_equal(ref_out_data, out_data)



def run_all_tests_functional_preproc():

    test_get_idx_whole_timeseries()
    test_get_idx_partial_timeseries()
    test_get_idx_partial_timeseries_overshoot()
    test_run_func_preproc()
    test_run_func_motion_correct()
    test_run_functional_brain_mask_3dAutoMask()
    test_run_functional_brain_mask_BET()
    test_run_mean_functional()   
    
    
    