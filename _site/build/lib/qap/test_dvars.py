
test_sub_dir = "test_data/1019436/session_1"


def test_remove_zero_variance_voxels():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    
    from qap.dvars import remove_zero_variance_voxels

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
                                  "dvars_data", \
                                  "no_zero_variance_voxels_mask.p"))
                                   
    func_img = nb.load(func_motion)
    mask_img = nb.load(func_mask)
    
    func_data = func_img.get_data()
    mask_data = mask_img.get_data()
                                    
    out_mask_data = remove_zero_variance_voxels(func_data, mask_data)
                                    
    with open(ref_out, "r") as f:
        ref_mask_data = pickle.load(f)
        
    # create a vector of True and False values
    bool_vector = ref_mask_data == out_mask_data

    assert bool_vector.all()



def test_load():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    
    from qap.dvars import load

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
                                  "dvars_data", \
                                  "loaded_func.p"))
                                    
    func_out_data = load(func_motion, func_mask)

    # to match the size of the reference output (shortened for file size
    # issues)
    func_out_data = func_out_data[0:20]
    
    with open(ref_out, "r") as f:
        ref_out_data = pickle.load(f)
        
    # create a vector of True and False values
    bool_vector = ref_out_data == func_out_data
    
    assert bool_vector.all()
    
    
    
def test_robust_stdev():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    
    from qap.dvars import robust_stdev

    func_data_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                         "rest_1", \
                                         "dvars_data", \
                                         "loaded_func.p"))
                                                                         
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "dvars_data", \
                                  "robust_stdev_output.p"))
                                    
    with open(func_data_file, "r") as f:
        func_data = pickle.load(f)
    
    with open(ref_out, "r") as f:
        ref_mask_data = pickle.load(f)
        
    
    func_out_data = robust_stdev(func_data)
    
        
    # create a vector of True and False values
    bool_vector = ref_mask_data == func_out_data
    
    assert bool_vector.all()



def test_ar1():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    
    from qap.dvars import ar1

    func_data_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                         "rest_1", \
                                         "dvars_data", \
                                         "loaded_func.p"))
                                                                         
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "rest_1", \
                                  "dvars_data", \
                                  "ar1_output.p"))
                                    
    with open(func_data_file, "r") as f:
        func_data = pickle.load(f)
    
    with open(ref_out, "r") as f:
        ref_out_data = pickle.load(f)
        
        
    func_out_data = ar1(func_data)
        
        
    # create a vector of True and False values
    bool_vector = ref_out_data == func_out_data
    
    assert bool_vector.all()                           
                                    
                                

def run_all_tests_dvars():

    test_remove_zero_variance_voxels()
    test_load()
    test_robust_stdev()
    test_ar1()