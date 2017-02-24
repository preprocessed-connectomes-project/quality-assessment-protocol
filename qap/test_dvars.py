
import pytest
test_sub_dir = "test_data"


@pytest.mark.quick
def test_remove_zero_variance_voxels():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    import numpy as np
    
    from qap.dvars import remove_zero_variance_voxels

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))
                                      
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_brain_mask.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "no_zero_variance_voxels_mask.p"))
                                   
    func_img = nb.load(func_reorient)
    mask_img = nb.load(func_mask)
    
    func_data = func_img.get_data()
    mask_data = mask_img.get_data()
                                    
    out_mask_data = remove_zero_variance_voxels(func_data, mask_data)
                                    
    with open(ref_out, "r") as f:
        ref_mask_data = pickle.load(f)
        
    np.testing.assert_array_equal(ref_mask_data, out_mask_data)


@pytest.mark.quick
def test_load():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    import numpy as np
    
    from qap.dvars import load

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))
                                      
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_brain_mask.nii.gz"))
                                    
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "loaded_func.p"))
                                    
    func_out_data = load(func_reorient, func_mask)

    # match the reference array
    func_out_data = func_out_data[0:10]
   
    with open(ref_out, "r") as f:
        ref_out_data = pickle.load(f)

    np.testing.assert_array_equal(ref_out_data, func_out_data)
    

@pytest.mark.quick
def test_robust_stdev():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    import numpy as np
    
    from qap.dvars import robust_stdev

    func_data_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                         "loaded_func.p"))
                                                                         
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "robust_stdev_output.p"))
                                    
    with open(func_data_file, "r") as f:
        func_data = pickle.load(f)
    
    with open(ref_out, "r") as f:
        ref_mask_data = pickle.load(f)        
    
    func_out_data = robust_stdev(func_data)
        
    np.testing.assert_array_equal(ref_mask_data, func_out_data)


@pytest.mark.quick
def test_ar1():

    import os
    import pickle
    import pkg_resources as p
    
    import nibabel as nb
    import numpy as np
    
    from qap.dvars import ar1

    func_data_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                         "loaded_func.p"))
                                                                         
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "ar1_output.p"))
                                    
    with open(func_data_file, "r") as f:
        func_data = pickle.load(f)
    
    with open(ref_out, "r") as f:
        ref_out_data = pickle.load(f)
        
    func_out_data = ar1(func_data)
        
    np.testing.assert_array_almost_equal(ref_out_data, func_out_data)                           
                                    
