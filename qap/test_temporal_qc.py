
import pytest
test_sub_dir = "test_data"


@pytest.mark.quick
def test_calculate_percent_outliers():

    from qap.temporal_qc import calculate_percent_outliers

    dataset = [1,1,2,3,27,34,45,49,54,55,67,294,345,352,356,593,632,675,763,\
        764,825,866,2954,4634,4856,5934,29954]

    percent_outliers, IQR = calculate_percent_outliers(dataset)

    out_tuple = (percent_outliers, IQR)

    assert out_tuple == (0.18518518518518517, 747.5)


@pytest.mark.quick
def test_fd_jenkinson():

    import os
    import numpy as np
    import pkg_resources as p
    
    from qap.temporal_qc import fd_jenkinson

    coord_xfm = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "coordinate_transformation.aff12.1D"))

    ref_meanfd = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                     "meanFD.1D"))                    

    meanfd = fd_jenkinson(coord_xfm, out_array=True)
    
    ref_meanfd_arr = np.genfromtxt(ref_meanfd)        
    
    np.testing.assert_array_equal(ref_meanfd_arr, meanfd)
    

@pytest.mark.quick  
def test_outlier_timepoints_no_mask():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import outlier_timepoints

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))
                                  
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "outlier_timepoints_output_nomask.p"))
                                    
    out_list = outlier_timepoints(func_reorient)

    with open(ref_out, "r") as f:
        ref_list = pickle.load(f)
        
    assert out_list == ref_list


@pytest.mark.quick  
def test_outlier_timepoints_with_mask():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import outlier_timepoints

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))

    func_brain_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "functional_brain_mask.nii.gz"))
                                  
    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "outlier_timepoints_output_withmask.p"))
                                    
    out_list = outlier_timepoints(func_reorient, mask_file=func_brain_mask)

    with open(ref_out, "r") as f:
        ref_list = pickle.load(f)
        
    assert out_list == ref_list    


@pytest.mark.quick
def test_quality_timepoints():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.temporal_qc import quality_timepoints

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))

    ref_out = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "quality_timepoints_output.p"))
                                    
    out_list = quality_timepoints(func_reorient)

    with open(ref_out, "r") as f:
        ref_list = pickle.load(f)  
    
    assert out_list == ref_list


@pytest.mark.quick
def test_global_correlation():

    import os
    import pkg_resources as p

    import numpy.testing as nt

    from qap.temporal_qc import global_correlation

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "func_reorient.nii.gz"))
                                  
    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_brain_mask.nii.gz"))

    gcor = global_correlation(func_reorient, func_mask)

    nt.assert_almost_equal(gcor, 0.13903011798720202, decimal=4)
