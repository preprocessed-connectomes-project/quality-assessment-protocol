
import pytest
test_sub_dir = "test_data"


@pytest.mark.quick
def test_get_idx_whole_timeseries():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, "End", 0)
    
    
    assert idx_tuple == (123,0)
    
    
@pytest.mark.quick
def test_get_idx_partial_timeseries():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))    
    
    idx_tuple = get_idx(func_scan, 100, 20)
    
    
    assert idx_tuple == (100,20)
        

@pytest.mark.quick
def test_get_idx_partial_timeseries_overshoot():

    import os
    import pkg_resources as p

    from qap.functional_preproc import get_idx
    
    func_scan = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "functional_scan.nii.gz"))     
    
    idx_tuple = get_idx(func_scan, 250, 20)
    
    
    assert idx_tuple == (123,20)
