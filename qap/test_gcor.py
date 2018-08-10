import pytest
from qap.qap_utils import NumpyEncoder
from qap.test_utils import get_test_dir


test_sub_dir = "test_data"


@pytest.mark.quick
def test_global_correlation():
    import json
    import os
    import pkg_resources as p

    import numpy.testing as nt

    from qap.gcor import calc_global_correlation

    working_directory = get_test_dir('gcor')

    func_reorient = p.resource_filename("qap", os.path.join(test_sub_dir,
                                                            "func_reorient.nii.gz"))

    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir,
                                                        "functional_brain_mask.nii.gz"))

    global_correlation = calc_global_correlation(func_reorient, func_mask, debug=True)

    json_out_filename = os.path.join(working_directory, 'test_out.json')

    with open(json_out_filename, 'w') as ofd:
        json.dump(global_correlation, ofd, indent=2, cls=NumpyEncoder)

    nt.assert_almost_equal(global_correlation['GCOR'], 0.13893567197118509, decimal=4)
