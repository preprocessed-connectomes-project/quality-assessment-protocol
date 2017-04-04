
import pytest
import unittest


class TestGetMaskedData(unittest.TestCase):

    def setUp(self):
        # init
        import os
        import pkg_resources as p
        import nibabel as nb
        from qap.qap_utils import get_masked_data as gmd
        self.gmd = gmd

        # inputs
        self.func_mean = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "func_reorient.nii.gz"))
        self.func_mask = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "functional_brain_mask.nii.gz"))
        self.masked_func = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "masked_func.nii.gz"))

        ref_masked = nb.load(self.masked_func)
        self.ref_masked_data = ref_masked.get_data()

    def test_masking_files(self):
        import numpy.testing as nt
        masked_data = self.gmd(self.func_mean, self.func_mask, files=True)
        msg = "%s\n%s" % (str(self.ref_masked_data.shape),
                          str(masked_data.shape))
        nt.assert_array_equal(self.ref_masked_data, masked_data, msg)


@pytest.mark.quick
def test_create_anatomical_background_mask():

    import numpy as np
    from qap.qap_utils import create_anatomical_background_mask

    # middle two rows are the "head", first and last rows are background,
    # including "noise" in the fourth row
    anat_data = np.asarray([
      [0, 0, 0, 0], [0, 0.25, 0.5, 0],
      [2, 2, 2, 2], [0, 0.25, 0.5, 0]
    ])

    fg_mask_data = np.asarray([
      [0, 0, 0, 0], [0, 1, 1, 0],
      [1, 1, 1, 1], [0, 0, 0, 0]    
    ])
    
    ref_bg_mask_data = np.asarray([
      [1, 1, 1, 1], [1, 0, 0, 1],
      [0, 0, 0, 0], [1, 1, 1, 1]
    ])

    bg_mask_data = create_anatomical_background_mask(anat_data, fg_mask_data,
        exclude_zeroes=False)

    np.testing.assert_array_equal(ref_bg_mask_data, bg_mask_data)


@pytest.mark.quick
def test_create_anatomical_background_mask_exclude_zeroes():

    import numpy as np
    from qap.qap_utils import create_anatomical_background_mask

    # middle two rows are the "head", first and last rows are background,
    # including "noise" in the fourth row
    anat_data = np.asarray([
      [0, 0, 0, 0], [0, 0.25, 0.5, 0],
      [2, 2, 2, 2], [0, 0.25, 0.5, 0]
    ])

    fg_mask_data = np.asarray([
      [0, 0, 0, 0], [0, 1, 1, 0],
      [1, 1, 1, 1], [0, 0, 0, 0]    
    ])
    
    ref_bg_mask_data = np.asarray([
      [0, 0, 0, 0], [0, 0, 0, 0],
      [0, 0, 0, 0], [0, 1, 1, 0]
    ])

    bg_mask_data = create_anatomical_background_mask(anat_data, fg_mask_data,
        exclude_zeroes=True)

    np.testing.assert_array_equal(ref_bg_mask_data, bg_mask_data)


@pytest.mark.quick
def test_create_anatomical_background_mask_failure():

    from qap.qap_utils import create_anatomical_background_mask

    # send them in as lists instead of NumPy arrays, which should error out
    anat_data = [
      [0, 0, 0, 0], [0, 0.25, 0.5, 0],
      [2, 2, 2, 2], [0, 0.25, 0.5, 0]
    ]

    fg_mask_data = [
      [0, 0, 0, 0], [0, 1, 1, 0],
      [1, 1, 1, 1], [0, 0, 0, 0]    
    ]
    
    ref_bg_mask_data = [
      [0, 0, 0, 0], [0, 0, 0, 0],
      [0, 0, 0, 0], [0, 1, 1, 0]
    ]

    with pytest.raises(Exception) as excinfo:
        bg_mask_data = create_anatomical_background_mask(fg_mask_data, \
        	anat_data)

    assert "must be a NumPy" in str(excinfo.value)
