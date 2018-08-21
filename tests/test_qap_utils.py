
import pytest
import unittest


class TestGetMaskedData(unittest.TestCase):

    def setUp(self):
        # init
        import os
        import pkg_resources as p
        import nibabel as nb
        import numpy.testing as nt
        from qap.qap_utils import get_masked_data as gmd
        self.gmd = gmd
        self.nt = nt
        self.nb = nb

        # inputs
        self.func_reorient = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "func_reorient.nii.gz"))
        self.func_mean = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "mean_functional.nii.gz"))
        self.func_mask = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "functional_brain_mask.nii.gz"))
        self.masked_func = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "masked_func.nii.gz"))
        self.masked_mean_func = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "masked_mean_func.nii.gz"))

        ref_masked = nb.load(self.masked_func)
        self.ref_masked_ts = ref_masked.get_data()

        ref_masked = nb.load(self.masked_mean_func)
        self.ref_masked_vol = ref_masked.get_data()

    def test_masking_timeseries(self):
        # inputs are filepath strings, for a 4D input file
        masked_data = self.gmd(self.func_reorient, self.func_mask)
        msg = "%s\n%s" % (str(self.ref_masked_ts.shape),
                          str(masked_data.shape))
        self.nt.assert_array_equal(self.ref_masked_ts, masked_data, msg)

    def test_masking_one_volume(self):
        # inputs are filepath strings, for a 3D input file
        masked_data = self.gmd(self.func_mean, self.func_mask)
        msg = "%s\n%s" % (str(self.ref_masked_vol.shape),
                          str(masked_data.shape))
        self.nt.assert_array_equal(self.ref_masked_vol, masked_data, msg)

    def test_masking_timeseries_data(self):
        # inputs are Numpy data arrays, for a 4D input file
        func_data = self.nb.load(self.func_reorient).get_data()
        mask_data = self.nb.load(self.func_mask).get_data()
        masked_data = self.gmd(func_data, mask_data)
        msg = "%s\n%s" % (str(self.ref_masked_vol.shape),
                          str(masked_data.shape))
        self.nt.assert_array_equal(self.ref_masked_ts, masked_data, msg)

    def test_masking_one_volume_data(self):
        # inputs are Numpy data arrays, for a 3D input file
        func_data = self.nb.load(self.func_mean).get_data()
        mask_data = self.nb.load(self.func_mask).get_data()
        masked_data = self.gmd(func_data, mask_data)
        msg = "%s\n%s" % (str(self.ref_masked_vol.shape),
                          str(masked_data.shape))
        self.nt.assert_array_equal(self.ref_masked_vol, masked_data, msg)


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


@pytest.mark.quick
def test_create_bundles():

    from qap.qap_utils import create_bundles
    
    scans = {
        'task-read': {},
        'task-listen': {},
        'rest': {},
    }

    sessions = {
        'ses-0001': scans,
        'ses-0002': scans,
    }

    subjects = {
        'sub-0001': sessions,
        'sub-0002': sessions,
        'sub-0003': sessions,
        'sub-0004': sessions,
    }

    sites = {
        'NYU': subjects,
        'CALTECH': subjects,
    }

    bundles = create_bundles(sites, 5)

    keys = {
        tranche_id: i
        for i, bundle in enumerate(bundles)
        for tranche_id, _ in bundle.items()
    }

    # Scans must be in the same bundle
    for (site_1, sub_1, ses_1, _), bundle_1 in keys.items():
        for (site_2, sub_2, ses_2, _), bundle_2 in keys.items():
            if site_1 == site_2 and sub_1 == sub_2 and ses_1 == ses_2:
                assert bundle_1 == bundle_2
