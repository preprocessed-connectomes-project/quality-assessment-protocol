
import pytest

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
      [0, 0, 0, 0], [0, 0, 0, 0],
      [0, 0, 0, 0], [0, 1, 1, 0]
    ])

    bg_mask_data, no_zeroes = \
        create_anatomical_background_mask(fg_mask_data, anat_data)

    np.testing.assert_array_equal(ref_bg_mask_data, bg_mask_data)
    assert no_zeroes == True


@pytest.mark.quick
def test_create_anatomical_background_mask_leave_zeroes():

    import numpy as np
    from qap.qap_utils import create_anatomical_background_mask

    # middle two rows are the "head", first and last rows are background,
    # including "noise" in the fourth row

    # more than 60% of anat image voxels are zero- will leave them in, note
    # ref background mask
    anat_data = np.asarray([
      [0, 0, 0, 0], [0, 0.25, 0.5, 0],
      [0, 2, 2, 0], [0, 0.25, 0.5, 0]
    ])

    fg_mask_data = np.asarray([
      [0, 0, 0, 0], [0, 1, 1, 0],
      [1, 1, 1, 1], [0, 0, 0, 0]    
    ])
    
    ref_bg_mask_data = np.asarray([
      [1, 1, 1, 1], [1, 0, 0, 1],
      [0, 0, 0, 0], [1, 1, 1, 1]
    ])

    bg_mask_data, no_zeroes = \
        create_anatomical_background_mask(fg_mask_data, anat_data)

    np.testing.assert_array_equal(ref_bg_mask_data, bg_mask_data)
    assert no_zeroes == False


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
        bg_mask_data, no_zeroes = \
            create_anatomical_background_mask(fg_mask_data, anat_data)

    assert "must be a NumPy" in str(excinfo.value)
