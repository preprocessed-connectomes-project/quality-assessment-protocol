
import pytest


@pytest.mark.quick
def test_pick_seg_type():

	from qap.anatomical_preproc_utils import pick_seg_type

	prob_maps = [
      "/path/to/segmap_0.nii.gz",
      "/path/to/segmap_1.nii.gz",
      "/path/to/segmap_2.nii.gz"
	]

	file_csf = pick_seg_type(prob_maps, "csf")
	file_gm = pick_seg_type(prob_maps, "gm")
	file_wm = pick_seg_type(prob_maps, "wm")

	assert file_csf == "/path/to/segmap_0.nii.gz"
	assert file_gm == "/path/to/segmap_1.nii.gz"
	assert file_wm == "/path/to/segmap_2.nii.gz"
