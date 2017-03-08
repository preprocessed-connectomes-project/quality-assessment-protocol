
import pytest
test_sub_dir = "test_data"


@pytest.mark.quick
def test_summary_mask():

    import os
    import pickle
    import pkg_resources as p

    from qap.spatial_qc import summary_mask
    from qap.qap_utils import load_image, load_mask

    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_reorient.nii.gz"))
                                   
    head_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "qap_head_mask.nii.gz"))

    anat_data = load_image(anat_reorient)
    mask_data = load_mask(head_mask, anat_reorient)

    summary_tuple = summary_mask(anat_data, mask_data)

    assert int(summary_tuple[0]) == 230
    assert int(summary_tuple[1]) == 233
    assert int(summary_tuple[2]) == 157221
    

@pytest.mark.quick
def test_check_datatype():

    import numpy as np
    
    from qap.spatial_qc import check_datatype

    sample_array = [[[0,1,2],[3,4,5]],[[6,7,8],[9,10,11]]]

    # put it into NumPy format
    sample_array = np.asarray(sample_array)

    sample_array_int32 = sample_array.astype('int32')
    sample_array_float32 = sample_array.astype('float32')

    output_int_in = check_datatype(sample_array_int32)
    output_float_in = check_datatype(sample_array_float32)

    assert_list = []

    if output_int_in.dtype == "int32":
        assert_list.append(1)
    else:
        assert_list.append(0)

    if output_float_in.dtype == "int32":
        assert_list.append(1)
    else:
        assert_list.append(0)

    assert assert_list == [1,1]
    

@pytest.mark.quick
def test_snr():

    import numpy.testing as nt
    from qap.spatial_qc import snr

    fg_mean = 419.6682196897284
    bg_std = 7.8977607536544143

    snr_out = snr(fg_mean, bg_std)

    nt.assert_almost_equal(snr_out, 53.13762125492869, decimal=4)


@pytest.mark.quick
def test_cnr():

    import numpy.testing as nt
    from qap.spatial_qc import cnr

    mean_gm = 382.08429670720869
    mean_wm = 641.94028605569667
    std_bg = 7.8977607536544143

    cnr_out = cnr(mean_gm, mean_wm, std_bg)

    nt.assert_almost_equal(cnr_out, 32.902489383240514, decimal=4)
    

@pytest.mark.quick
def test_cortical_contrast():

    import numpy.testing as nt
    from qap.spatial_qc import cortical_contrast

    mean_gm = 382.08429670720869
    mean_wm = 641.94028605569667

    cort_out = cortical_contrast(mean_gm, mean_wm)

    nt.assert_almost_equal(cort_out, 0.5075190453873176, decimal=4)


@pytest.mark.quick
def test_fber():

    import os
    import pickle
    import pkg_resources as p
    
    import numpy.testing as nt

    from qap.spatial_qc import fber
    from qap.qap_utils import load_image, load_mask

    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_reorient.nii.gz"))
                                   
    head_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "qap_head_mask.nii.gz"))

    skull_only_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "skull_only_mask.nii.gz"))

    anat_data = load_image(anat_reorient)
    mask_data = load_mask(head_mask, anat_reorient)
    bg_data = 1 - mask_data

    head_data = load_mask(skull_only_mask, anat_reorient)

    fber_out = fber(anat_data, head_data, bg_data)

    nt.assert_almost_equal(fber_out, 341.72165992685609, decimal=4)


@pytest.mark.quick
def test_efc():

    import os
    import pickle
    import pkg_resources as p
    
    import numpy.testing as nt

    from qap.spatial_qc import efc
    from qap.qap_utils import load_image

    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_reorient.nii.gz"))
                                   
    anat_data = load_image(anat_reorient)

    efc_out = efc(anat_data)

    nt.assert_almost_equal(efc_out, 0.36522517588147252, decimal=3)


@pytest.mark.quick
def test_artifacts_no_qi2():

    import os
    import pickle
    import pkg_resources as p

    import numpy.testing as nt

    from qap.spatial_qc import artifacts
    from qap.qap_utils import load_image, load_mask

    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_reorient.nii.gz"))
                                   
    head_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "qap_head_mask.nii.gz"))

    bg_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_bg_mask.nii.gz"))

    anat_data = load_image(anat_reorient)
    mask_data = load_mask(head_mask, anat_reorient)
    bg_data = load_mask(bg_mask, anat_reorient)

    art_out = artifacts(anat_data, mask_data, bg_data, calculate_qi2=False)

    nt.assert_almost_equal(art_out[0], 0.10064793870393487, decimal=4)


@pytest.mark.skip()
@pytest.mark.quick
def test_artifacts_with_qi2():

    # this will fail until the code in 'if calculate_qi2' is updated

    import os
    import pickle
    import pkg_resources as p

    from qap.spatial_qc import artifacts
    from qap.qap_utils import load_image, load_mask

    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))
                                   
    head_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "qap_head_mask", \
                                    "mprage_resample_thresh_maths_maths_" \
                                    "maths.nii.gz"))

    anat_data = load_image(anat_reorient)
    mask_data = load_mask(head_mask, anat_reorient)

    art_out = artifacts(anat_data, mask_data, calculate_qi2=True)

    # not the actual expected output, needs verification
    assert art_out == 0


@pytest.mark.quick
def test_fwhm_out_vox():

    import os
    import pkg_resources as p
    
    import numpy.testing as nt

    from qap.spatial_qc import fwhm

    anat_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_reorient.nii.gz"))
    
    mask_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "qap_head_mask.nii.gz"))
    

    fwhm_out = fwhm(anat_file, mask_file, out_vox=True)

    nt.assert_almost_equal(fwhm_out[0], 3.7207333333333334, decimal=4)
    nt.assert_almost_equal(fwhm_out[1], 3.8991000000000002, decimal=4)
    nt.assert_almost_equal(fwhm_out[2], 4.4016999999999999, decimal=4)
    nt.assert_almost_equal(fwhm_out[3], 3.997033333333333, decimal=4)


@pytest.mark.quick
def test_fwhm_no_out_vox():

    import os
    import pkg_resources as p
    
    import numpy.testing as nt

    from qap.spatial_qc import fwhm

    anat_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_reorient.nii.gz"))
    
    mask_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "qap_head_mask.nii.gz"))

    fwhm_out = fwhm(anat_file, mask_file, out_vox=False)

    nt.assert_almost_equal(fwhm_out[0], 11.1622, decimal=4)
    nt.assert_almost_equal(fwhm_out[1], 11.6973, decimal=4)
    nt.assert_almost_equal(fwhm_out[2], 13.2051, decimal=4)
    nt.assert_almost_equal(fwhm_out[3], 11.991099999999999, decimal=4)


@pytest.mark.quick
def test_ghost_direction():

    import os
    import pkg_resources as p
    
    import numpy.testing as nt

    from qap.spatial_qc import ghost_direction
    from qap.qap_utils import load_image, load_mask

    mean_epi = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                   "mean_functional.nii.gz"))
                                   
    func_brain_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "functional_brain_mask" \
                                          ".nii.gz"))

    mean_epi_data = load_image(mean_epi)
    funcmask_data = load_mask(func_brain_mask, mean_epi)

    gsr_out_x = ghost_direction(mean_epi_data, funcmask_data, "x")
    gsr_out_y = ghost_direction(mean_epi_data, funcmask_data, "y")
    gsr_out_z = ghost_direction(mean_epi_data, funcmask_data, "z")

    gsr_out_all = (gsr_out_x, gsr_out_y, gsr_out_z)

    nt.assert_almost_equal(gsr_out_all[0], -0.018987976014614105, decimal=4)
    nt.assert_almost_equal(gsr_out_all[1], 0.020795321092009544, decimal=4)
    nt.assert_almost_equal(gsr_out_all[2], 0.06708560138940811, decimal=4)  
    
