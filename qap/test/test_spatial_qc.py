
test_sub_dir = "test_data/1019436/session_1"


def test_summary_mask():

    import os
    import pickle
    import pkg_resources as p

    from qap.spatial_qc import summary_mask
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

    summary_tuple = summary_mask(anat_data, mask_data)

    assert summary_tuple == (419.6682196897284, 313.74669703973018, 3349975)
    
    
    
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
    
    
    
def test_snr():

    from qap.spatial_qc import snr

    fg_mean = 419.6682196897284
    bg_std = 7.8977607536544143

    snr_out = snr(fg_mean, bg_std)

    assert snr_out == 53.137621254928689



def test_cnr():

    from qap.spatial_qc import cnr

    mean_gm = 382.08429670720869
    mean_wm = 641.94028605569667
    std_bg = 7.8977607536544143

    cnr_out = cnr(mean_gm, mean_wm, std_bg)

    assert cnr_out == 32.902489383240514
    


def test_fber():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.spatial_qc import fber
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

    fber_out = fber(anat_data, mask_data)

    assert fber_out == 672.5508111505



def test_efc():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.spatial_qc import efc
    from qap.qap_utils import load_image

    anat_reorient = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                        "anat_1", \
                                        "anatomical_reorient", \
                                        "mprage_resample.nii.gz"))
                                   
    anat_data = load_image(anat_reorient)

    efc_out = efc(anat_data)

    assert efc_out == 0.430141308439828



def test_artifacts_no_qi2():

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

    art_out = artifacts(anat_data, mask_data, calculate_qi2=False)

    assert art_out == (0.063601601, None)



def test_artifacts():

    ''' this will fail until the code in 'if calculate_qi2' is updated '''

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

    ''' not the actual expected output, needs verification '''
    assert art_out == 0



def test_fwhm_out_vox():

    import os
    import pkg_resources as p
    
    from qap.spatial_qc import fwhm

    anat_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_reorient", \
                                    "mprage_resample.nii.gz"))
    
    mask_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "qap_head_mask", \
                                    "mprage_resample_thresh_maths_" \
                                    "maths_maths.nii.gz"))
    

    fwhm_out = fwhm(anat_file, mask_file, out_vox=True)

    assert fwhm_out == (3.5655491499068379, 3.9207807010902229, \
                        4.0859699705086703, 3.8512161978600825)



def test_fwhm_no_out_vox():

    import os
    import pkg_resources as p
    
    from qap.spatial_qc import fwhm

    anat_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "anatomical_reorient", \
                                    "mprage_resample.nii.gz"))
    
    mask_file = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "qap_head_mask", \
                                    "mprage_resample_thresh_maths_" \
                                    "maths_maths.nii.gz"))

    fwhm_out = fwhm(anat_file, mask_file, out_vox=False)

    assert fwhm_out == (3.56555, 3.9207800000000002, 4.0859500000000004, \
                        3.85121)



def test_ghost_direction_x():

    import os
    import pickle
    import pkg_resources as p
    
    from qap.spatial_qc import ghost_direction
    from qap.qap_utils import load_image, load_mask

    mean_epi = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                   "rest_1", \
                                   "mean_functional", \
                                   "rest_calc_tshift_resample_volreg_" \
                                   "tstat.nii.gz"))
                                   
    func_brain_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "rest_1", \
                                          "functional_brain_mask", \
                                          "rest_calc_tshift_resample_volreg" \
                                          "_mask.nii.gz"))

    mean_epi_data = load_image(mean_epi)
    funcmask_data = load_mask(func_brain_mask, mean_epi)

    gsr_out = ghost_direction(mean_epi_data, funcmask_data, "x")

    assert gsr_out == -0.016198311



def test_ghost_direction_y():

    import os
    import pickle
    import pkg_resources as p

    from qap.spatial_qc import ghost_direction
    from qap.qap_utils import load_image, load_mask

    mean_epi = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                   "rest_1", \
                                   "mean_functional", \
                                   "rest_calc_tshift_resample_volreg_" \
                                   "tstat.nii.gz"))
                                   
    func_brain_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                          "rest_1", \
                                          "functional_brain_mask", \
                                          "rest_calc_tshift_resample_volreg" \
                                          "_mask.nii.gz"))

    mean_epi_data = load_image(mean_epi)
    funcmask_data = load_mask(func_brain_mask, mean_epi)

    gsr_out = ghost_direction(mean_epi_data, funcmask_data, "y")

    assert gsr_out == 0.016582608



def run_all_tests_spatial_qc():

    test_summary_mask()
    test_check_datatype()
    test_snr()
    test_cnr()
    test_fber()
    test_efc()
    test_artifacts_no_qi2()
    #test_artifacts()
    test_fwhm_out_vox()
    test_fwhm_no_out_vox()
    test_ghost_direction_x()
    test_ghost_direction_y()

    
    
