

test_sub_dir = "test_data/1019436/session_1"


def test_select_thresh():

    import os
    import pkg_resources as p
    
    from qap.qap_workflows_utils import select_thresh
    
    input_skull = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                      "anat_1", \
                                      "anatomical_reorient", \
                                      "mprage_resample.nii.gz"))

    thresh_out = select_thresh(input_skull)


    assert thresh_out == 206



def test_slice_head_mask():

    import os
    import pkg_resources as p

    import nibabel as nb
    import numpy as np
    
    from qap.qap_workflows_utils import slice_head_mask


    infile = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                 "anat_1", \
                                 "anatomical_reorient", \
                                 "mprage_resample.nii.gz"))
                                 
    transform = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "qap_headmask_c3d_xfm", \
                                    "qc_fsl_affine_xfm.mat"))
    
    standard = p.resource_filename("qap", os.path.join("test_data", \
                                   "MNI152_T1_2mm.nii.gz"))
                                   

    slice_mask_path = slice_head_mask(infile, transform, standard)

    # get file info
    slice_mask_img = nb.load(slice_mask_path)

    slice_mask_data = slice_mask_img.get_data()

    # test file
    test_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "qap_headmask_slice_head_mask", \
                                    "mprage_resample_slice_mask.nii.gz"))
    
    test_mask_img = nb.load(test_mask)

    test_mask_data = test_mask_img.get_data()
    

    os.system("rm %s" % slice_mask_path)


    # create a vector of True and False values
    bool_vector = slice_mask_data == test_mask_data

    assert bool_vector.all()
    
    
    
def test_qap_anatomical_spatial():

    import os
    import pkg_resources as p
    
    from qap.qap_workflows_utils import qap_anatomical_spatial

    anat = p.resource_filename("qap", os.path.join(test_sub_dir, \
                               "anat_1", \
                               "anatomical_reorient", \
                               "mprage_resample.nii.gz"))
                               
    head_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "anat_1", \
                                    "qap_head_mask", \
                                    "mprage_resample_thresh_maths_maths" \
                                    "_maths.nii.gz"))
                                    
    gm_path = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "anatomical_gm_mask", \
                                  "segment_seg_1.nii.gz"))
                                  
    wm_path = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                  "anat_1", \
                                  "anatomical_wm_mask", \
                                  "segment_seg_2.nii.gz"))
                                  
    csf_path = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                   "anat_1", \
                                   "anatomical_csf_mask", \
                                   "segment_seg_0.nii.gz"))
    
    subject = "1019436"
    session = "session_1"
    scan = "anat_1"

    qc = qap_anatomical_spatial(anat, head_mask, gm_path, wm_path, csf_path, \
                                    subject, session, scan)

    assert (len(qc.keys()) == 27) and (None not in qc.values())



def test_qap_functional_spatial():

    import os
    import pkg_resources as p
    
    from qap.qap_workflows_utils import qap_functional_spatial

    mean_func = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "mean_functional", \
                                    "rest_calc_tshift_resample_volreg_" \
                                    "tstat.nii.gz"))

    func_mask = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_brain_mask", \
                                    "rest_calc_tshift_resample_volreg_" \
                                    "mask.nii.gz"))
                                    
   
    subject = "1019436"
    session = "session_1"
    scan = "rest_1"

    direction = "y"

    qc = qap_functional_spatial(mean_func, func_mask, direction, subject, \
                                    session, scan)


    assert (len(qc.keys()) == 17) and (None not in qc.values())



def test_qap_functional_temporal():

    import os
    import pkg_resources as p
    
    from qap.qap_workflows_utils import qap_functional_temporal

    mask_path = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "functional_brain_mask", \
                                    "rest_calc_tshift_resample_volreg" \
                                    "_mask.nii.gz"))

    func_path = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                    "rest_1", \
                                    "func_motion_correct", \
                                    "rest_calc_tshift_resample_" \
                                    "volreg.nii.gz"))
    
    matrix = p.resource_filename("qap", os.path.join(test_sub_dir, \
                                 "rest_1", \
                                 "coordinate_transformation", \
                                 "rest_calc_tshift_resample.aff12.1D"))
    
    subject = "1019436"
    session = "session_1"
    scan = "rest_1"

    qc = qap_functional_temporal(func_path, mask_path, matrix, subject, \
                                     session, scan)

    assert (len(qc.keys()) == 10) and (None not in qc.values())



def run_all_tests_qap_workflows_utils():

    test_select_thresh()
    test_slice_head_mask()
    test_qap_anatomical_spatial()
    test_qap_functional_spatial()
    test_qap_functional_temporal()
  

