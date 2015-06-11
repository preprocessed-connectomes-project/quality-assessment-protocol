
import os
import sys

base_QAP_dir = "/tdata/QAP/script_revamp/QAP"

base_test_dir = os.path.join(base_QAP_dir, "qc_test")

test_sub_dir = os.path.join(base_test_dir, "pipeline_folder", "2014113")



def qap_spatial(anatomical_reorient, head_mask_path, anatomical_gm_mask, anatomical_wm_mask, anatomical_csf_mask, subject_id, session_id, scan_id, out_vox=True):

    import os
    import sys

    from qap.spatial_qc import summary_mask, snr, cnr, fber, efc, \
                               artifacts, fwhm
    from qap.qap_utils import load_image, load_mask

    # Load the data
    anat_data = load_image(anatomical_reorient)
    fg_mask = load_mask(head_mask_path, anatomical_reorient)
    bg_mask = 1 - fg_mask

    gm_mask = load_mask(anatomical_gm_mask, anatomical_reorient)
    wm_mask = load_mask(anatomical_wm_mask, anatomical_reorient)
    csf_mask = load_mask(anatomical_csf_mask, anatomical_reorient)


    # Initialize QC
    qc              = dict()

    qc['subject'] = subject_id

    qc['session'] = session_id

    qc['scan'] = scan_id

   
    # FBER
    qc['fber'] = fber(anat_data, fg_mask)
    
    # EFC
    qc['efc'] = efc(anat_data)
    
    # Artifact
    qc['qi1'], _    = artifacts(anat_data, fg_mask, calculate_qi2=False)
    
    # Smoothness in voxels
    tmp             = fwhm(anatomical_reorient, head_mask_path, out_vox=out_vox)
    qc['fwhm_x'], qc['fwhm_y'], qc['fwhm_z'], qc['fwhm'] = tmp
    
    
    # Summary Measures
    qc['fg_mean'], qc['fg_std'], qc['fg_size']      = summary_mask(anat_data, fg_mask)
    qc['bg_mean'], qc['bg_std'], qc['bg_size']      = summary_mask(anat_data, bg_mask)
    
    qc['gm_mean'], qc['gm_std'], qc['gm_size']      = (None, None, None)
    qc['wm_mean'], qc['wm_std'], qc['wm_size']      = (None, None, None)
    qc['csf_mean'], qc['csf_std'], qc['csf_size']   = (None, None, None)
    qc['cnr']   = None
    qc['snr']   = None
    
    # More Summary Measures
    qc['gm_mean'], qc['gm_std'], qc['gm_size']      = summary_mask(anat_data, gm_mask)
    qc['wm_mean'], qc['wm_std'], qc['wm_size']      = summary_mask(anat_data, wm_mask)
    qc['csf_mean'], qc['csf_std'], qc['csf_size']   = summary_mask(anat_data, csf_mask)

    # SNR
    qc['snr']       = snr(qc['fg_mean'], qc['bg_std'])

    # CNR
    qc['cnr']   = cnr(qc['gm_mean'], qc['wm_mean'], qc['bg_std'])
    

    return qc



def qap_spatial_epi(mean_epi, func_brain_mask, subject_id, session_id, scan_id, out_vox=True):

    import os
    import sys

    from qap.spatial_qc import summary_mask, snr, fber, efc, fwhm, ghost_all
    from qap.qap_utils import load_image, load_mask

    # Load the data
    anat_data = load_image(mean_epi)
    fg_mask = load_mask(func_brain_mask, mean_epi)
    bg_mask = 1 - fg_mask

    # Initialize QC
    qc              = dict()

    qc['subject'] = subject_id

    qc['session'] = session_id

    qc['scan'] = scan_id

   
    # FBER
    qc['fber'] = fber(anat_data, fg_mask)
    
    # EFC
    qc['efc'] = efc(anat_data)
    
    
    # Smoothness in voxels
    tmp             = fwhm(mean_epi, func_brain_mask, out_vox=out_vox)
    qc['fwhm_x'], qc['fwhm_y'], qc['fwhm_z'], qc['fwhm'] = tmp
    
    # Ghosting
    tmp         = ghost_all(anat_data, fg_mask)
    qc['ghost_x'], qc['ghost_y'], qc['ghost_z'] = tmp

    
    # Summary Measures
    qc['fg_mean'], qc['fg_std'], qc['fg_size']      = summary_mask(anat_data, fg_mask)
    qc['bg_mean'], qc['bg_std'], qc['bg_size']      = summary_mask(anat_data, bg_mask)
    

    qc['snr']   = None
    

    # SNR
    qc['snr']       = snr(qc['fg_mean'], qc['bg_std'])
    

    return qc



def qap_temporal(func_motion_correct, func_brain_mask, coord_xfm_matrix, subject_id, session_id, scan_id, motion_threshold=1.0):

    import sys

    #sys.path.insert(0,"/home/ubuntu/pcp-qap/qap/qclib")

    from qap.temporal_qc import mean_dvars_wrapper, summarize_fd, \
                                mean_outlier_timepoints, \
                                mean_quality_timepoints

    # DVARS
    mean_dvars  = mean_dvars_wrapper(func_motion_correct, func_brain_mask)

    # Mean FD (Jenkinson)
    (mean_fd, num_fd, percent_fd) = summarize_fd(coord_xfm_matrix, threshold=motion_threshold)

    # 3dTout
    mean_outlier= mean_outlier_timepoints(func_motion_correct, func_brain_mask)

    # 3dTqual
    mean_quality= mean_quality_timepoints(func_motion_correct)

    # Compile
    qc = {
        "subject":  subject_id,
        "session":  session_id,
        "scan":     scan_id, 
        "dvars":    mean_dvars, 
        "mean_fd":  mean_fd, 
        'num_fd':   num_fd, 
        'perc_fd':  percent_fd, 
        "outlier":  mean_outlier,
        "quality":  mean_quality
    }
    

    return qc
    
    
    
def append_to_csv(sub_qap_dict, outfile, append):

    import os
    import csv

    fields = sub_qap_dict.keys()

    # put these at the forefront of the list of header items, to make the
    # output CSV's more readable

    fields = sorted(fields)

    if "subject" in fields:
        fields.remove("subject")
        fields.insert(0, "subject")

    if "session" in fields:
        fields.remove("session")
        fields.insert(1, "session")

    if "scan" in fields:
        fields.remove("scan")
        fields.insert(2, "scan")

    if append:

        if not os.path.isfile(outfile):

            with open(outfile, "a") as out_f:

                csv_writer = csv.DictWriter(out_f, fields)

                csv_writer.writeheader()

                csv_writer.writerow(sub_qap_dict)

        else:

            with open(outfile, "a") as out_f:

                csv_writer = csv.DictWriter(out_f, fields)

                csv_writer.writerow(sub_qap_dict)

    else:

            with open(outfile, "wt") as out_f:

                csv_writer = csv.DictWriter(out_f, fields)

                csv_writer.writeheader()

                csv_writer.writerow(sub_qap_dict)


    return outfile



def test_qap_spatial():

    import os

    anat = os.path.join(test_sub_dir, "anatomical_reorient/mprage_resample.nii.gz")
    head_mask = os.path.join(test_sub_dir, "2014113_qc_head_mask.nii.gz")
    gm_path = os.path.join(test_sub_dir, "anatomical_gm_mask/_gm_threshold_0.7/segment_prob_1_maths_maths_maths.nii.gz")
    wm_path = os.path.join(test_sub_dir, "anatomical_wm_mask/_wm_threshold_0.98/segment_prob_2_maths_maths_maths.nii.gz")
    csf_path = os.path.join(test_sub_dir, "anatomical_csf_mask/_csf_threshold_0.98/segment_prob_0_maths_maths_maths.nii.gz")
    subject = "2014113"
    session = "session_1"
    scan = "rest_1"

    qc = qap_spatial(anat, head_mask, gm_path, wm_path, csf_path, subject, session, scan)
    print qc
    assert (len(qc.keys()) == 30) and (None not in qc.values())



def test_qap_spatial_epi():

    import os

    func_mask = os.path.join(test_sub_dir, "functional_brain_mask/_scan_rest_1_rest/rest_calc_tshift_resample_volreg_mask.nii.gz")
    motion_correct = os.path.join(test_sub_dir, "motion_correct/_scan_rest_1_rest/rest_calc_tshift_resample_volreg.nii.gz")
    subject = "2014113"
    session = "session_1"
    scan = "rest_1"

    _,anat_path = mkstemp(suffix=".nii.gz", prefix="tmp_mean_epi_")
    cmd = "fslmaths %s -Tmean %s" % (motion_correct, anat_path)
    os.system(cmd)

    qc = run_qap_spatial_epi(anat_path, func_mask, subject, session, scan)

    assert (len(qc.keys()) == 30) and (None not in qc.values())



def test_qap_temporal():

    import os

    mask_path = os.path.join(test_sub_dir, "functional_brain_mask/_scan_rest_1_rest/rest_calc_tshift_resample_volreg_mask.nii.gz")
    func_path = os.path.join(test_sub_dir, "motion_correct/_scan_rest_1_rest/rest_calc_tshift_resample_volreg.nii.gz")
    matrix = os.path.join(test_sub_dir, "coordinate_transformation/_scan_rest_1_rest/rest_calc_tshift_resample.aff12.1D")
    subject = "2014113"
    scan = "rest_1"

    qc = run_qc_temporal(func_path, mask_path, matrix, subject=subject,
                         scan=scan)

    assert (len(qc.keys()) == 8) and (None not in qc.values())



def run_all_tests():

    test_qap_spatial()
    test_qap_spatial_epi()
    test_qap_temporal()



