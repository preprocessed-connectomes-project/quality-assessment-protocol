"""
Serves to wrap together spatial_qc.py and temporal_qc.py with relevant inputs.
"""

from spatial_qc import summary_mask, snr, cnr, fber, efc, artifacts, fwhm, ghost_all
from temporal_qc import mean_dvars_wrapper, mean_fd_wrapper, mean_outlier_timepoints, mean_quality_timepoints


###
# Spatial
###

# anat_type: anat or epi
def run_qc_spatial(anat_path, fg_path, anat_type="anat", 
                   gm_path=None, wm_path=None, csf_path=None, 
                   subject=None, session=None, scan=None):
    # Lead the data
    anat_data       = load_image(anat_path)
    fg_mask         = load_mask(fg_path, anat_path)
    bg_mask         = 1 - fg_mask
    if anat_type != "epi":
        gm_mask         = load_mask(gm_path, anat_path)
        wm_mask         = load_mask(wm_path, anat_path)
        csf_mask        = load_mask(csf_path, anat_path)
    
    # Initialize QC
    qc              = dict()
    if subject:
        qc['subject'] = subject
    if session:
        qc['session'] = session
    if scan:
        qc['scan']    = scan
    
    # Summary Measures
    qc['fg_mean'], qc['fg_std'], qc['fg_size']      = summary_mask(anat_data, fg_mask)
    qc['bg_mean'], qc['bg_std'], qc['bg_size']      = summary_mask(anat_data, bg_mask)
    if anat_type != "epi":
        qc['gm_mean'], qc['gm_std'], qc['gm_size']      = summary_mask(anat_data, gm_mask)
        qc['wm_mean'], qc['wm_std'], qc['wm_size']      = summary_mask(anat_data, wm_mask)
        qc['csf_mean'], qc['csf_std'], qc['csf_size']   = summary_mask(anat_data, csf_mask)
    else:
        qc['gm_mean'], qc['gm_std'], qc['gm_size']      = (None, None, None)
        qc['wm_mean'], qc['wm_std'], qc['wm_size']      = (None, None, None)
        qc['csf_mean'], qc['csf_std'], qc['csf_size']   = (None, None, None)
    
    # SNR
    if anat_type != "epi":
        qc['snr']       = snr(qc['gm_mean'], qc['bg_std'])
    else:
        qc['snr']       = snr(qc['fg_mean'], qc['bg_std'])
    
    # CNR
    if anat_type != "epi":
        qc['cnr']   = cnr(qc['gm_mean'], qc['wm_mean'], qc['bg_std'])
    else:
        qc['cnr']   = None
    
    # FBER
    qc['fber']      = fber(anat_data, fg_mask)
    
    # EFC
    qc['efc']       = efc(anat_data)
    
    # Artifact
    if anat_type != "epi":
        qc['qi1'], _    = artifacts(anat_data, fg_mask, calculate_qi2=False)
    
    # Smoothness in voxels
    tmp             = fwhm(anat_path, fg_path, out_vox=True)
    qc['fwhm_x'], qc['fwhm_y'], qc['fwhm_z'], qc['fwhm'] = tmp
    
    # Ghosting
    if anat_type == "epi":
        tmp         = ghost_all(anat_data, fg_mask)
        qc['ghost_x'], qc['ghost_y'], qc['ghost_z'] = tmp
    else:
        qc['ghost_x'], qc['ghost_y'], qc['ghost_z'] = None, None, None
    
    return qc


###
# Temporal
###

def run_qc_temporal(func_path, mask_path, 
                    motion_matrix_path, 
                    subject=None, session=None, scan=None, 
                    verbose=False):
    # DVARS
    if verbose: print "...dvars"
    mean_dvars  = mean_dvars_wrapper(func_path, mask_path)
    
    # Mean FD (Jenkinson)
    if verbose: print "...mean FD"
    mean_fd     = mean_fd_wrapper(motion_matrix)
    
    # 3dTout
    if verbose: print "...3dTout"
    mean_outlier= mean_outlier_timepoints(func_path, mask_path)
    
    # 3dTqual
    if verbose: print "...3dTqual"
    mean_quality= mean_quality_timepoints(func_path)
    
    # Compile
    qc = {
        "subject":  subject, 
        "session":  session, 
        "scan":     scan, 
        "dvars":    mean_dvars, 
        "fd":       mean_fd, 
        "outlier":  mean_outlier, 
        "quality":  mean_quality
    }
    
    return qc
