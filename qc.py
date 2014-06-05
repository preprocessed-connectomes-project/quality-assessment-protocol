"""
Serves to wrap together spatial_qc.py and temporal_qc.py with relevant inputs.

Two layers of wrappers exist:
1. Takes in anatomical/functional inputs and runs QC
2. Takes in CPAC inputs and calls #1.
"""

import os, sys, re
import numpy as np
import nibabel as nib
from glob import glob
from tempfile import mkstemp

from utils import gen_file_map
from spatial_qc import summary_mask, snr, cnr, fber, efc, artifacts, fwhm, ghost_all
from temporal_qc import mean_dvars_wrapper, mean_fd_wrapper, mean_outlier_timepoints, mean_quality_timepoints


###
# Loading Data Functions
###

def load_image(image_file):
    img         = nib.load(image_file)
    dat         = img.get_data()
    return(dat)

def load_mask(mask_file, ref_file):
    mask_img    = nib.load(mask_file)
    mask_dat    = mask_img.get_data()
    ref_img     = nib.load(ref_file)
    
    # Check that the specified mask is binary.
    mask_vals   = np.unique(mask_dat)
    if (mask_vals.size != 2) or not (mask_vals == [0, 1]).all():
        print("Error: Mask is not binary")
        raise Exception("")
    
    # Verify that the mask and anatomical images have the same dimensions.
    if ref_img.shape != mask_img.shape:
        print("Error: Mask and anatomical image are different dimensions")
        raise Exception("")

    # Verify that the mask and anatomical images are in the same space (have the samme affine matrix)
    if (mask_img.get_affine() == ref_img.get_affine()).all == False:
        print("Error: Mask and anatomical image are not in the same space")
        raise Exception("")
    
    return mask_dat

def load_func(func_file, mask_file):
    func_dat    = load_image(func_file)
    mask        = load_mask(mask_file)
    func        = func_dat.astype(np.float)
    if len(func.shape) != 4:
        raise Exception("Input functional %s should be 4-dimensional" % func_file)
    func        = func[mask.nonzero()].T # will have ntpts x nvoxs
    return(func)


###
# Spatial
###

def run_qc_spatial(anat_type, anat_path, 
                   fg_path, gm_path=None, wm_path=None, csf_path=None, 
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

class CpacPaths(object):
    def __init__(self, pipeline_dir, subject_id, session_id, sink_dir=None):
        # Set paths for CPAC
        if sink_dir is None:
            sink_dir    = os.path.dirname(pipeline_dir)
        
        sub_dir         = os.path.join(pipeline_dir, "%s_session_%i" % (subject_id, session_id))
        if not os.path.exists(sub_dir):
            raise Exception("Error: subject '%s' directory '%s' not found" % (subject_id, sub_dir))
        
        path_files      = os.path.join(sub_dir, "path_files_here", "*.txt")
        file_map        = gen_file_map(sink_dir, path_files)
        file_map        = [ v[0] for k,v in file_map.iteritems() if len(v) > 0 ]
        self.file_map   = file_map
        
        return
    
    def extract_all(self, resource_id):
        paths = [ v[3] for v in self.file_map if v[2] == resource_id ]
        return paths
    
    def extract(self, resource_id):
        paths = self.extract_all(resource_id)
        if len(paths) > 1:
            raise Exception("Error: there is more than one path for %s" % resource_id)
        elif len(paths) == 0:
            raise Exception("Error: no path found for %s" % resource_id)        
        return paths[0]

class CpacPaths2(object):
    def __init__(self, pipeline_dir, subject_id, session_id, sink_dir=None):
        # Set paths for CPAC
        if sink_dir is None:
            sink_dir    = os.path.dirname(pipeline_dir)
        
        sub_dir         = os.path.join(pipeline_dir, "%s_session_%i" % (subject_id, session_id))
        if not os.path.exists(sub_dir):
            raise Exception("Error: subject '%s' directory '%s' not found" % (subject_id, sub_dir))
        
        self.sub_dir    = sub_dir
        
        return
    
    def extract(self, resource_id, out_file="*"):
        raw_path = os.path.join(self.sub_dir, resource_id, out_file)
        paths = glob(raw_path)
        if len(paths) > 1:
            raise Exception("Error: there is more than one path for %s" % resource_id)
        elif len(paths) == 0:
            raise Exception("Error: no path found for %s" % resource_id)        
        return paths[0]

def cpac_qc_spatial(subject_id, pipeline_dir, session_id=1, sink_dir=None):
    """
    For now assume T1
    """
    print "Running spatial QC for %s - %i" % (subject_id, session_id)
    
    # Generate head mask (for now assume it exists and use that)
    def find_or_generate_fg_mask():
        pipeline    = os.path.basename(pipeline_dir)
        base_dir    = os.path.dirname(pipeline_dir)
        fg_mask_path= os.path.join(base_dir, "sym_links", pipeline, 
            "_compcor_ncomponents_5_linear1.global1.motion1.quadratic1.compcor1.CSF_0.96_GM_0.7_WM_0.96", 
            "%s_session_%i" % (subject_id, session_id), "scan", "anat", "qc_head_mask.nii.gz")
        
        if not os.path.exists(fg_mask_path):
            print "Error: head mask '%s' doesn't exist" % fg_mask_path
            raise Exception("")
        
        return fg_mask_path
        
    try:
        # Get all the paths
        #cpac_paths  = CpacPaths(pipeline_dir, subject_id, session_id, sink_dir)
        cpac_paths  = CpacPaths2(pipeline_dir, subject_id, session_id, sink_dir)
        anat_path   = cpac_paths.extract("anatomical_reorient")
        fg_path     = cpac_paths.extract("anatomical_head", "30_head_clean_mask.nii.gz")
        #fg_path     = find_or_generate_fg_mask()
        gm_path     = cpac_paths.extract("anatomical_gm_mask", "*/*.nii.gz")
        wm_path     = cpac_paths.extract("anatomical_wm_mask", "*/*.nii.gz")
        csf_path    = cpac_paths.extract("anatomical_csf_mask", "*/*.nii.gz")
        
        qc = run_qc_spatial("t1", anat_path, fg_path, 
                            gm_path, wm_path, csf_path, 
                            subject=subject_id, session=session_id)
    except:
        e = sys.exc_info()[1]
        print "Error: subject %s with session %i failed" % (subject_id, session_id)
        print "Details: %s" % e
        qc = None
        #raise
    
    return qc

def cpac_qc_spatial_epi(subject_id, pipeline_dir, session_id=1, scan_id=1, sink_dir=None):
    """
    For now assume T1
    """
    print "Running spatial QC for %s - %i" % (subject_id, session_id)
    
    anat_path = None
    
    try:
        # Get all the paths
        #cpac_paths  = CpacPaths(pipeline_dir, subject_id, session_id, sink_dir)
        #anat_path   = cpac_paths.extract("mean_functional", "_scan_rest_1_rest/*.nii.gz")
        #fg_path     = find_or_generate_fg_mask()
        cpac_paths  = CpacPaths2(pipeline_dir, subject_id, session_id, sink_dir)
        fg_path     = cpac_paths.extract("functional_brain_mask", "_scan_rest_%i_rest/*.nii.gz" % scan_id)
        
        # Create mean EPI with background noise
        func_path   = cpac_paths.extract("motion_correct", "_scan_rest_%i_rest/*.nii.gz" % scan_id)
        _,anat_path = mkstemp(suffix=".nii.gz", prefix="tmp_mean_epi_")
        cmd         = "fslmaths %s -Tmean %s" % (func_path, anat_path)
        os.system(cmd)
        
        #import code
        #code.interact(local=locals())
        
        qc = run_qc_spatial("epi", anat_path, fg_path, 
                            subject=subject_id, session=session_id, scan=scan_id)
    except:
        e = sys.exc_info()[1]
        print "Error: subject %s with session %i failed" % (subject_id, session_id)
        print "Details: %s" % e
        qc = None
        #raise
    finally:
        if anat_path:
            os.remove(anat_path)
    
    return qc


def test_cpac_qc_spatial():
    from qc import cpac_qc_spatial
    # Try out with ABIDE
    subject_id      = "0050002"
    session_id      = 1
    pipeline_dir    = '/data/Projects/ABIDE_Initiative/CPAC/Output_2013-11-22/pipeline_MerrittIsland'
    qc              = cpac_qc_spatial(subject_id, pipeline_dir, session_id)
    
    #cols            = ["subject", "session", "snr", "cnr", ]
    
    
    import pandas as pd
    pd.DataFrame(qc)


def test_run_many_subs():
    import numpy as np
    from qc import cpac_many_subs
    
    session_id      = 1
    pipeline_dir    = '/data/Projects/ABIDE_Initiative/CPAC/Output_2013-11-22/pipeline_MerrittIsland'
    
    sublist_file = "/data/Projects/ABIDE_Initiative/scripts/for_grant/10_raw_anat/z_sublist.txt"
    sublist = np.loadtxt(sublist_file)
    sublist = [ "%07i" % s for s in sublist ]
    
    df1 = cpac_many_subs(sublist[:3], pipeline_dir, session_id, ncores=3)
    df2 = cpac_many_subs(sublist[:3], pipeline_dir, session_id, ncores=1)
   
     
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
    #if verbose: print "...mean FD"
    #mean_fd     = mean_fd_wrapper(motion_matrix)
    mean_fd     = None
    
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

class CpacFunctionalPaths(object):
    def __init__(self, pipeline_dir, subject_id, session_id, scan_id, sink_dir=None):
        # Set paths for CPAC
        if sink_dir is None:
            sink_dir    = os.path.dirname(pipeline_dir)
        
        sub_dir         = os.path.join(pipeline_dir, "%s_session_%i" % (subject_id, session_id))
        if not os.path.exists(sub_dir):
            raise Exception("Error: subject '%s' directory '%s' not found" % (subject_id, sub_dir))
        
        path_files      = os.path.join(sub_dir, "path_files_here", "*.txt")
        file_map        = gen_file_map(sink_dir, path_files)
        file_map        = [ v[0] for k,v in file_map.iteritems() if len(v) > 0 ]
        self.file_map   = file_map
        
        self.scan_id    = str(scan_id)
        
        return
    
    def extract_all(self, resource_id, scan_id=None):
        if scan_id is None:
            scan_id = self.scan_id
        else:
            scan_id = str(scan_id)
        paths = [ v[3] for v in self.file_map if v[2] == resource_id ]
        paths = [ p for p in paths if re.search("/_scan.*%s.*/" % scan_id, p) ]
        return paths
    
    def extract(self, resource_id):
        paths = self.extract_all(resource_id)
        if len(paths) > 1:
            raise Exception("Error: there is more than one path for %s" % resource_id)
        elif len(paths) == 0:
            raise Exception("Error: no path found for %s" % resource_id)
        return paths[0]

class CpacFunctionalPaths2(object):
    def __init__(self, pipeline_dir, subject_id, session_id, scan_id, sink_dir=None):
        # Set paths for CPAC
        if sink_dir is None:
            sink_dir    = os.path.dirname(pipeline_dir)
        
        sub_dir         = os.path.join(pipeline_dir, "%s_session_%i" % (subject_id, session_id))
        if not os.path.exists(sub_dir):
            raise Exception("Error: subject '%s' directory '%s' not found" % (subject_id, sub_dir))
        self.sub_dir    = sub_dir
        
        self.scan_id    = str(scan_id)
        
        return
    
    def extract(self, resource_id):
        raw_path = os.path.join(self.sub_dir, resource_id, "_scan*_%s_*" % self.scan_id, "*")
        paths = glob(raw_path)
        if len(paths) > 1:
            raise Exception("Error: there is more than one path for %s" % resource_id)
        elif len(paths) == 0:
            raise Exception("Error: no path found for %s" % resource_id)
        return paths[0]

def cpac_qc_temporal(subject_id, pipeline_dir, session_id=1, scan_id=1, sink_dir=None):
    """
    """
    print "Running temporal QC for %s - %i" % (subject_id, session_id)
        
    try:
        # Get all the paths
        #cpac_paths  = CpacFunctionalPaths(pipeline_dir, subject_id, session_id, scan_id, sink_dir)
        cpac_paths  = CpacFunctionalPaths2(pipeline_dir, subject_id, session_id, scan_id, sink_dir)
        func_path   = cpac_paths.extract("motion_correct")
        mask_path   = cpac_paths.extract("functional_brain_mask")
        #matrix_path = extract_path("?")
        matrix_path = ""
        
        qc = run_qc_temporal(func_path, mask_path, 
                            matrix_path, 
                            subject=subject_id, session=session_id, scan=scan_id, verbose=False)
    except:
        e = sys.exc_info()[1]
        print "Error: subject %s with session %i failed" % (subject_id, session_id)
        print "Details: %s" % e
        qc = None
        #raise
    
    return qc


###
# Wrapper (Combines)
###

def cpac_many_subs(qc_type, subjects, pipeline_dir, ncores=1, **kwrds):
    """
    cpac_many_subs("spatial", subjects, pipeline_dir)
    """
    import pandas as pd
    from multiprocessing import Pool
    from functools import partial
    
    if qc_type == "spatial":
        apply_qc = cpac_qc_spatial  
    elif qc_type == "spatial_epi":
        apply_qc = cpac_qc_spatial_epi
    elif qc_type == "temporal":
        apply_qc = cpac_qc_temporal
    else:
        print "Error: qc_type '%s' not recognized" % qc_type
        raise Exception("")
    
    if ncores == 1:
        qcs     = [ apply_qc(subj, pipeline_dir, **kwrds) for subj in subjects ]
    else:
        partial_apply_qc = partial(apply_qc, 
                                    pipeline_dir=pipeline_dir, 
                                    **kwrds)
        pool    = Pool(processes=ncores)
        qcs     = pool.map(partial_apply_qc, subjects)
        pool.close()
        pool.join()
    
    bad_subs    = [ subjects[i] for i in xrange(len(qcs)) if qcs[i] is None ]
    qcs         = [ qc for qc in qcs if qc is not None ]
    
    df = pd.DataFrame(qcs)
    
    return (bad_subs, df)
