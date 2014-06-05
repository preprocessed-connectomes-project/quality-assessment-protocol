"""
TODO
"""

import os
import sys
import numpy as np
import nibabel as nb
import pandas as pd
import scipy.ndimage as nd
import scipy.stats as stats
from tempfile import mkdtemp

# DVARS
from dvars import mean_dvars_wrapper

# MeanFD
from CPAC.generate_motion_statistics import calculate_FD_J as fd_jenkinson
def mean_fd_wrapper(in_file, fd_out_file=None):
    # Run in temporary working directory
    tmpdir = mkdtemp()
    os.setwd(tmpdir)
    
    # Compute
    out_file    = fd_jenkinson(in_file)
    fd          = np.loadtxt(out_file)
    mean_fd     = fd.mean()
    
    # Clean-Up
    if fd_out_file:
        os.rename(out_file, fd_out_file)
    else:
        os.remove(out_file)
    os.remove(tmpdir)
    
    return mean_fd

# 3dTout
def outlier_timepoints(func_file, mask_file, out_fraction=True):
    """
    Calculates the number of 'outliers' in a 4D functional dataset, 
    at each time-point.
    
    Will call on AFNI's 3dToutcount.
    
    Parameters
    ----------
    func_file: str
        Path to 4D functional file (could be motion corrected or not??)
    mask_file: str
        Path to functional brain mask
    out_fraction: bool (default: True)
        Whether the output should be a count (False) or fraction (True)
        of the number of masked voxels which are outliers at each time point.
    
    Returns
    -------
    outliers: list
    """
    import commands
    
    opts    = []
    if out_fraction:
        opts.append("-fraction")
    opts.append("-mask %s" % mask_file)
    opts.append(func_file)
    str_opts= " ".join(opts)
    
    # TODO:
    # check if should use -polort 2 (http://www.na-mic.org/Wiki/images/8/86/FBIRNSupplementalMaterial082005.pdf)
    # or -legendre to remove any trend
    cmd     = "3dToutcount %s" % str_opts
    out     = commands.getoutput(cmd)
    
    # Extract time-series in output
    lines   = out.splitlines()
    ## remove general information
    lines   = [ l for l in lines if l[:2] != "++" ]
    ## string => floats
    outliers= [ float(l.strip()) for l in lines ] # note: don't really need strip
    
    return outliers

def mean_outlier_timepoints(*args, **kwrds):
    outliers        = outlier_timepoints(*args, **kwrds)
    mean_outliers   = np.mean(outliers)
    return mean_outliers


# 3dTqual
def quality_timepoints(func_file, automask=True):
    """
    Calculates a 'quality index' for each timepoint in the 4D functional dataset.
    Low values are good and indicate that the timepoint is not very different from the norm.
    """
    import subprocess
    
    opts    = []
    if automask:
        opts.append("-automask")
    opts.append(func_file)
    str_opts= " ".join(opts)
    
    cmd     = "3dTqual %s" % str_opts
    p       = subprocess.Popen(cmd.split(" "), 
                                stdout=subprocess.PIPE, 
                                stderr=subprocess.PIPE)
    out,err = p.communicate()
    
    #import code
    #code.interact(local=locals())
    
    # Extract time-series in output
    lines   = out.splitlines()
    ## remove general information
    lines   = [ l for l in lines if l[:2] != "++" ]
    ## string => floats
    outliers= [ float(l.strip()) for l in lines ] # note: don't really need strip
    
    return outliers

def mean_quality_timepoints(*args, **kwrds):
    qualities       = quality_timepoints(*args, **kwrds)
    mean_qualities  = np.mean(qualities)
    return mean_qualities
