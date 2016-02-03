
import os
import sys
import numpy as np
import nibabel as nb
import pandas as pd
import scipy.ndimage as nd
import scipy.stats as stats
from tempfile import mkdtemp
import shutil

# DVARS
from dvars import mean_dvars_wrapper


def pass_floats(output_string):

    # for parsing AFNI output strings

    lines = output_string.splitlines()

    values_list = []

    for l in lines:

        try:
            val = float(l)
            values_list.append(val)
        except:
            pass


    return values_list



def calculate_percent_outliers(values_list):

    import numpy as np
    from workflow_utils import raise_smart_exception

    try:

        # calculate the IQR
        sorted_values = sorted(values_list)

        third_qr, first_qr = np.percentile(sorted_values, [75, 25])
        IQR = third_qr - first_qr

        # calculate percent outliers
        third_qr_threshold = third_qr + (1.5 * IQR)
        first_qr_threshold = first_qr - (1.5 * IQR)

        high_outliers = \
            [val for val in sorted_values if val > third_qr_threshold]
        low_outliers = \
            [val for val in sorted_values if val < first_qr_threshold]

        total_outliers = high_outliers + low_outliers

        percent_outliers = \
            float(len(total_outliers)) / float(len(sorted_values))

    except:
        
        raise_smart_exception(locals())


    return percent_outliers, IQR



def fd_jenkinson(in_file, rmax=80., out_file=None):
    '''
    @ Krsna
    May 2013
    compute
    1) Jenkinson FD from 3dvolreg's *.affmat12.1D file from -1Dmatrix_save
    option input: subject ID, rest_number, name of 6 parameter motion
    correction file (an output of 3dvolreg) output: FD_J.1D file
    Assumptions: 1) subject is available in BASE_DIR
    2) 3dvolreg is already performed and the 1D motion parameter and 1D_matrix
    file file is present in sub?/rest_? called as --->'lfo_mc_affmat.1D'

    Method to calculate Framewise Displacement (FD) calculations
    (Jenkinson et al., 2002)
    Parameters; in_file : string
                rmax : float
                The default radius (as in FSL) of a sphere represents the brain
    Returns; out_file : string
    NOTE: infile should have one 3dvolreg affine matrix in one row -
    NOT the motion parameters
    '''

    import numpy as np
    import os
    import os.path as op
    from shutil import copyfile
    import sys
    import math
    from qap.workflow_utils import raise_smart_exception

    if out_file is None:
        fname, ext = op.splitext(op.basename(in_file))
        out_file = op.abspath('%s_fdfile%s' % (fname, ext))

    # if in_file (coordinate_transformation) is actually the rel_mean output
    # of the MCFLIRT command, forward that file
    if 'rel.rms' in in_file:
        copyfile(in_file, out_file)
        return out_file

    try:
        pm_ = np.genfromtxt(in_file)
    except:
        raise_smart_exception(locals())

    original_shape = pm_.shape
    pm = np.zeros((pm_.shape[0], pm_.shape[1] + 4))
    pm[:, :original_shape[1]] = pm_
    pm[:, original_shape[1]:] = [0.0, 0.0, 0.0, 1.0]

    # rigid body transformation matrix
    T_rb_prev = np.matrix(np.eye(4))

    flag = 0
    X = [0]  # First timepoint
    for i in range(0, pm.shape[0]):
        # making use of the fact that the order of aff12 matrix is "row-by-row"
        T_rb = np.matrix(pm[i].reshape(4, 4))

        if flag == 0:
            flag = 1
        else:
            M = np.dot(T_rb, T_rb_prev.I) - np.eye(4)
            A = M[0:3, 0:3]
            b = M[0:3, 3]

            FD_J = math.sqrt(
                (rmax * rmax / 5) * np.trace(np.dot(A.T, A)) + np.dot(b.T, b))
            X.append(FD_J)

        T_rb_prev = T_rb

    try:
        np.savetxt(out_file, np.array(X))
    except:
        raise_smart_exception(locals())


    return out_file



def outlier_timepoints(func_file, out_fraction=True):

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
    import re
    import numpy as np
    from workflow_utils import raise_smart_exception

    opts = []
    if out_fraction:
        opts.append("-fraction")
    #opts.append("-mask %s" % mask_file)
    opts.append(func_file)
    str_opts = " ".join(opts)

    # TODO:
    # check if should use -polort 2 (http://www.na-mic.org/Wiki/images/8/86/FBIRNSupplementalMaterial082005.pdf)
    # or -legendre to remove any trend
    cmd = "3dToutcount %s" % str_opts

    try:
        out = commands.getoutput(cmd)
    except:
        err = "[!] QAP says: Something went wrong with running AFNI's " \
              "3dToutcount."
        raise_smart_exception(locals(),err)

    # Extract time-series in output
    #lines = out.splitlines()

    # remove general information and warnings
    #outliers = [float(l) for l in lines if re.match("[0-9]+$", l.strip())]
    outliers = pass_floats(out)


    return outliers



def mean_outlier_timepoints(*args, **kwrds):

    outliers = outlier_timepoints(*args, **kwrds)

    # calculate the outliers of the outliers! AAHH!
    percent_outliers, IQR = calculate_percent_outliers(outliers)

    mean_outliers = np.mean(outliers)
    
    return mean_outliers, percent_outliers, IQR



def quality_timepoints(func_file):

    """
    Calculates a 'quality index' for each timepoint in the 4D functional
    dataset. Low values are good and indicate that the timepoint is not very
    different from the norm.
    """

    import subprocess
    from workflow_utils import raise_smart_exception

    opts = []
    opts.append(func_file)
    str_opts = " ".join(opts)

    cmd = "3dTqual %s" % str_opts

    try:
        p = subprocess.Popen(cmd.split(" "),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()
    except:
        err = "[!] QAP says: Something went wrong with running AFNI's " \
              "3dTqual."
        raise_smart_exception(locals(),err)

    # Extract time-series in output
    #lines = out.splitlines()
    # remove general information
    #lines = [l for l in lines if l[:2] != "++"]
    # string => floats
    #quality = [float(l.strip())
    #            for l in lines]  # note: don't really need strip

    quality = pass_floats(out)

    # get percent outliers and IQR
    percent_outliers, IQR = calculate_percent_outliers(quality)


    return quality, percent_outliers, IQR



def mean_quality_timepoints(*args, **kwrds):
    qualities, percent_outliers, IQR = quality_timepoints(*args, **kwrds)
    mean_qualities = np.mean(qualities)
    return mean_qualities, percent_outliers, IQR



def global_correlation(func_motion, func_mask):

    import scipy
    import numpy as np
    from dvars import load

    zero_variance_func = load(func_motion, func_mask)

    list_of_ts = zero_variance_func.transpose()

    # get array of z-scored values of each voxel in each volume of the
    # timeseries
    demeaned_normed = []

    for ts in list_of_ts:

        demeaned_normed.append(scipy.stats.mstats.zscore(ts))

    demeaned_normed = np.asarray(demeaned_normed)

    # make an average of the normalized timeseries, into one averaged
    # timeseries, a vector of N volumes
    volume_list = demeaned_normed.transpose()

    avg_ts = []

    for voxel in volume_list:

        avg_ts.append(voxel.mean())

    avg_ts = np.asarray(avg_ts)

    # calculate the global correlation
    gcor = (avg_ts.transpose().dot(avg_ts)) / len(avg_ts)

    return gcor
