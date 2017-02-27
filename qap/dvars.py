import numpy as np


def remove_zero_variance_voxels(func_timeseries, mask):
    """Modify a head mask to exclude timeseries voxels which have zero 
    variance.

    :type func_timeseries: Nibabel data
    :param func_timeseries: The 4D functional timeseries.
    :type mask: Nibabel data
    :param mask: The binary head mask.
    :rtype: Nibabel data
    :return: The binary head mask, but with voxels of zero variance excluded.
    """

    # this needs optimization/refactoring
    for i in range(0, len(func_timeseries)):
        for j in range(0, len(func_timeseries[0])):
            for k in range(0, len(func_timeseries[0][0])):
                var = func_timeseries[i][j][k].var()
                if int(var) == 0:
                    mask[i][j][k] = mask[i][j][k] * 0

    return mask


def load(func_file, mask_file, check4d=True):
    """Load the functional timeseries data from a NIFTI file into Nibabel data
    format, check/validate the data, and remove voxels with zero variance.

    :type func_file: str
    :param func_file: Filepath to the NIFTI file containing the 4D functional
                      timeseries.
    :type mask_file: str
    :param mask_file: Filepath to the NIFTI file containing the binary
                      functional brain mask.
    :type check4d: bool
    :param check4d: (default: True) Check the timeseries data to ensure it is
                    four dimensional.
    :rtype: Nibabel data
    :return: The validated functional timeseries data with voxels of zero
             variance excluded.
    """

    import nibabel as nib
    from qap_utils import raise_smart_exception

    try:
        func_img = nib.load(func_file)
        mask_img = nib.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask = mask_img.get_data()
    func = func_img.get_data().astype(np.float)

    if check4d and len(func.shape) != 4:
        err = "Input functional %s should be 4-dimensional" % func_file
        raise_smart_exception(locals(),err)

    mask_var_filtered = remove_zero_variance_voxels(func, mask)

    func = func[mask_var_filtered.nonzero()].T # will have ntpts x nvoxs
    
    return func


def robust_stdev(func):
    """Compute robust estimation of standard deviation.

    :type func: Nibabel data
    :param func: The functional timeseries data.
    :rtype: float
    :return: The standard deviation value.
    """

    lower_qs = np.percentile(func, 25, axis=0)
    upper_qs = np.percentile(func, 75, axis=0)
    stdev = (upper_qs - lower_qs)/1.349
    return stdev


def ar_nitime(x, order=1, center=False):
    """Derive a model of the noise present in the functional timeseries for 
    the calculation of the standardized DVARS.

    - Borrowed from nipy.algorithms.AR_est_YW. aka "from nitime import
      algorithms as alg".

    :type x: Nibabel data
    :param x: The vector of one voxel's timeseries.
    :type order: int
    :param order: (default: 1) Which lag of the autocorrelation of the
                  timeseries to use in the calculation.
    :type center: bool
    :param center: (default: False) Whether to center the timeseries (to
                   demean it).
    :rtype: float
    :return: The modeled noise value for the current voxel's timeseries.
    """

    from nitime.lazy import scipy_linalg as linalg
    import nitime.utils as utils
    if center:
        x = x.copy()
        x = x - x.mean()
    r_m = utils.autocorr(x)[:order + 1]
    Tm = linalg.toeplitz(r_m[:order])
    y = r_m[1:]
    ak = linalg.solve(Tm, y)
    return ak[0]


def ar1(func, method=ar_nitime):
    """Apply the 'ar_nitime' function across the centered functional 
    timeseries.

    :type func: Nibabel data
    :param func: The functional timeseries data.
    :type method: Python function
    :param method: (default: ar_nitime) The algorithm to use to calculate AR1.
    :rtype: NumPy array
    :return: The vector of AR1 values.
    """
    func_centered = func - func.mean(0)
    ar_vals = np.apply_along_axis(method, 0, func_centered)
    return ar_vals


def calc_dvars(func_file, mask_file, output_all=False):
    """Calculate the standardized DVARS metric.

    :type func_file: str
    :param func_file: The filepath to the NIFTI file containing the functional
                       timeseries.
    :type mask_file: str
    :param mask_file: The filepath to the NIFTI file containing the binary
                      functional brain mask.
    :type output_all: bool
    :param output_all: (default: False) Whether to output all versions of
                       DVARS measure (non-standardized, standardized and
                       voxelwise standardized).
    :rtype: NumPy array
    :return: The output DVARS values vector.
    """

    from qap_utils import raise_smart_exception

    # load data
    func = load(func_file, mask_file)

    # Robust standard deviation
    func_sd = robust_stdev(func)
    
    # AR1
    func_ar1 = ar1(func)

    # Predicted standard deviation of temporal derivative
    func_sd_pd = np.sqrt(2 * (1 - func_ar1)) * func_sd
    diff_sd_mean= func_sd_pd.mean()

    # Compute temporal difference time series 
    func_deriv = np.diff(func, axis=0)

    # DVARS
    # (no standardization)
    dvars_plain = func_deriv.std(1, ddof=1) # TODO: Why are we not ^2 this & getting the sqrt?
    # standardization
    dvars_stdz  = dvars_plain/diff_sd_mean
    # voxelwise standardization
    diff_vx_stdz = func_deriv/func_sd_pd
    dvars_vx_stdz = diff_vx_stdz.std(1, ddof=1)
    
    if output_all:
        try:
            out = np.vstack((dvars_stdz, dvars_plain, dvars_vx_stdz))
        except:
            raise_smart_exception(locals())
    else:
        try:
            out = dvars_stdz.reshape(len(dvars_stdz), 1)
        except:
            raise_smart_exception(locals())
    
    return out
