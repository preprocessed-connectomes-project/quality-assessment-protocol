
def raise_smart_exception(local_vars, msg=None):
    """Raise an exception with more information about the traceback, and
    enforce inclusion of the locals().

    :type local_vars: dict
    :param local_vars: Input for locals().
    :type msg: str
    :param msg: (default: None) The custom error message for the exception in
                question.
    """

    import traceback
    e = "\n\nLocal variables:\n%s\n\n%s\n\n" \
        % (str(local_vars), str(traceback.format_exc()))
    if msg:
        e = "%s\n\n%s\n\n" % (e, str(msg))
    raise Exception(e)


def read_nifti_image(nifti_infile):
    """Read a NIFTI file into Nibabel-format image data.

    :type nifti_infile: str
    :param nifti_infile: The filepath of the NIFTI image to read in.
    :rtype: Nibabel image
    :return: Image data in Nibabel format.
    """

    import nibabel as nb

    try:
        nifti_img = nb.load(nifti_infile)
    except:
        err = "\n\n[!] Could not load the NIFTI image using Nibabel:\n" \
              "%s\n\n" % nifti_infile
        raise_smart_exception(locals(), err)

    return nifti_img


def get_masked_data(data, mask):
    """Extract data from a NIFTI file and apply a mask to it.

    :type data: NumPy array or str
    :param data: The data (either in a NumPy array, or the filepath to a NIFTI
                 file).
    :type mask: NumPy array or str
    :param mask: The mask data (either in a NumPy array, or the filepath to a
                 NIFTI file).
    :type files: bool
    :param files: Whether or not the input parameters are filepath strings or
                  NumPy array objects.
    """

    import numpy as np
    import nibabel as nb

    if isinstance(data, basestring):
        img = nb.load(data)
        data = img.get_data()
    if isinstance(mask, basestring):
        mask_img = nb.load(mask)
        mask = mask_img.get_data()

    data = np.asarray(data, dtype=float)
    mask = np.asarray(mask, dtype=float)

    if len(data.shape) > 3:
        # if timeseries
        masked_data = np.asarray([volume * mask.T for volume in data.T]).T
    else:
        masked_data = np.asarray(data * mask)

    return masked_data


def create_zero_variance_mask(func_timeseries, mask):
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

    import numpy as np
    import nibabel as nib

    try:
        func_img = nib.load(func_file)
        mask_img = nib.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask = mask_img.get_data()
    func = func_img.get_data().astype(np.float)

    if check4d and len(func.shape) != 4:
        err = "Input functional %s should be 4-dimensional" % func_file
        raise_smart_exception(locals(), err)

    mask_var_filtered = create_zero_variance_mask(func, mask)

    func = func[mask_var_filtered.nonzero()].T  # will have ntpts x nvoxs

    return func


def load_timeseries(func_file, mask_file, check4d=True):
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

    import numpy as np
    import nibabel as nb

    try:
        func_img = nb.load(func_file)
        mask_img = nb.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask = mask_img.get_data()
    func = func_img.get_data().astype(np.float)

    if check4d and len(func.shape) != 4:
        err = "Input functional %s should be 4-dimensional" % func_file
        raise_smart_exception(locals(), err)

    mask_var_filtered = create_zero_variance_mask(func, mask)

    func_masked = get_masked_data(func, mask_var_filtered)

    return func, func_masked


def robust_stdev(func):
    """Compute robust estimation of standard deviation.

    :type func: Nibabel data
    :param func: The functional timeseries data.
    :rtype: float
    :return: The standard deviation value.
    """

    import numpy as np

    lower_qs = np.percentile(func, 25, axis=3)
    upper_qs = np.percentile(func, 75, axis=3)
    stdev = (upper_qs - lower_qs) / 1.349
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


def calc_ar1(func, method=ar_nitime):
    """Apply the 'ar_nitime' function across the centered functional
    timeseries.

    :type func: Nibabel data
    :param func: The functional timeseries data.
    :type method: Python function
    :param method: (default: ar_nitime) The algorithm to use to calculate AR1.
    :rtype: NumPy array
    :return: The vector of AR1 values.
    """

    import numpy as np

    ar1 = np.zeros(func.shape[0:3])

    for x in range(0, len(func)):
        for y in range(0, len(func[x])):
            for z in range(0, len(func[x][y])):
                try:
                    ar1[x][y][z] = method(func[x][y][z], center=True)
                except:
                    continue

    return ar1


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

    import numpy as np

    # load data
    # func is shape (x, y, z, t)
    # func_masked is shape (t, all voxels flattened)
    func, func_masked = load_timeseries(func_file, mask_file)

    # get binary mask of functional brain
    binary_mask_data = read_nifti_image(mask_file).get_data()

    # Robust standard deviation
    func_sd = robust_stdev(func)

    # AR1
    func_ar1 = calc_ar1(func)

    # Predicted standard deviation of temporal derivative
    func_sd_pd = np.sqrt(2 * (1 - func_ar1)) * func_sd

    # Compute temporal difference squared time series
    func_deriv = np.diff(func, axis=3)
    func_deriv_sq = func_deriv ** 2

    # Flip func_deriv_sq from (x, y, z, t) to (t, x, y, z)
    flipped_deriv = []
    for x in range(0, len(func_deriv_sq.T)):
        flipped_deriv.append(func_deriv_sq.T[x].T)
    flipped_deriv = np.asarray(flipped_deriv)

    # in (t, x, y, z), find the mean of nonzero voxels for all (x, y, z)
    # voxels within each "t" (each volume)
    diff_var = []
    for ts in flipped_deriv:
        diff_var.append(ts[ts.nonzero()].mean())
    diff_var = np.asarray(diff_var)

    diff_sd = get_masked_data(func_sd_pd, binary_mask_data)
    diff_sd_mean = diff_sd[diff_sd.nonzero()].mean()

    # calculate standardized DVARS
    dvars_stdz = np.sqrt(diff_var) / diff_sd_mean

    if output_all:
        # voxelwise standardization
        diff_vx_stdz = func_deriv / func_sd_pd
        dvars_vx_stdz = diff_vx_stdz.std(1, ddof=1)
        try:
            out = np.vstack((dvars_stdz, diff_var, dvars_vx_stdz))
        except Exception as e:
            raise_smart_exception(locals(), e)
    else:
        try:
            out = dvars_stdz.reshape(len(dvars_stdz), 1)
        except Exception as e:
            raise_smart_exception(locals(), e)

    return out


def main():

    import os
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("func_file", type=str,
                        help="filepath to the functional timeseries file")
    parser.add_argument("mask_file", type=str,
                        help="filepath to the functional mask file")

    args = parser.parse_args()

    dvars_vector = calc_dvars(args.func_file, args.mask_file)

    filepath = os.path.abspath(args.func_file)
    #data_id = "_".join(filepath.split("/")[-1].split("_")[0:2])

    data_id = "_".join(filepath.split("/")[-2:]).rstrip(".nii.gz")

    out_string = data_id
    for num in dvars_vector:
        out_string = out_string + "," + str(num[0])

    with open(os.path.join(os.getcwd(), "dvars.csv"), "a") as f:
        print >>f, out_string


if __name__ == "__main__":
    main()