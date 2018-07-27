
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
    t_m = linalg.toeplitz(r_m[:order])
    y = r_m[1:]
    ak = linalg.solve(t_m, y)
    return ak[0]


def scale_then_center(x, scale=1):
    """
    Multiply vector x by scalar scale and then subtract the mean.

    :param x: vector to be operated on.
    :param scale: scalar to multiply against vector
    :return: scaled then centered vector
    """
    scaled_x = scale * x
    mu = scaled_x.mean()
    return scaled_x - mu


def quantile(x, q=0.5):
    """
    calculate percentile q of x, this is designed to be more consistent with the Matlab function quantile

    :param x: vector of length > 1
    :param q: value between 0 and 1 that indicates the desired percentile
    :return: the percentile, if lies between consecutive numbers of x, then interpolation is used
       (pct_ndx - min_ndx)*x[min_ndx] + (max_ndx - pct_ndx)*x[max_ndx], if q=0 the min is returned,
       if q = 1 max is returned, and if q=0.5 the median is returned
    """
    import numpy as np

    if q < 0 or q > 1:
        raise ValueError('q {0} is out of range, should be a fraction from 0 to 1'.format(q))

    if len(x) < 2:
        raise ValueError('quantile expects a vector, not a scalar')

    if q == 1:
        out = np.max(x)
    elif q == 0.5:
        out = np.median(x)
    elif q == 0:
        out = np.min(x)

    else:

        t_x = np.copy(x)
        t_x.sort()

        q_ndx = q * (len(t_x))
        min_ndx = np.floor(q_ndx)
        max_ndx = np.ceil(q_ndx)

        if max_ndx == min_ndx:
            out = t_x[q_ndx]
        else:
            out = (q_ndx - min_ndx)*t_x[int(min_ndx)] + (max_ndx-q_ndx)*t_x[int(max_ndx)]

    return out


def calc_dse(functional_data_filename, mask_filename, debug=False):
    """ Calculate the standardized DVARS metric and derivatives from:

    Afyouni, S. , Nichols, T. (2018). Insight and inference for DVARS, NeuroImage, 172, 291-312.
        https://doi.org/10.1016/j.neuroimage.2017.12.098.

    Nichols, T. (2012, Oct 28). Standardizing DVARS. Retrieved from
        http://blogs.warwick.ac.uk/nichols/entry/standardizing_dvars.

    which extend the definition of DVARS popularized in:

    Power, J. D., Barnes, K. A, Snyder, A. Z., Schlaggar, B. L., & Petersen, S. E. (2012). Spurious but systematic
        correlations in functional connectivity MRI networks arise from subject motion. NeuroImage, 59(3), 2142-54.
        https://doi.org/10.1016/j.neuroimage.2011.10.018


    DVARS is a measure of the amount of frame-to-frame intensity variation in an fMRI dataset. This information can
    be used to evaluate the overall quality of a dataset and to identify frames to be removed from the data.

    The calculation of DVARS, Avar, Delta Pct DVAR and statistical tests are consistent with the DVARS_calc.m from
        https://github.com/asoroosh/DVARS.

    Relative DVARS is consistent with the original FSL based script DVARS.sh from
        http://go.warwick.ac.uk/tenichols/scripts/fsl/#DVARS

    Relative DVARS values vary between the Matlab and FSL implementations, but the results are perfectly correlated.

    :type functional_data_filename: str
    :param functional_data_filename: The filepath to the NIFTI file containing the functional
                       timeseries.
    :type mask_filename: str
    :param mask_filename: The filepath to the NIFTI file containing the binary
                      functional brain mask.
    :type debug: bool
    :param debug: a flag indicating whether extra values should be added to the output dictionary
                      to assist with debugging
    :rtype: dictionary
    :return: dictionary containing values calculated from time course data.

        'DSE_Total_Variance': The total variance in the data after scaling and demeaning (Avar)
        'DSE_Autocorrelation': A representation of the amount of autocorrelation in the data (1 - abs(Svar-Dvar)/Avar).
            If there is a lot of autocorrelation in the data, Svar >> Dvar and this value will be close to 0. With less
            autocorrelation (or if the data has been cleaned), Svar and Dvar will become closer, and this value will be
            close to 1.
        'DSE_Relative_DVARS': The original DVARS time series scaled to be comparable across datasets and scanners.
        'DSE_Delta_Percent_DVAR': An estimate of the difference of a frame's DVAR and the expected DVAR for the dataset.
            This is useful for identifying frames that are outliers.
        'DSE_Delta_Percent_DVAR_90_percentile': The 90th percentile of Delta_Percent_DVAR, gives an estimate of how much
            DVAR would remain if the worse 10% of frames were removed.
        'DSE_P_Values': The likelihood that a DVARS value as extreme as the one occurred could happen by chance.
            Statistical test to determine whether a frame is an outlier.
        'DSE_Outliers': The frames that are considered outliers using the criteria p value < 0.05 and the
            practical threshold Delta_Percent_DVAR > 5%.
        'DSE_Outlier_Count': The number of outliers determined by the above criteria.

        If debug is True the following parameters will also be included in the dictionary. These variables are named
        to make them consistent with the variables output by the DSE Matlab toolbox
        'DSE_DBG_dim': The dimensions of the functional imaging data after it has been conformed (# voxels x # frames)
        'DSE_DBG_GrandMean0': The grand mean of the functional imaging data prior to scaling
        'DSE_DBG_Scale2Norm': The scalar value that was used in scaling the functional imaging data
        'DSE_DBG_GrandMean': The grand mean of the functional imaging data after scaling
        'DSE_DBG_DVARS2': DVARS squared timeseries, used in calculating the statistical tests.
        'DSE_DBG_E': Robust mean of DVARS2, used in calculating chi^2 statistic.
        'DSE_DBG_PowerTrans': Order of the power transform used for the robust calculation of the DVARS2 variance.
        'DSE_DBG_S': Robust standard deviation of DVARS2, used in calculating chi^2 statistic.
        'DSE_DBG_nu': Spatial degrees of freedom for the chi^2 test.
        'DSE_DBG_c': Chi^2 statistic used to test if a frame is an outlier.
        'DSE_DBG_Mu0': Robust estimate of the expected value of DVARS2.

    """

    import numpy as np
    import nibabel as nb
    from scipy import stats

    dvars_out_dictionary = {}

    # load and check functional data and conform it to our needs (float, 2D vx x time array)
    functional_data_image = nb.load(functional_data_filename)

    if len(functional_data_image.shape) != 4:
        raise ValueError("DVARS expects the input functional data to be 4D.")

    functional_data = functional_data_image.get_data().astype(np.float).reshape(
        np.prod(functional_data_image.shape[0:3]), functional_data_image.shape[3])

    # load, check, and conform the mask
    mask_image = nb.load(mask_filename)

    if len(mask_image.shape) != 3:
        raise ValueError("DVARS expects the input mask to be 3D.")

    if mask_image.shape != functional_data_image.shape[:3]:
        raise ValueError("Mask image shape ({0}) does not match "
                         "functional data shape ({1}).".format(mask_image.shape,
                                                               functional_data_image.shape))

    mask_data = mask_image.get_data().reshape(np.prod(mask_image.shape))

    # binarize the mask, allows us to work with non-binary input images
    # assume that nonzero is in -brain and zero is out of brain
    mask_data[np.isclose(mask_data, 0)] = 0
    mask_data[mask_data != 0] = 1

    # update mask to exclude voxels that have zero variance in the functional data
    mask_data[np.isclose(functional_data.var(1), 0)] = 0

    # mask functional data
    functional_data = functional_data[mask_data == 1, :]
    if debug is True:
        dvars_out_dictionary["DSE_DBG_dim"] = functional_data.shape

    # scale values so that grand mean == 100 (keeping consistent with the DSE Matlab toolbox)
    if debug is True:
        dvars_out_dictionary["DSE_DBG_GrandMean0"] = functional_data.mean()
        dvars_out_dictionary["DSE_DBG_Scale2Norm"] = 100 / functional_data.mean()
    functional_data = (100.0 / functional_data.mean()) * functional_data
    if debug is True:
        dvars_out_dictionary["DSE_DBG_GrandMean"] = functional_data.mean()

    # center the data by subtracting the mean
    functional_data = np.apply_along_axis(lambda x: scale_then_center(x), 1, functional_data)

    # == calculate DVARS, avar, dvar, and svar based measures

    # Compute temporal difference squared time series
    functional_avar_timeseries = np.square(functional_data).mean(0)
    functional_avar = functional_avar_timeseries.mean()

    # store variable for output
    dvars_out_dictionary['DSE_Total_Variance'] = functional_avar

    functional_dvar_timeseries = .25 * np.square(functional_data[:, 1:] - functional_data[:, 0:-1]).mean(0)
    functional_svar_timeseries = .25 * np.square(functional_data[:, 1:] + functional_data[:, 0:-1]).mean(0)

    # 1 - abs(%dvar - %svar)
    dvars_out_dictionary['DSE_Autocorrelation'] = 1.0 - (
                abs(functional_dvar_timeseries.mean() - functional_svar_timeseries.mean()) / functional_avar)

    # time series for identifying outliers
    functional_delta_pct_dvar = 100 * (
                functional_dvar_timeseries - np.median(functional_dvar_timeseries)) / functional_avar
    dvars_out_dictionary['DSE_Delta_Percent_DVAR'] = functional_delta_pct_dvar
    dvars_out_dictionary['DSE_Delta_Percent_DVAR_90_percentile'] = np.percentile(functional_delta_pct_dvar, 90)

    # p-values for identifying outliers
    functional_dvars_square_timeseries = np.square(functional_data[:, 0:-1] - functional_data[:, 1:]).mean(0)
    if debug is True:
        dvars_out_dictionary['DSE_DBG_DVARS2'] = functional_dvars_square_timeseries

    dvars_robust_mean = np.median(functional_dvars_square_timeseries)
    if debug is True:
        dvars_out_dictionary['DSE_DBG_E'] = dvars_robust_mean

    transform_power = 1.0 / 3.0
    if debug is True:
        dvars_out_dictionary["DSE_DBG_PowerTrans"] = transform_power

    power_transformed_dvars2 = np.power(functional_dvars_square_timeseries, transform_power)
    power_transformed_dvars2_median = np.median(power_transformed_dvars2)
    half_iqr_sd_dvars2 = 2.0 * (
                np.median(power_transformed_dvars2) - quantile(power_transformed_dvars2, .25)) / 1.349

    # (1/dd*M_Z^(1/dd-1)*H_IQRsd(Z))^2
    dvars_robust_variance = np.power((1 / transform_power) *
                                     np.power(power_transformed_dvars2_median, (1/transform_power)-1) *
                                     half_iqr_sd_dvars2,
                                     2)
    if debug is True:
        dvars_out_dictionary['DSE_DBG_S'] = np.sqrt(dvars_robust_variance)

    dvars_chi_square_stat = 2 * dvars_robust_mean / dvars_robust_variance * functional_dvars_square_timeseries
    if debug is True:
        dvars_out_dictionary['DSE_DBG_c'] = dvars_chi_square_stat

    dvars_chi_square_df = 2 * np.power(dvars_robust_mean, 2) / dvars_robust_variance
    if debug is True:
        dvars_out_dictionary['DSE_DBG_nu'] = dvars_chi_square_df

    dvars_chi_square_p_values = stats.chi2.sf(dvars_chi_square_stat, dvars_chi_square_df)
    dvars_out_dictionary['DSE_P_Values'] = dvars_chi_square_p_values

    dvars_out_dictionary['DSE_Outliers'] = [True if p < 0.05 and d > 5 else False for p, d in
                                            zip(dvars_chi_square_p_values, functional_delta_pct_dvar)]
    dvars_out_dictionary['DSE_Outlier_Count'] = sum(dvars_out_dictionary['DSE_Outliers'])

    # == Transform to relative DVARS using the method in FSL script DVARS.sh

    # Calculate Robust standard deviation, use the "higher" interpolation method, which seems to be
    # consistent with the way that fslmaths works
    functional_voxel_lower_quartiles = np.percentile(functional_data, 25, axis=1, interpolation="higher")
    functional_voxel_upper_quartiles = np.percentile(functional_data, 75, axis=1, interpolation="higher")
    functional_robust_stdev = (functional_voxel_upper_quartiles - functional_voxel_lower_quartiles) / 1.349
    if debug is True:
        dvars_out_dictionary['DSE_DBG_Mu0'] = functional_robust_stdev.mean()

    # Calculate AR1 coefficients
    functional_ar1_coefficients = np.apply_along_axis(lambda x: ar_nitime(x, center=True), 1, functional_data)

    # Calculate standard deviation of temporal derivative
    functional_robust_stdev_predicted = np.sqrt(2 * (1 - functional_ar1_coefficients)) * functional_robust_stdev

    # Transform
    functional_relative_dvars = np.sqrt(functional_dvars_square_timeseries) / functional_robust_stdev_predicted.mean()
    dvars_out_dictionary['DSE_Relative_DVARS'] = functional_relative_dvars

    return dvars_out_dictionary
