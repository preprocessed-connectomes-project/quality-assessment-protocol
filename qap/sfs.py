
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


def calc_sfs(functional_data_filename, mask_filename, noise_voxel_standard_deviation_percentile=98, mask_erosions=None,
             noise_voxel_mask_filename=None, motion_regressors_filename=None, friston_twentyfour=True,
             detrend_polynomial_order=None, number_of_principle_components=5, debug=False, working_directory=None):
    """ Calculate signal fluctuation sensitivity


    Adaptation of the methods described in:

    DeDora DJ, Nedic S, Katti P, Arnab S, Wald LL, Takahashi A, Van Dijk KRA, Strey HH and Mujica-Parodi LR (2016)
        Signal Fluctuation Sensitivity: An Improved Metric for Optimizing Detection of Resting-State fMRI Networks.
        Front. Neurosci. 10:180. doi: 10.3389/fnins.2016.00180


    Signal Fluctuation Sensitivity is an attempt to estimate the signal to noise ratio of resting state fMRI data. The
    main problem is differentiating signal from noise. In the paper, noise is estimated from the standard deviation of
    a ROI in CSF. Signal is the stand deviation of a ROI in a area of the brain that is being tested. The ratio of
    signal to noise is then multiplied by the ratio of the target region's mean to the noise ROI mean. For this function
    SFS is estimated voxel-by-voxel, so the target regions are each voxel and the noise is estimated from some voxels
    that are determined to be "noise" voxels.

    There are some problems with this method as proposed, chief among them is the using standard deviation of 'signal'
    voxels by standard deviation of noise voxels as the estimator of SNR and second is the selection of noise voxels.
    According to the Behzadi CompCorr paper and several other observations in the literature, voxels that contain
    intensity fluctuations due to physiological noise or head motion tend to have the highest variance. On one hand this
    provides a means of identifying voxels for our noise signal, on the other it means that the ratio of standard
    deviations of signal and noise voxels should be less than one for all voxels that are not the noise voxels, and will
    be one for those. But not all voxels are contaminated by the noise, so it would be more appropriate to only consider
    the amount of noise that is present in the signal voxel. We do this by changing the model to use the ratio of
    standard deviations between the residual of a voxels time series, after the noise has been removed by linear
    regression, and the standard deviation of the fit of the noise to the time series. This can also allow for the
    inclusion of other noise models, such as head motion or RETROICOR regressors.

    The problem of selecting noise voxels can be solved using the tCompCorr method proposed in Behzadi 2007. In this
    method the voxels with the top 2% variance are selected as noise. These voxels tend to be CSF, veins, and voxels on
    the boundary of the brain that are corrupted by head motion. This last set of voxels is troubling because they are
    very likely to contain grey matter as well, and we do not want to include grey matter signal in our 'noise' model.
    We propose two ways of addressing this, one is to erode the whole brain mask so that these boundary voxels are not
    included for consideration and the other is to use a mask of white matter and CSF for the nuisance voxels.

    Once noise voxels are identified, a set of noise regressors are calculated from the top 5 components of a PCA of the
    noise voxel time series, following the recommendations of Behzadi 2007.

    A drawback of these methods of defining noise voxels is that they may miss the noise fluctuations induced by head
    motion, which primarily contaminate  boundary voxels. This may be resolved by including head motion regressors
    in the noise model.

    Behzadi, Y., Restom, K., Liau, J., & Liu, T. T. (2007). A Component Based Noise Correction Method (CompCor) for
        BOLD and Perfusion Based fMRI. NeuroImage, 37(1), 90–101. http://doi.org/10.1016/j.neuroimage.2007.04.042

    :type functional_data_filename: str
    :param functional_data_filename: The filepath to the NIFTI file containing the functional
                       time series.
    :type mask_filename: str
    :param mask_filename: The file path to the NIFTI file containing the binary
                      functional brain mask.

    :type noise_voxel_standard_deviation_percentile: int
    :param noise_voxel_standard_deviation_percentile: The default algorithm for selecting noise voxels is to choose the
         voxels whose standard deviation in a top fraction of all voxels' standard deviations. This value determines
         the top fraction that should be used. Expressed as an integer between 0 and 100, where 95 (for example)
         indicates the top 5% of voxels (the 95th percentile), etc...

    :type mask_erosions: int
    :param mask_erosions: The number of times to erode the mask before using it to choose noise voxels by standard
        deviation. This is used to avoid including grey matter signal from voxels at the boundary of the brain in the
        estimate noise signal. 3 to 5 should be OK. If None, no erosion will be performed.


    :type noise_voxel_mask_filename: str
    :param noise_voxel_mask_filename: String containing the full path of the nifti file to use as mask for defining
         the location of the noise voxels. If not included then will use the voxels whose variance is in the top 2%
         of all voxel variances.

    :type motion_regressors_filename: str
    :param motion_regressors_filename: string containing the full path to a tab seperated values file containing
         motion regressors to include in the noise model.

    :type friston_twentyfour: bool
    :param friston_twentyfour: flag indicating whether or not Friston's 24 parameter motion model should be used for
         the motion regressor, this is only relevant if a motion_regressors_filename is provided

    :type detrend_polynomial_order: int
    :param detrend_polynomial_order: integer representing the order of polynomial to include in the noise model. This is
         typically used to correct for fluctuations due to scanner drift and slow head motion. 2 should be good. If
         None, polynomial detrending will not be included.

    :type number_of_principle_components: int
    :param number_of_principle_components: the number of principle components to extract from the noise voxel time
          series

    :type debug: bool
    :param debug: a flag indicating whether extra values should be added to the output dictionary
                      to assist with debugging

    :type working_directory:str
    :param working_directory: path to write output files to, if not set, will write to current directory

    :rtype: dictionary
    :return: dictionary containing values calculated from time course data.

        'SFS_STD_File_Path' : Path to nifti file containing voxelwise standard deviation
        'SFS_Noise_Mask_File_Path' : Path to nifti file containing mask used for selecting noise voxels
        'SFS_Nose_Voxels_File_Path' : Path to nifti file with noise voxels indicated with a '1' and '0' otherwise
        'SFS_File_Path': Path to nifti file containing voxelwise SFS
        'SFS_Mean': The mean SFS across the brain voxels.
        'SFS_Std' : The standard deviation of SFS across brain voxels.
        'SFS_Median': The median SFS across the brain voxels.
        'SFS_25pct' : The 25th percentile of SFS.
        'SFS_75pct' : The 75th percentile of SFS.

        If debug is True the following parameters will also be included in the dictionary. These variables are named
        to make them consistent with the variables output by the DSE Matlab toolbox
        'SFS_DBG_Dim': The dimensions of the functional imaging data after it has been conformed (# voxels x # frames)
        'SFS_DBG_Dim_Masked': The dimensions of the functional imaging data after it has been conformed and masked (#
            voxels x # frames)
        'SFS_DBG_GrandMean0': The grand mean of the functional imaging data prior to scaling
        'SFS_DBG_Scale2Norm': The scalar value that was used in scaling the functional imaging data
        'SFS_DBG_GrandMean': The grand mean of the functional imaging data after scaling
        'SFS_DBG_Std_Thresh': The standard deviation at the threshold percentile.


    """

    import numpy as np
    import nibabel as nb
    from scipy.linalg import svd
    import scipy.ndimage.morphology as morph
    import os

    if working_directory is 'None':
        working_directory = os.getcwd()

    sfs_out_dictionary = {}

    # load and check functional data and conform it to our needs (float, 2D vx x time array)
    functional_data_image = nb.load(functional_data_filename)

    if len(functional_data_image.shape) != 4:
        raise ValueError("SFS expects the input functional data to be 4D.")

    functional_data = functional_data_image.get_data().astype(np.float).reshape(
        np.prod(functional_data_image.shape[0:3]), functional_data_image.shape[3])

    # load, check, and conform the mask
    mask_image = nb.load(mask_filename)

    if len(mask_image.shape) != 3:
        raise ValueError("SFS expects the input mask to be 3D.")

    if mask_image.shape != functional_data_image.shape[:3]:
        raise ValueError("Mask image shape ({0}) does not match "
                         "functional data shape ({1}).".format(mask_image.shape,
                                                               functional_data_image.shape))

    mask_data = mask_image.get_data().reshape(np.prod(mask_image.shape))
    if debug is True:
        sfs_out_dictionary["SFS_DBG_Dim"] = functional_data.shape

    # binarize the mask, allows us to work with non-binary input images
    # assume that nonzero is in -brain and zero is out of brain
    mask_data[np.isclose(mask_data, 0)] = 0
    mask_data[mask_data != 0] = 1

    if debug is True:
        sfs_out_dictionary["SFS_DBG_Mask_Voxel_Count"] = mask_data.sum()

    # update mask to exclude voxels that have zero variance in the functional data
    mask_data[np.isclose(functional_data.var(1), 0)] = 0

    if debug is True:
        sfs_out_dictionary["SFS_DBG_Variance_Filtered_Mask_Voxel_Count"] = mask_data.sum()

    # mask functional data
    functional_data = functional_data[mask_data == 1, :]
    if debug is True:
        sfs_out_dictionary["SFS_DBG_Dim_Masked"] = functional_data.shape

    # scale values so that grand mean == 100 (keeping consistent with the DSE Matlab toolbox)
    functional_data_scale = 100 / functional_data.mean()
    if debug is True:
        sfs_out_dictionary["DSE_DBG_GrandMean0"] = functional_data.mean()
        sfs_out_dictionary["DSE_DBG_Scale2Norm"] = 100 / functional_data.mean()
    functional_data = functional_data_scale * functional_data

    if debug is True:
        sfs_out_dictionary["DSE_DBG_GrandMean"] = functional_data.mean()

    # save the mean to use in SFS calculations
    functional_data_mean = functional_data.mean(1)

    if debug is True:
        mean_output_filename = os.path.join(working_directory, "sfs_mean.nii.gz")
        mean_output_data = mask_data.copy()
        mean_output_data[mask_data == 1] = functional_data_mean
        mean_output_image = nb.Nifti1Image(mean_output_data.reshape(mask_image.shape),  mask_image.affine)
        mean_output_image.to_filename(mean_output_filename)
        sfs_out_dictionary["SFS_Mean_Image_Filename"] = mean_output_filename

    # center the data by subtracting the mean
    functional_data = np.apply_along_axis(lambda x: scale_then_center(x), 1, functional_data)

    # now select noise voxels
    if noise_voxel_mask_filename:
        noise_voxel_mask_image = nb.load(noise_voxel_mask_filename)

        if len(noise_voxel_mask_image.shape) != 3:
            raise ValueError("SFS expects the noise voxel mask to be 3D.")

        if noise_voxel_mask_image.shape != functional_data_image.shape[0:3]:
            raise ValueError("Noise voxel mask shape ({0}) does not match "
                             "functional data shape ({1}).".format(mask_image.shape,
                                                                   functional_data_image.shape))

        # conform and binarize the mask
        noise_voxel_mask = noise_voxel_mask_image.get_data().reshape(np.prod(noise_voxel_mask_image.shape))
        noise_voxel_mask[np.isclose(noise_voxel_mask, 0)] = 0
        noise_voxel_mask[noise_voxel_mask != 0] = 1

    else:

        if mask_erosions:
            if not isinstance(mask_erosions, int) and mask_erosions <= 0:
                raise ValueError("Mask erosions ({0}, {1}) should be an "
                                 "integer greater than 0".format(mask_erosions, type(mask_erosions)))

            noise_voxel_mask = morph.binary_erosion(mask_data.reshape(mask_image.shape),
                                                    iterations=mask_erosions).reshape(mask_data.shape)
        else:
            noise_voxel_mask = mask_data.copy()

        if debug is True:
            noise_voxel_mask_output_filename = os.path.join(working_directory,
                                                            "sfs_noise_voxel_mask_pre_threshold.nii.gz")
            noise_voxel_mask_output_image = nb.Nifti1Image(noise_voxel_mask.reshape(mask_image.shape).astype(np.uint8),
                                                           mask_image.affine)
            noise_voxel_mask_output_image.to_filename(noise_voxel_mask_output_filename)
            sfs_out_dictionary["SFS_Noise_Voxel_Mask_Pre_Thresh_Filename"] = noise_voxel_mask_output_filename

        noise_voxels = noise_voxel_mask[mask_data == 1]

        if debug is True:
            sfs_out_dictionary["SFS_DBG_Eroded_Mask_Voxel_Count"] = noise_voxel_mask.sum()

        # calculate the voxels with standard deviation in the top percentile of all voxels
        if not noise_voxel_standard_deviation_percentile or \
                not isinstance(noise_voxel_standard_deviation_percentile, int) or \
                noise_voxel_standard_deviation_percentile < 0 or \
                noise_voxel_standard_deviation_percentile > 100:

            raise ValueError(
                "Noise voxel standard deviation percentile ({0}, {1}) should be an integer between 0 and 100".format(
                    noise_voxel_standard_deviation_percentile, type(noise_voxel_standard_deviation_percentile)))

        functional_data_standard_deviation = functional_data.std(1)

        if debug is True:
            standard_deviation_output_filename = os.path.join(working_directory, "sfs_standard_deviation.nii.gz")
            standard_deviation_output_data = mask_data.copy()
            standard_deviation_output_data[mask_data == 1] = functional_data_standard_deviation

            standard_deviation_output_image = nb.Nifti1Image(standard_deviation_output_data.reshape(mask_image.shape),
                                                             mask_image.affine)

            standard_deviation_output_image.to_filename(standard_deviation_output_filename)
            sfs_out_dictionary["SFS_Standard_Deviation_Image_Filename"] = standard_deviation_output_filename

        critical_value = np.percentile(functional_data_standard_deviation[noise_voxels == 1],
                                       noise_voxel_standard_deviation_percentile)
        if debug is True:
            sfs_out_dictionary["SFS_DBG_Critical_Value"] = critical_value

        noise_voxels = noise_voxels * (functional_data_standard_deviation >= critical_value)

        noise_voxel_mask[mask_data == 1] = noise_voxels.astype(np.uint8)

    noise_voxel_mask_output_filename = os.path.join(working_directory, "sfs_noise_voxel_mask.nii.gz")
    noise_voxel_mask_output_image = nb.Nifti1Image(noise_voxel_mask.reshape(mask_image.shape).astype(np.uint8),
                                                   mask_image.affine)
    noise_voxel_mask_output_image.to_filename(noise_voxel_mask_output_filename)
    sfs_out_dictionary["SFS_Noise_Voxel_Mask_Filename"] = noise_voxel_mask_output_filename

    if debug is True:
        sfs_out_dictionary["SFS_DBG_Noise_Voxel_Count"] = noise_voxel_mask.sum()

    noise_voxels = noise_voxel_mask[mask_data == 1]

    # now we do the regressions that will be used to calculate SFS.

    # if the user specified polynomial detrending and motion regressors, lets remove those first so that they will
    # not contaminate the nuisance regressor. The nuisance regressors will then be regressed from the functional data.
    # At worst case this two step procedure will would be exactly the same as entering everything in a single model.
    # At best it will reduce the collinearity in the design matrix that might make things easier to estimate.

    # noise_time_series = functional_data[noise_voxel_mask[mask_data == 1] == 1]
    nuisance_regressors = None
    if motion_regressors_filename or detrend_polynomial_order:

        if detrend_polynomial_order:
            if not isinstance(detrend_polynomial_order,
                              int) or detrend_polynomial_order < 0 or detrend_polynomial_order > 5:
                raise ValueError('Detrend polynomial order ({0}, {1}) should be an integer in the range [0,10)'.format(
                    detrend_polynomial_order, type(detrend_polynomial_order)))

            # includes 0th order for the mean
            nuisance_regressors = np.ones((functional_data.shape[1], detrend_polynomial_order+1))
            for order in range(0, detrend_polynomial_order+1):
                nuisance_regressors[:, order] = np.arange(0, functional_data.shape[1], dtype=np.float32) ** (order)

        if motion_regressors_filename:
            temp_regressors = np.loadtxt(motion_regressors_filename)
            if friston_twentyfour is True:
                # squares
                temp_regressors = np.append(temp_regressors, temp_regressors ** 2, axis=1)
                # zero padded one-lags of std motion and squares
                temp_regressors = np.append(temp_regressors,
                                            np.append(np.zeros((1, temp_regressors.shape[1]), dtype=np.float64),
                                                      temp_regressors[1:, :], axis=0), axis=1)

            if nuisance_regressors is not None:
                nuisance_regressors = np.append(nuisance_regressors, temp_regressors, axis=1)
            else:
                nuisance_regressors = temp_regressors

        # perform the OLS fit
        nuisance_betas = np.linalg.inv(nuisance_regressors.transpose().dot(nuisance_regressors)).dot(
            nuisance_regressors.transpose()).dot(functional_data.transpose())

        functional_residuals = functional_data - nuisance_regressors.dot(nuisance_betas).transpose()

        # get noise voxels
        noise_voxel_time_series = functional_residuals[noise_voxels == 1, :]

    else:
        noise_voxel_time_series = functional_data[noise_voxels == 1, :]

    # now get the top 5 principle components of the noise_voxel_time_series
    (left_principle_components, singular_values, right_principle_components_transposed) = svd(noise_voxel_time_series,
                                                                                            full_matrices=False)

    noise_components = right_principle_components_transposed[0:number_of_principle_components, :].transpose()

    if nuisance_regressors is not None:
        # if there are nuisance regressors add them back in to the design matrix
        nuisance_regressors = np.append(nuisance_regressors, noise_components, axis=1)
    else:
        nuisance_regressors = noise_components

    if debug:
        sfs_out_dictionary["SFS_Nuisance_Regressors"] = nuisance_regressors.transpose()

    # redo OLS regression, this time the results will be used to calculate SFS
    nuisance_betas = np.linalg.inv(nuisance_regressors.transpose().dot(nuisance_regressors)).dot(
        nuisance_regressors.transpose()).dot(functional_data.transpose())

    modeled_var = nuisance_regressors.dot(nuisance_betas).transpose().var(1)
    residual_var = (functional_data - nuisance_regressors.dot(nuisance_betas).transpose()).var(1)
    noise_voxel_mean = functional_data_mean[noise_voxels == 1].mean()

    signal_fluctuation_sensitivity = functional_data_mean / noise_voxel_mean * residual_var / modeled_var

    signal_fluctuation_sensitivity_output_filename = os.path.join(working_directory, "sfs.nii.gz")
    signal_fluctuation_sensitivity_output_data = mask_data.copy()
    signal_fluctuation_sensitivity_output_data[mask_data == 1] = signal_fluctuation_sensitivity

    signal_fluctuation_sensitivity_output_image = nb.Nifti1Image(
        signal_fluctuation_sensitivity_output_data.reshape(mask_image.shape),
        mask_image.affine)

    signal_fluctuation_sensitivity_output_image.to_filename(signal_fluctuation_sensitivity_output_filename)
    sfs_out_dictionary["SFS_SImage_Filename"] = signal_fluctuation_sensitivity_output_filename

    sfs_out_dictionary["SFS_Mean"] = np.mean(signal_fluctuation_sensitivity)
    sfs_out_dictionary["SFS_Std"] = np.std(signal_fluctuation_sensitivity)
    sfs_out_dictionary["SFS_median"] = np.median(signal_fluctuation_sensitivity)
    sfs_out_dictionary["SFS_25pct"] = np.percentile(signal_fluctuation_sensitivity, 25)
    sfs_out_dictionary["SFS_75pct"] = np.percentile(signal_fluctuation_sensitivity, 75)

    return sfs_out_dictionary
