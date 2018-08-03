def calc_global_correlation(functional_data_filename, mask_filename, debug=False):
    """ Calculate the global correlation (GCOR) of the functional timeseries.

    - From "Correcting Brain-Wide Correlation Differences in Resting-State
      fMRI", Ziad S. Saad et al. More info here:
        https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3749702

    :type functional_data_filename: str
    :param functional_data_filename: Full path to the functional data to use in the calculation

    :type mask_filename: str
    :param mask_filename: full path to the mask to use to distinguish in-brain and out-of-brain voxels
        in the calculation

    :param debug: whether debug information should be written to the output dictionary

    :rtype: dictionary
    :return: dictionary containing values calculated from time course data.

        'GCOR': global correlation value

        If debug is True the following parameters will also be included in the dictionary.
        'GCOR_DBG_Dim': The dimensions of the functional imaging data after it has been conformed (# voxels x # frames)
        'GCOR_DBG_Dim_Masked': The dimensions of the functional imaging data after it has been conformed and masked (#
            voxels x # frames)
        'GCOR_DBG_Mask_Voxel_Count': The number of voxels included in the brain mask
        'GCOR_DBG_Variance_Filtered_Mask_Voxel_Count': number of voxels after zero variance voxels have been removed
        'GCOR_DBG_VoxelCorrelations': vector containing the correlation between each voxel's time series and the mean
            time series

    """
    from qap.temporal_qc import variance_scale_then_mean_center
    import numpy as np
    import nibabel as nb

    gcor_out_dictionary = {}

    # load and check functional data and conform it to our needs (float, 2D vx x time array)
    functional_data_image = nb.load(functional_data_filename)

    if len(functional_data_image.shape) != 4:
        raise ValueError("GCOR expects the input functional data to be 4D.")

    functional_data = functional_data_image.get_data().astype(np.float).reshape(
        np.prod(functional_data_image.shape[0:3]), functional_data_image.shape[3])

    # load, check, and conform the mask
    mask_image = nb.load(mask_filename)

    if len(mask_image.shape) != 3:
        raise ValueError("GCOR expects the input mask to be 3D.")

    if mask_image.shape != functional_data_image.shape[:3]:
        raise ValueError("Mask image shape ({0}) does not match "
                         "functional data shape ({1}).".format(mask_image.shape,
                                                               functional_data_image.shape))

    mask_data = mask_image.get_data().reshape(np.prod(mask_image.shape))
    if debug is True:
        gcor_out_dictionary["GCOR_DBG_Dim"] = functional_data.shape

    # binarize the mask, allows us to work with non-binary input images
    # assume that nonzero is in -brain and zero is out of brain
    mask_data[np.isclose(mask_data, 0)] = 0
    mask_data[mask_data != 0] = 1

    if debug is True:
        gcor_out_dictionary["GCOR_DBG_Mask_Voxel_Count"] = mask_data.sum()

    # update mask to exclude voxels that have zero variance in the functional data
    mask_data[np.isclose(functional_data.var(1), 0)] = 0

    if debug is True:
        gcor_out_dictionary["GCOR_DBG_Variance_Filtered_Mask_Voxel_Count"] = mask_data.sum()

    # mask functional data
    functional_data = functional_data[mask_data == 1, :]
    if debug is True:
        gcor_out_dictionary["GCOR_DBG_Dim_Masked"] = functional_data.shape

    # variance normalize and mean center the data
    functional_data = np.apply_along_axis(func1d=variance_scale_then_mean_center, axis=1, arr=functional_data,
                                          var_scale=True, mean_center=True)

    # calculate the mean vector
    functional_data_mean_time_series = functional_data.mean(0)

    voxel_correlations = functional_data.dot(functional_data_mean_time_series) / functional_data.shape[1]

    if debug is True:
        gcor_out_dictionary['GCOR_DBG_voxel_correlations'] = voxel_correlations

    gcor_out_dictionary['GCOR'] = voxel_correlations.mean()

    return gcor_out_dictionary
