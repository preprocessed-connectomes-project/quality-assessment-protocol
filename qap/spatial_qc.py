
def summary_mask(anat_data, mask_data):
    """Will calculate the three values (mean, stdev, and size) and return them
    as a tuple.

    :type anat_data: NumPy array
    :param anat_data: The anatomical scan data.
    :type mask_data: NumPy array
    :param mask_data: The binary mask to mask the anatomical data with.
    :rtype: tuple
    :return: The summary values (mean, standard deviation, size) of the scan.
    """
    
    import numpy as np
    
    anat_masked = anat_data[mask_data == 1].astype(np.float)
    mean = anat_masked.mean()
    std = anat_masked.std(ddof=1)
    size = len(anat_masked)
    
    return anat_masked, mean, std, size


def check_datatype(background):
    """Process the image data to only include non-negative integer values.

    :type background: NumPy array
    :param background: The voxel values of teh background (outside of the head
                       ) of the anatomical image.
    :rtype: NumPy array
    :return: The input array with floats converted to integers and
             negative values set to zero.
    """

    import numpy as np

    # If this is float then downgrade the data to an integer with some checks
    if np.issubdtype(background.dtype, float):
        background2 = background.astype('int32')    
        # Ensure downgrading datatype didn't really change the values
        if np.abs(background2.astype('float32') - background).mean() > 0.05:
            print("WARNING: Downgraded float to an int but values are different by more than 0.05")
        background = background2
        del background2
   
    # We only allow integer values for now
    if not np.issubdtype(background.dtype, int):
        print("QI1 can not be calculated for data that is not integer or floating point: {0}".format(background.dtype))
        raise TypeError
        
    # convert any negative voxel values to zero, provide warning
    for vox in background.flatten():
        if vox < 0:
            print("\nWARNING: Negative voxel values in anatomical scan converted to zero.\n")
            background = background.clip(0)
            break

    return background


def convert_negatives(img_data):
    """Convert any negative voxel values to zero and provide a warning.

    :type img_data: NumPy array
    :param img_data: The anatomical image's voxel values.
    :rtype: NumPy array
    :return: The input array with negative values set to zero.
    """

    for vox in img_data.flatten():
        if vox < 0:
            print("\nWARNING: Negative voxel values in anatomical scan converted to zero.\n")
            img_data = img_data.clip(0)
            break

    return img_data


def snr(mean_fg, std_bg):
    """Calculate the Signal-to-Noise Ratio (SNR) of an image.

    - For anatomical images:
        SNR = (mean GM intensity) / (std of background intensities)
    - For functional images:
        SNR = (mean brain intensity) / (std of background intensities)

    :type mean_fg: float
    :param mean_fg: The mean value of voxel intensities in the foreground
                    (either within the head or a particular tissue) of the
                    image.
    :type std_bg: float
    :param std_bg: The standard deviation of the voxel intensities of the
                   background (outside of the head) voxels.
    :rtype: float
    :return: The signal-to-noise ratio (SNR).
    """

    snr_val = mean_fg / std_bg

    return snr_val


def cnr(mean_gm, mean_wm, std_bg):
    """Calculate Contrast-to-Noise Ratio (CNR) of an image.

    - CNR = |(mean GM intensity) - (mean WM intensity)| /
                                          (std of background intensities)

    :type mean_gm: float
    :param mean_gm: The mean value of the gray matter voxels.
    :type mean_wm: float
    :param mean_wm: The mean value of the whiet matter voxels.
    :type std_bg: float
    :param std_bg: The standard deviation of the voxel intensities of the
                   background (outside the head) voxels.
    :rtype: float
    :return: The contrast-to-noise (CNR) ratio.
    """

    import numpy as np
    cnr_val = np.abs(mean_gm - mean_wm)/std_bg

    return cnr_val


def cortical_contrast(mean_gm, mean_wm):
    """Calculate the vertex-wise cortical contrast.

    - cortical contrast = (mean WM intensity) - (mean GM intensity) /
                            ( (mean WM intensity + mean GM intensity) / 2 )

    :type mean_gm: float
    :param mean_gm: The mean value of the gray matter voxels.
    :type mean_wm: float
    :param mean_wm: The mean value of the white matter voxels.
    :rtype: float
    :return: The cortical contrast value.
    """

    cort_con = (mean_wm - mean_gm) / ((mean_wm + mean_gm) / 2)

    return cort_con

    
def calc_fber(anat_data, skull_mask_data, bg_mask_data):
    """Calculate the Foreground-to-Background Energy Ratio (FBER) of an image.

    - FBER = (mean foreground energy) / (mean background energy)

    :type anat_data: NumPy array
    :param anat_data: The anatomical/spatial data of the image.
    :type skull_mask_data: NumPy array
    :param skull_mask_data: The binary mask defining the head.
    :type bg_mask_data: NumPy array
    :param bg_mask_data: The binary mask defining the background (outside of
                         the head).
    :rtype: float
    :return: The foreground-to-background energy ratio (FBER).
    """

    import numpy as np

    mean_fg = (np.abs(anat_data[skull_mask_data == 1]) ** 2).sum() / (skull_mask_data.sum())
    mean_bg = (np.abs(anat_data[bg_mask_data == 1]) ** 2).sum() / (bg_mask_data.size - bg_mask_data.sum())
    fber = mean_fg / mean_bg

    return fber


def calc_efc(anat_data):
    """Calculate the Entropy Focus Criterion of the image.

    - EFC based on Atkinson 1997, IEEE TMI
    - We normalize the original equation by the maximum entropy so our EFC
      can be easily compared across images with different dimensions.

    :type anat_data: Nibabel data
    :param anat_data: The anatomical image data.
    :rtype: float
    :return: The entropy focus criterion (EFC) value.
    """

    import numpy as np
        
    # let's get rid of those negative values
    anat_data = convert_negatives(anat_data)
        
    # Calculate the maximum value of the EFC (which occurs any time all 
    # voxels have the same value)
    efc_max = 1.0 * np.prod(anat_data.shape) * (1.0 / np.sqrt(np.prod(anat_data.shape))) * np.log(
        1.0 / np.sqrt(np.prod(anat_data.shape)))
    
    # Calculate the total image energy
    b_max = np.sqrt((anat_data**2).sum())
    
    # Calculate EFC (add 1e-16 to the image data to keep log happy)
    efc = (1.0 / efc_max) * np.sum((anat_data / b_max) * np.log((anat_data + 1e-16) / b_max))
    
    if np.isnan(efc): 
        print("NaN found for efc ({0:3.2f},{1:3.2f})".format(efc_max, b_max))
    
    return efc


def artifacts(anatomical_reorient, qap_head_mask_path, exclude_zeroes=False):
    """Create a binary mask of the anatomical background voxels corrupted by
    artifacts via a binary opening operation.

    :param anatomical_reorient: String filepath to the deobliqued, reoriented
                                anatomical scan NIFTI file.
    :param qap_head_mask_path: String filepath to the binary head mask NIFTI
                               file, including the region in front of and
                               below the mouth and nose.
    :param exclude_zeroes: Boolean choosing whether to exclude artificial
                           zero voxels.
    :return: String filepath to the background artifacts NIFTI file.
    """

    import os
    import numpy as np
    import nibabel as nb
    import scipy.ndimage as nd
    from qap.qap_utils import load_image, load_mask, \
        create_anatomical_background_mask

    # Load the data
    anat_data, anat_img = load_image(anatomical_reorient, return_img=True)
    fg_mask = load_mask(qap_head_mask_path, anatomical_reorient)

    # bg_mask is the inversion of the "qap_head_mask"
    bg_mask = create_anatomical_background_mask(anat_data, fg_mask,
                                                exclude_zeroes)

    # Create an image containing only background voxels (everything 
    # outside bg_mask set to 0)
    background = anat_data.copy()
    background[bg_mask != 1] = 0

    # make sure the datatype is an int
    background = check_datatype(background)

    # Find the background threshold (the most frequently occurring value 
    # excluding 0)
    bg_counts = np.bincount(background.flatten())
    bg_threshold = np.argmax(bg_counts[1:]) + 1

    # Apply this threshold to the background voxels to identify voxels
    # contributing artifacts. 
    background[background <= bg_threshold] = 0
    background[background != 0] = 1

    # Create a structural element to be used in an opening operation.
    structure_element = np.zeros((3, 3, 3))
    structure_element[0, 1, 1] = 1
    structure_element[1, 1, :] = 1
    structure_element[1, :, 1] = 1
    structure_element[2, 1, 1] = 1

    # Perform an opening operation on the background data.
    background = nd.binary_opening(background, structure=structure_element)

    bg_mask_img = nb.Nifti1Image(bg_mask, anat_img.affine, anat_img.header)
    # this writes it into the node's folder in the working directory
    bg_mask_file = os.path.join(os.getcwd(), "qap-bg-head-mask.nii.gz")
    bg_mask_img.to_filename(bg_mask_file)

    fav_bg_img = nb.Nifti1Image(background.astype(int), anat_img.affine, anat_img.header)
    # this writes it into the node's folder in the working directory
    fav_bg_file = os.path.join(os.getcwd(), "fav-artifacts-background.nii.gz")
    fav_bg_img.to_filename(fav_bg_file)

    return fav_bg_file, bg_mask_file


def fav(fav_artifacts_bg, anatomical_reorient, bg_head_mask, calculate_quality_index_2=False):
    """Calculate FAV, the fraction of total voxels that contain artifacts.

    - Detect artifacts in the anatomical image using the method described in
      Mortamet et al. 2009 (MRM).
    - Optionally, also calculates QI2, the distance between the distribution
      of noise voxel (non-artifact background voxels) intensities, and a
      Ricean distribution.

    :param fav_artifacts_bg: String filepath to the mask NIFTI file outlining
                             the background voxels containing artifacts.
    :param anatomical_reorient: String filepath to the anatomical scan NIFTI
                                file.
    :param bg_head_mask: String filepath to the binary background mask NIFTI
                         file.
    :param calculate_quality_index_2: Boolean switch determining whether to calculate Qi2.
    :return: Float values of Qi1 and Qi2 (if Qi2 not calculated, return None).
    """

    import numpy as np
    import nibabel
    from qap.qap_utils import load_image

    # Load the data and make sure that it is 0 or 1
    background = nibabel.load(fav_artifacts_bg).get_data()
    background[np.isclose(background, 0)] = 0
    background[background != 0] = 1

    bg_mask = nibabel.load(bg_head_mask).get_data()
    bg_mask[np.isclose(bg_mask, 0)] = 0
    bg_mask[bg_mask != 0] = 1

    # Count the number of voxels that remain after the opening operation. 
    # These are artifacts.
    quality_index_1 = float(background.sum()) / float(bg_mask.sum())

    if calculate_quality_index_2:
        import scipy.stats as ss
        # Load the data
        anat_data = load_image(anatomical_reorient, return_img=False)

        background_artifacts = load_image(fav_artifacts_bg, return_img=False)

        # Now lets focus on the noise, which is everything in the background
        # that was not identified as artifact
        print("{0} background noise voxels {1}".format((bg_mask - background_artifacts).sum(),
              (bg_mask - background_artifacts).sum() / float(bg_mask.sum())))
        background_noise = anat_data[(bg_mask - background_artifacts) == 1]
        print("{0}".format(np.sum(background_noise > 0)))

        # calculate the histogram of the noise and its derivative
        histogram = np.histogram(background_noise, bins='auto')[0]
        histogram = 1.0 * histogram / float(histogram.sum())
        delta_histogram = histogram[1:] - histogram[:-1]

        # find the first value on the right tail, i.e. tail with negative
        # slope, i.e. dH < 0 that is less than or equal to half of the
        # histograms max
        first_negative_slope = np.nonzero(delta_histogram < 0)[0][0]
        half_max_right_tail = np.nonzero(histogram[first_negative_slope:] < (histogram.max()/2))[0][0]

        # divide by the standard deviation
        background_noise_z_statistic = background_noise / background_noise.std()
        background_chi_params = ss.chi.fit(background_noise_z_statistic)
        # print background_chi_params

        # now generate values that are consistent with the histogram
        yx = [float(y) / background_noise.std() for y in range(0, histogram.size)]
        rvs = ss.chi.pdf(yx, background_chi_params[0], loc=background_chi_params[1], scale=background_chi_params[2])

        # now we can calculate the goodness of fit
        gof = np.average(np.absolute(histogram[half_max_right_tail:] - rvs[half_max_right_tail:]))
        quality_index_2 = quality_index_1 + gof

    else:
        quality_index_2 = None

    return quality_index_1, quality_index_2


def fwhm(anat_file, mask_file, out_vox=False):
    """Calculate the FWHM of the input image using AFNI's 3dFWHMx.

    - Uses AFNI 3dFWHMx. More details here:
        https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dFWHMx.html

    :type anat_file: str
    :param anat_file: The filepath to the anatomical image NIFTI file.
    :type mask_file: str
    :param mask_file: The filepath to the binary head mask NIFTI file.
    :type out_vox: bool
    :param out_vox: (default: False) Output the FWHM as number of voxels
                    instead of mm (the default).
    :rtype: tuple
    :return: A tuple of the FWHM values (x, y, z, and combined).
    """

    import nibabel as nib
    import numpy as np
    from scipy.special import cbrt

    import subprocess

    fwhm_string_list = ["3dFWHMx", "-combined", "-mask", mask_file,
                        "-input", anat_file]
    try:
        retcode = subprocess.check_output(fwhm_string_list)
    except Exception as e:
        err = "\n\n[!] Something went wrong with AFNI's 3dFWHMx. Error " \
              "details: %s\n\n" % e
        raise Exception(err)

    # extract output
    vals = np.array(retcode.split()[-4:], dtype=np.float)
    
    if out_vox:
        # get pixel dimensions
        img = nib.load(anat_file)
        hdr = img.get_header()
        pixdim = hdr['pixdim'][1:4]
    
        # convert to voxels
        pixdim = np.append(pixdim, cbrt(pixdim.prod()))
        # get the geometrix mean
        vals = vals / pixdim
    
    return tuple(vals)


def skew_and_kurt(anat_data):

    from scipy import stats

    skew = stats.skew(anat_data)
    kurt = stats.kurtosis(anat_data)

    return skew, kurt


def ghost_direction(epi_data, mask_data, direction="y", ref_file=None,
                    out_file=None):
    """Calculate the Ghost to Signal Ratio of EPI images.

    - GSR from Giannelli 2010. More details here:
        https://www.ncbi.nlm.nih.gov/pubmed/21081879
    - This should be used for EPI images where the phase encoding direction
      is known.

    :type epi_data: Nibabel data
    :param epi_data: The mean of the functional timeseries.
    :type mask_data: Nibabel data
    :param mask_data: The functional brain binary mask data.
    :type direction: str
    :param direction: (default: 'y') The phase encoding direction of the EPI
                      image.
    :type ref_file: str
    :param ref_file: (default: None) If you are saving the Nyquist ghost mask,
                      this is the filepath of the reference file to use to
                     populate the header of the ghost mask NIFTI file.
    :type out_file: str
    :param out_file: (default: None) If you are saving the Nyquist ghost mask,
                      this is the filepath to the ghost mask NIFTI file.
    :rtype: float
    :return: The ghost-to-signal ratio (GSR) value.
    """
    
    import numpy as np
    
    # first we need to make a nyquist ghost mask, we do this by circle 
    # shifting the original mask by N/2 and then removing the intersection
    # with the original mask
    n2_mask_data = np.zeros_like(mask_data)
    
    # rotate by n/2
    if direction == "x":
        n2 = int(np.floor(mask_data.shape[0]/2))
        n2_mask_data[:n2, :, :] = mask_data[n2:(n2*2), :, :]
        n2_mask_data[n2:(n2*2), :, :] = mask_data[:n2, :, :]
    elif direction == "y":
        n2 = int(np.floor(mask_data.shape[1]/2))
        n2_mask_data[:, :n2, :] = mask_data[:, n2:(n2*2), :]
        n2_mask_data[:, n2:(n2*2), :] = mask_data[:, :n2, :]
    elif direction == "z":
        n2 = int(np.floor(mask_data.shape[2]/2))
        n2_mask_data[:, :, :n2] = mask_data[:, :, n2:(n2*2)]
        n2_mask_data[:, :, n2:(n2*2)] = mask_data[:, :, :n2]
    else:
        raise ValueError("Unknow direction {0}, should be x, y, or z".format(direction))

    # now remove the intersection with the original mask
    n2_mask_data = n2_mask_data * (1 - mask_data)
    
    # now create a non-ghost background region, that contains 2s
    n2_mask_data = n2_mask_data + 2 * (1 - n2_mask_data - mask_data)
    
    # Save mask
    if ref_file is not None and out_file is not None:
        import nibabel as nib
        ref = nib.load(ref_file)
        out = nib.Nifti1Image(n2_mask_data, ref.get_affine(), ref.get_header()) 
        out.to_filename(out_file)
   
    # now we calculate the Ghost to signal ratio, but here we define signal
    # as the entire foreground image
    gsr = (epi_data[n2_mask_data == 1].mean() - epi_data[n2_mask_data == 2].mean())/epi_data[n2_mask_data == 0].mean()

    return gsr


def ghost_all(epi_data, mask_data):
    """Call the 'ghost_direction' function on all possible phase encoding 
    directions.

    :type epi_data: Nibabel data
    :param epi_data: The mean of the functional timeseries.
    :type mask_data: Nibabel data
    :param mask_data: The functional brain binary mask data.
    :rtype: tuple
    :return: The ghost-to-signal ratios (GSR) of each phase encoding direction.
    """
    import numpy as np
    
    directions = ["x", "y"]
    ghost_signal_ratios = [ghost_direction(epi_data, mask_data, d) for d in directions] + [np.NaN]
    
    return tuple(ghost_signal_ratios)
