
def load_image(image_file):
    """Load a raw scan image from a NIFTI file and check it.

    Keyword arguments:
      image_file -- [string] path to the image, usually a structural or 
                    functional scan

    Returns:
      dat -- [Nibabel data] image data in Nibabel format
    """

    import nibabel as nib
    import numpy as np

    from workflow_utils import raise_smart_exception

    try:
        img = nib.load(image_file)
    except:
        raise_smart_exception(locals())

    dat = img.get_data()

    # Ensure that data is cast as at least 32-bit
    if np.issubdtype(dat.dtype, float):
        dat = dat.astype('float32')
        # Check for negative values
        if (dat < 0).any():
            print "found negative values, setting to zero (see file: %s)" \
                  % image_file
            dat[dat<0] = 0

    elif np.issubdtype(dat.dtype, int):
        dat = dat.astype('int32')

    elif np.issubdtype(dat.dtype, np.uint8):
        dat = dat.astype(np.uint8)

    else:
        msg = "Error: Unknown datatype %s" % dat.dtype
        raise_smart_exception(locals(),msg)

    return dat


def load_mask(mask_file, ref_file):
    """Load a mask from a NIFTI file and check the shape and dimensions.

    Keyword arguments:
      mask_file -- [string] filepath to binarized mask file
      ref_file -- [string] filepath to anatomical file the mask is meant for

    Returns:
      mask_dat -- [Nibabel data] the mask data in Nibabel format
    """

    import nibabel as nib
    import numpy as np

    from workflow_utils import raise_smart_exception

    try:
        mask_img = nib.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask_dat = mask_img.get_data()
    ref_img = nib.load(ref_file)

    # Check that the specified mask is binary.
    mask_vals   = np.unique(mask_dat)
    if (mask_vals.size != 2) or not (mask_vals == [0, 1]).all():
        err = "Error: Mask is not binary, has %i unique val(s) of %s " \
              "(see file %s)" % (mask_vals.size, mask_vals, mask_file)
        raise_smart_exception(locals(),err)

    # Verify that the mask and anatomical images have the same dimensions.
    if ref_img.shape != mask_img.shape:
        err = "Error: Mask and anatomical image are different dimensions " \
              "for %s" % mask_file
        raise_smart_exception(locals(),err)

    # Verify that the mask and anatomical images are in the same space
    # (have the samme affine matrix)
    if (mask_img.get_affine() == ref_img.get_affine()).all == False:
        err = "Error: Mask and anatomical image are not in the same space " \
              "for %s vs %s" % (mask_file, ref_file)
        raise_smart_exception(locals(),err)

    return mask_dat


def create_anatomical_background_mask(anatomical_data, fg_mask_data, 
    exclude_zeroes=False):
    """Create a mask of the area outside the head in an anatomical scan by
    inverting a provided foreground mask.

    Keyword arguments:
      anatomical_data -- [Numpy array] an array of the raw anatomical data
      fg_mask_data -- [Numpy array] an array of binary foreground mask data
      exclude_zeroes -- [boolean] (default: False) flag to exclude pure zero 
                        values when creating the background mask

    Returns:
      bg_mask_data -- [Nibabel data] background mask data in Nibabel format
    """

    import numpy as np
    from workflow_utils import raise_smart_exception

    # invert the foreground mask
    try:
        bg_mask_data = 1 - fg_mask_data
    except Exception as e:
        err = "\n\n[!] Input data must be a NumPy array object, and not a " \
              "list.\n\nError details: %s\n\n" % e
        raise_smart_exception(locals(),err)

    if exclude_zeroes:
        # modify the mask to exclude zeroes in the background of the
        # anatomical image, as these are often introduced artificially and can
        # skew the QAP metric results
        bool_anat_data = anatomical_data > 0
        bg_mask_data = bg_mask_data * bool_anat_data

    return bg_mask_data
