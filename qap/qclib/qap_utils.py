def load_image(image_file):

    import nibabel as nib
    import numpy as np

    '''
    inputs
        image_file: path to an image, usually a structural or functional scan

    outputs
        dat: the image data in Nibabel format
    '''
   
    img = nib.load(image_file)
    dat = img.get_data()

    # Ensure that data is cast as at least 32-bit
    if np.issubdtype(dat.dtype, float):

        dat = dat.astype('float32')
        
        # Check for negative values
        if (dat < 0).any():
            print "found negative values, setting to zero (see file: %s)" % image_file
            dat[dat<0] = 0
    
    elif np.issubdtype(dat.dtype, int):

        dat = dat.astype('int32')
    
    else:

        msg = "Error: Unknown datatype %s" % dat.dtype
        raise Exception(msg)
    
    return dat



def load_mask(mask_file, ref_file):

    import nibabel as nib
    import numpy as np

    '''
    inputs
        mask_file: binarized mask file
        ref_file: anatomical file the mask is meant for

    outputs
        mask_dat: the mask file in Nibabel format
    '''

    mask_img = nib.load(mask_file)
    mask_dat = mask_img.get_data()
    ref_img = nib.load(ref_file)
    
    # Check that the specified mask is binary.
    mask_vals   = np.unique(mask_dat)
    if (mask_vals.size != 2) or not (mask_vals == [0, 1]).all():
        err = "Error: Mask is not binary, has %i unique val(s) of %s " \
              "(see file %s)" % (mask_vals.size, mask_vals, mask_file)
        raise Exception(err)
    
    # Verify that the mask and anatomical images have the same dimensions.
    if ref_img.shape != mask_img.shape:
        err = "Error: Mask and anatomical image are different dimensions " \
              "for %s" % mask_file
        raise Exception(err)

    # Verify that the mask and anatomical images are in the same space (have the samme affine matrix)
    if (mask_img.get_affine() == ref_img.get_affine()).all == False:
        err = "Error: Mask and anatomical image are not in the same space " \
              "for %s vs %s" % (mask_file, ref_file)
        raise Exception(err)
    
    return mask_dat
