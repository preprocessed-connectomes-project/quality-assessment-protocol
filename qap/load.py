import os
import numpy as np
import nibabel as nib

###
# Loading Data Functions
###

def load_image(image_file):
    """
    Loads an analyze/nifti image, generally for 3D images. 
    Casts input as 32-bit float (if float) or 32-bit uint (if int).
    
    Paramaters
    ----------
    image_file: str
        path to data
    
    Returns
    -------
    dat: nparray
        Image as numpy array (i.e., 3D array)
    """
    img         = nib.load(image_file)
    dat         = img.get_data()
    # Ensure that data is cast as at least 32-bit
    if np.issubdtype(dat.dtype, float):
        dat = dat.astype('float32')
        # Check for negative values
        if (dat < 0).any():
            print "found negative values, setting to zero (see file: %s)" % image_file
            dat[dat<0] = 0
    elif np.issubdtype(dat.dtype, int):
        dat = dat.astype('uint32')
    else:
        msg = "Error: Unknown datatype %s" % dat.dtype
        print msg
        raise Exception(msg)
    return(dat)

def load_mask(mask_file, ref_file=None):
    """
    Loads the mask for the data. Also performs some checks.
    
    Paramaters
    ----------
    mask_file: str
        path to mask data
    ref_file: str
        path to reference data, generally data that will be masked (optional)
    
    Returns
    -------
    mask_dat: nparray
        Image as numpy array (i.e., 3D array)
    """
    mask_img    = nib.load(mask_file)
    mask_dat    = mask_img.get_data()
    
    # Check that the specified mask is binary.
    mask_vals   = np.unique(mask_dat)
    if (mask_vals.size != 2) or not (mask_vals == [0, 1]).all():
        print("Error: Mask is not binary")
        raise Exception("")
    
    if ref_file is not None:
        ref_img     = nib.load(ref_file)
        # Verify that the mask and anatomical images have the same dimensions.
        if ref_img.shape != mask_img.shape:
            raise Exception("Error: Mask and anatomical image are different dimensions")
        # Verify that the mask and anatomical images are in the same space (have the samme affine matrix)
        if (mask_img.get_affine() == ref_img.get_affine()).all == False:
            raise Exception("Error: Mask and anatomical image are not in the same space")
    
    return mask_dat

def load_func(func_file, mask_file):
    """
    Loads and masks functional time-series data. Data is returned as a 2D 
    matrix with number of time-points (rows) x number of voxels (columns).
    
    Paramaters
    ----------
    func_file: str
        path to functional data
    mask_file: str
        path to mask data
    
    Returns
    -------
    func: nparray
        Functional time-series as 2D numpy array with ntpts x nvoxs
    """
    func_dat    = load_image(func_file)
    mask        = load_mask(mask_file)
    func        = func_dat.astype(np.float)
    if len(func.shape) != 4:
        raise Exception("Input functional %s should be 4-dimensional" % func_file)
    func        = func[mask.nonzero()].T # will have ntpts x nvoxs
    return(func)

def calc_mean_func(func_file):
    """
    Computes the mean functional data for use in spatial metrics. Note that 
    data is not masked here.
    
    Paramaters
    ----------
    func_file: str
        path to functional data
    
    Returns
    -------
    mean_func: nparray
        Mean functional 3D data
    """
    func_dat    = load_image(func_file)
    func        = func_dat.astype(np.float)
    if len(func.shape) != 4:
        raise Exception("Input functional %s should be 4-dimensional" % func_file)
    mean_func   = func.mean(3)  # double check
    return(mean_func)
