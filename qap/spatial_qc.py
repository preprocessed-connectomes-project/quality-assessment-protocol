
import os
import sys
import numpy as np
import nibabel as nb
import pandas as pd
import scipy.ndimage as nd
import scipy.stats as stats


def summary_mask(anat_data, mask_data):
    """Will calculate the three values (mean, stdev, and size) and return them
    as a tuple.
    
    Keyword Arguments:
      anat_data -- [Numpy array] the anatomical scan data
      mask_data -- [Numpy array] the binary mask to mask the anatomical data 
                                 with
    
    Returns:
      (mean, std, size) -- [Python tuple] the summary values of the scan
    """
    
    import numpy as np
    
    anat_masked = anat_data[mask_data == 1].astype(np.float)
    mean        = anat_masked.mean()
    std         = anat_masked.std(ddof=1)
    size        = len(anat_masked)
    
    return (mean, std, size)


def check_datatype(background):
    """Process the image data to only include non-negative integer values.

    Keyword Arguments:
      background -- [Numpy array] the voxel values of the background (outside 
                    the head) of the anatomical image

    Returns:
      background -- [Numpy array] the input array with floats converted to 
                    integers and negative values set to zero
    """

    import numpy as np

    # If this is float then downgrade the data to an integer with some checks
    if np.issubdtype(background.dtype, float):
        background2 = background.astype('int32')    
        # Ensure downgrading datatype didn't really change the values
        if np.abs(background2.astype('float32') - background).mean() > 0.05:
            print "WARNING: Downgraded float to an int but values are " \
                  "different by more than 0.05"
        background = background2
        del background2
   
    # We only allow integer values for now
    if not np.issubdtype(background.dtype, int):
        print "QI1 can not be calculated for data that is not integer or " \
              "floating point: %s" % background.dtype
        raise TypeError
        
    # convert any negative voxel values to zero, provide warning
    for vox in background.flatten():
        if vox < 0:
            print "\nWARNING: Negative voxel values in anatomical scan " \
                  "converted to zero.\n"
            background = background.clip(0)
            break

    return background


def convert_negatives(img_data):
    """Convert any negative voxel values to zero and provide a warning.

    Keyword Arguments:
      img_data -- [Numpy array] the anatomical image's voxel values

    Return:
      img_data -- [Numpy array] the input array with negative values set to 
                  zero
    """

    for vox in img_data.flatten():
        if vox < 0:
            print "\nWARNING: Negative voxel values in anatomical scan " \
                  "converted to zero.\n"
            img_data = img_data.clip(0)
            break

    return img_data


def snr(mean_fg, std_bg):
    """Calculate the Signal-to-Noise Ratio (SNR) of an image.

    Keyword Arguments:
      mean_fg -- [float] the mean value of voxel intensities in the foreground
                 (either within the head or a particular tissue) of the image
      std_bg -- [float] the standard deviation of the voxel intensities of the
                background (outside the head) voxels
    
    Returns:
      snr -- [float] the signal-to-noise ratio

    Notes:
      - For anatomical images:
          SNR = (mean GM intensity) / (std of background intensities)
      - For functional images:
          SNR = (mean brain intensity) / (std of background intensities)    
    """

    snr     = mean_fg / std_bg

    return snr


def cnr(mean_gm, mean_wm, std_bg):
    """Calculate Contrast-to-Noise Ratio (CNR) of an image.
    
    Keyword Arguments:
      mean_gm -- [float] the mean value of the gray matter voxels
      mean_wm -- [float] the mean value of the white matter voxels
      std_bg -- [float] the standard deviation of the voxel intensities of the
                background (outside the head) voxels

    Returns:
      cnr -- [float] the contrast-to-noise ratio

    Notes:
      - CNR = |(mean GM intensity) - (mean WM intensity)| / 
                                            (std of background intensities)    
    """

    import numpy as np
    cnr     = np.abs(mean_gm - mean_wm)/std_bg

    return cnr


def cortical_contrast(mean_gm, mean_wm):
    """Calculate the vertex-wise cortical contrast.

    Keyword Arguments:
      mean_gm -- [float] the mean value of the gray matter voxels
      mean_wm -- [float] the mean value of the white matter voxels

    Returns:
      cort_con -- [float] the cortical contrast value
   
    cortical contrast = (mean WM intensity) - (mean GM intensity) /
                            ( (mean WM intensity + mean GM intensity) / 2 )
    """

    cort_con = (mean_wm - mean_gm) / ((mean_wm + mean_gm) / 2)

    return cort_con

    
def fber(anat_data, skull_mask_data, bg_mask_data):
    """Calculate the Foreground-to-Background Energy Ratio (FBER) of an image.

    Keyword Arguments:
      anat_data -- [Numpy array] the anatomical/spatial data of the image
      skull_mask_data -- [Numpy array] the binary mask defining the head
      bg_mask_data -- [Numpy array] the binary mask defining the background 
                      (outside of the head)

    Returns:
      fber -- [float] the foreground-to-background energy ratio

    Notes:
      - FBER = (mean foreground energy) / (mean background energy)
    """

    import numpy as np

    mean_fg = (np.abs(anat_data[skull_mask_data == 1]) ** 2).sum() / (skull_mask_data.sum())
    mean_bg = (np.abs(anat_data[bg_mask_data == 1]) ** 2).sum() / (bg_mask_data.size - bg_mask_data.sum())
    fber    = mean_fg / mean_bg

    return fber


def efc(anat_data):
    """Calculate the Entropy Focus Criterion of the image.

    Keyword Arguments:
      anat_data -- [Nibabel data] the anatomical image data

    Returns:
      efc -- [float] the entropy focus criterion value

    Notes:
      - EFC based on Atkinson 1997, IEEE TMI
      - We normalize the original equation by the maximum entropy so our EFC
        can be easily compared across images with different dimensions.
    """

    import numpy as np
        
    # let's get rid of those negative values
    anat_data = convert_negatives(anat_data)
        
    # Calculate the maximum value of the EFC (which occurs any time all 
    # voxels have the same value)
    efc_max = 1.0 * np.prod(anat_data.shape) * (1.0 / np.sqrt(np.prod(anat_data.shape))) * \
                np.log(1.0 / np.sqrt(np.prod(anat_data.shape)))
    
    # Calculate the total image energy
    b_max   = np.sqrt((anat_data**2).sum())
    
    # Calculate EFC (add 1e-16 to the image data to keep log happy)
    efc     = (1.0 / efc_max) * np.sum((anat_data / b_max) * np.log((anat_data + 1e-16) / b_max))
    
    if np.isnan(efc): 
        print "NaN found for efc (%3.2f,%3.2f)" % (efc_max,b_max)
    
    return efc


def artifacts(anat_data, fg_mask_data, bg_mask_data, calculate_qi2=False):
    """Calculates QI1, the fraction of total voxels that within artifacts.

    Keyword Arguments:
      anat_data -- [Nibabel data] the anatomical image data
      fg_mask_data -- [Nibabel data] the binary mask of the head
      bg_mask_data -- [Nibabel data] the binary mask of the background
      calculate_qi2 -- [boolean] (default: False) whether to calculate Qi2

    Returns:
      (QI1,QI2) -- [Python tuple] the Qi1 and Qi2 values (Qi2 = None if not 
                   calculated)

    Notes:
      - Detect artifacts in the anatomical image using the method described in
        Mortamet et al. 2009 (MRM).
      - Optionally, also calculates QI2, the distance between the distribution 
        of noise voxel (non-artifact background voxels) intensities, and a 
        Ricean distribution.
    """

    import numpy as np

    # Create an image containing only background voxels (everything 
    # outside bg_mask set to 0)
    background = anat_data.copy()
    background[bg_mask_data != 1] = 0
    
    # make sure the datatype is an int
    background = check_datatype(background)
       
    # Find the background threshold (the most frequently occurring value 
    # excluding 0)
    bg_counts       = np.bincount(background.flatten())
    bg_threshold    = np.argmax(bg_counts[1:]) + 1

    # Apply this threshold to the background voxels to identify voxels
    # contributing artifacts. 
    background[background <= bg_threshold] = 0
    background[background != 0] = 1

    # Create a structural element to be used in an opening operation.
    struct_elmnt    = np.zeros((3,3,3))
    struct_elmnt[0,1,1] = 1
    struct_elmnt[1,1,:] = 1
    struct_elmnt[1,:,1] = 1
    struct_elmnt[2,1,1] = 1

    # Perform an opening operation on the background data.
    background      = nd.binary_opening(background, structure=struct_elmnt)

    # Count the number of voxels that remain after the opening operation. 
    # These are artifacts.
    QI1             = background.sum() / float(bg_mask_data.sum())
    
    ''' "bg" in code below not defined- need to ascertain what that should '''
    '''      be, and correct it- unit test for this part disabled for now  '''
    if calculate_qi2:
        # Now lets focus on the noise, which is everything in the background
        # that was not identified as artifact
        bgNoise     = anat_data[(fg_mask_data-bg)==1]

        # calculate the histogram of the noise and its derivative
        H           = np.bincount(bgNoise)
        H           = 1.0*H/H.sum()
        dH          = H[1:]-H[:-1]

        # find the first value on the right tail, i.e. tail with negative
        # slope, i.e. dH < 0 that is less than or equal to half of the
        # histograms max
        firstNegSlope = np.nonzero(dH<0)[0][0]
        halfMaxRightTail = np.nonzero(H[firstNegSlope:]<(H.max()/2))[0][0]

        # divide by the standard deviation
        bgNoiseZ    = bgNoise / bgNoise.std()
        bgChiParams = ss.chi.fit(bgNoiseZ)
        #print bgChiParams
    
        # now generate values that are consistent with the histogram
        yx          = range(0,H.size)/bgNoise.std()
        rvs         = ss.chi.pdf(yx,bgChiParams[0],loc=bgChiParams[1],scale=bgChiParams[2])

        # now we can calculate the goodness of fit
        gof         = np.average(np.absolute(H[halfMaxRightTail:]-rvs[halfMaxRightTail:]))
        QI2         = QI1+gof
    else:
        QI2         = None

    return (QI1,QI2)


def fwhm(anat_file, mask_file, out_vox=False):
    """Calculate the FWHM of the input image using AFNI's 3dFWHMx.
    
    Keyword Arguments:
      anat_file -- [string] the filepath to the anatomical image NIFTI file
      mask_file -- [string] the filepath to the binary head mask NIFTI file
      out_vox -- [boolean] (default: False) output the FWHM as # of voxels 
                 instead of mm (the default)

    Returns:
      tuple(vals) -- [Python tuple] a tuple of the FWHM values (x, y, z, and
                     combined)

    Notes:
      - Uses AFNI 3dFWHMx. More details here:
          https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dFWHMx.html
    """

    import commands
    import nibabel as nib
    import numpy as np
    from scipy.special import cbrt
    
    # call AFNI command to get the FWHM in x,y,z and combined
    cmd     = "3dFWHMx -combined -mask %s -input %s" % (mask_file, anat_file)
    out     = commands.getoutput(cmd)
    
    # extract output
    line    = out.splitlines()[-1].strip()
    vals    = np.array(line.split(), dtype=np.float)
    
    if out_vox:
        # get pixel dimensions
        img     = nib.load(anat_file)
        hdr     = img.get_header()
        pixdim  = hdr['pixdim'][1:4]
    
        # convert to voxels
        pixdim  = np.append(pixdim, cbrt(pixdim.prod()))
        # get the geometrix mean
        vals    = vals / pixdim
    
    return tuple(vals)


def ghost_direction(epi_data, mask_data, direction="y", ref_file=None,
                    out_file=None):
    """Calculate the Ghost to Signal Ratio of EPI images.

    Keyword Arguments:
      epi_data -- [Nibabel data] the mean of the functional timeseries
      mask_data -- [Nibabel data] the functional brain binary mask data
      direction -- [string] (default: 'y') the phase encoding direction of the
                   EPI image
      ref_file -- [string] (default: None) if you are saving the Nyquist ghost
                  mask, this is the filepath of the reference file to use the
                  header for the ghost mask NIFTI file
      out_file -- [string] (default: None) if you are saving the Nyquist ghost
                  mask, this is the filepath to the ghost mask NIFTI file to 
                  be written

    Returns:
      gsr -- [float] the ghost-to-signal ratio value

    Notes:
      - GSR from Giannelli 2010. More details here:
          https://www.ncbi.nlm.nih.gov/pubmed/21081879
      - This should be used for EPI images where the phase encoding direction
        is known.
    """
    
    import numpy as np
    
    # first we need to make a nyquist ghost mask, we do this by circle 
    # shifting the original mask by N/2 and then removing the intersection
    # with the original mask
    n2_mask_data    = np.zeros_like(mask_data)
    
    ## rotate by n/2
    if direction == "x":
        n2                      = np.floor(mask_data.shape[0]/2)
        n2_mask_data[:n2,:,:]   = mask_data[n2:(n2*2),:,:]
        n2_mask_data[n2:(n2*2),:,:]   = mask_data[:n2,:,:]
    elif direction == "y":
        n2                      = np.floor(mask_data.shape[1]/2)
        n2_mask_data[:,:n2,:]   = mask_data[:,n2:(n2*2),:]
        n2_mask_data[:,n2:(n2*2),:]   = mask_data[:,:n2,:]
    elif direction == "z":
        n2                      = np.floor(mask_data.shape[2]/2)
        n2_mask_data[:,:,:n2]   = mask_data[:,:,n2:(n2*2)]
        n2_mask_data[:,:,n2:(n2*2)]   = mask_data[:,:,:n2]        
    else:
        raise Exception("Unknown direction %s, should be x, y, or z" \
                        % direction)
    
    ## now remove the intersection with the original mask
    n2_mask_data    = n2_mask_data * (1-mask_data)
    
    ## now create a non-ghost background region, that contains 2s
    n2_mask_data    = n2_mask_data + 2*(1-n2_mask_data-mask_data)
    
    # Save mask
    if ref_file is not None and out_file is not None:
        import nibabel as nib
        ref = nib.load(ref_file)
        out = nib.Nifti1Image(n2_mask_data, ref.get_affine(), ref.get_header()) 
        out.to_filename(out_file)
   
    # now we calculate the Ghost to signal ratio, but here we define signal
    # as the entire foreground image
    gsr             = (epi_data[n2_mask_data==1].mean() - epi_data[n2_mask_data==2].mean())/epi_data[n2_mask_data==0].mean()

    
    return gsr


def ghost_all(epi_data, mask_data):
    """Call the 'ghost_direction' function on all possible phase encoding 
    directions.

    Keyword Arguments:
      epi_data -- [Nibabel data] the mean of the functional timeseries
      mask_data -- [Nibabel data] the functional brain binary mask data

    Returns:
      tuple(gsrs + [None]) -- [Python tuple] the ghost-to-signal ratios of 
                              each phase encoding direction
    """
    
    directions = ["x", "y"]
    gsrs = [ ghost_direction(epi_data, mask_data, d) for d in directions ]
    
    return tuple(gsrs + [None])


