# Quality Assessment Protocol

Here we detail running our different QC measures. We assume that the [quality-assessment-protocol repository](https://github.com/preprocessed-connectomes-project/quality-assessment-protocol) is in your python path as `qap`.

## Load Data

You can use some wrapper functions to load and mask the data. QC measures use five different files as inputs: anatomical data, anatomical head mask, functional time-series data, mean functional data, and functional head mask. Here, we show how to load or calculate each of those files. Ideally, the anatomical and functional data should be as 'raw' as possible. With our CPAC pipeline, we use anatomical data that has only been reoriented and use functional data that has been slice time and motion corrected.

	from qap import load_image, load_mask, load_func

	anat_file			= "test_anat.nii.gz"
	anat_mask_file	= "test_anat_mask.nii.gz"
	anat_data			= load_image(anat_file)
	anat_mask			= load_mask(anat_mask_file, anat_file)
	
	func_file			= "test_func.nii.gz"
	func_mask_file	= "test_anat_mask.nii.gz"
	func_data			= load_func(func_file, func_mask_file)
	

The mean functional data can be generated from the functional time-series via the following code, which will generate the mean functional image without any masking (e.g., preserving values outside the head, which are needed for various spatial metrics).

	from qap import calc_mean_func, load_mask
	
	func_file 		= "test_func.nii.gz"
	func_mask_file	= "test_anat_mask.nii.gz"
	mean_func 		= calc_mean_func(func_file)
	func_mask			= load_mask(func_mask_file)
	
Note that we will be using the above anatomical and functional outputs for our metrics described below.


## Spatial Anatomical

We calculate the measures below with the associated variable in square brackets followed by a brief description and possible link to a reference.

* **Contrast to Noise Ratio (CNR) [anat_cnr]:** Calculated as the mean of the gray matter values minus the mean of the white matter values, divided by the standard deviation of the air values [^1].
* **Entopy Focus Criterion (EFC) [anat_efc]:** Shannon’s entropy is used to summarize the principal directions distribution, higher energy indicating the distribution is more uniform (i.e., less noisy) [^2].
* **Foreground to Background Energy Ratio (FBER) [anat_fber]:** Mean energy of image values (i.e., mean of squares) within the head relative to outside the head.
* **Smoothness of Voxels (FWHM) [anat_fwhm]:** The full-width half maximum (FWHM) of the spatial distribution of the image intensity values [^3].
* **Artifact Detection (Qi1) [anat_qi1]:** The proportion of voxels with intensity corrupted by artifacts normalized by the number of voxels in the background [^4].
* **Signal-to-Noise Ratio (SNR) [anat_snr]:** The mean of image values within gray matter divided by the standard deviation of the image values within air (i.e., outside the head) [^1].

### CNR

We first compute the mean image values within gray matter (mean_gm) and within white-matter (mean_wm), and we will divide the difference of those values by standard deviation within the air (std_bg). Note that the air is considered anything outside the head mask.

	from qap import load_mask, summary_mask, cnr
	
	# Load the different masks
	fg_mask_data = anat_mask_data
	bg_mask_data = 1 - fg_mask_data
	gm_mask_data = load_mask(gm_mask_file, anat_mask_file)
	wm_mask_data = load_mask(wm_mask_file, anat_mask_file)
	
	# Calculate mean grey-matter, mean white-matter 
	# and standard deviation air
	mean_gm,_,_	= summary_mask(anat_data, gm_mask_data)
	mean_wm,_,_	= summary_mask(anat_data, wm_mask_data)
	_,std_bg,_	= summary_mask(anat_data, bg_mask_data)
	
	# SNR
	anat_cnr 		= cnr(mean_gm, mean_wm, std_bg)

### EFC

	from qap import etc
	anat_efc		= efc(anat_data)

### FBER

	from qap import fber
	anat_fber 	= fber(anat_data, anat_mask_data):

### FWHM

While the other spatial metrics exclusively use python code, the smoothness measures makes use of AFNI's command-line tool `3dFWHMx`. The option `out_vox` to the function `fwhm` indicates if the fwhm output is given in mm (`out_vox=False`) or in number of voxels (`out_vox=True`).

	from qap import fwhm
	anat_fwhm 	= fwhm(anat_file, anat_mask_file, out_vox=False)

### Qi1

	from qap import artifacts
	anat_qi1		= artifacts(anat_data, anat_mask_data, calculate_qi2=False)


### SNR

We first compute the mean image values within gray matter (mean_gm) and the standard deviation within the air (std_bg). Note that the air is considered anything outside the head mask.

	from qap import load_mask, summary_mask, snr
	
	# Load the different masks
	fg_mask_data = anat_mask_data
	bg_mask_data = 1 - fg_mask_data
	gm_mask_data = load_mask(gm_mask_file, anat_mask_file)
	
	# Calculate mean grey-matter and standard deviation air
	mean_gm,_,_	= summary_mask(anat_data, gm_mask_data)
	_,std_bg,_	= summary_mask(anat_data, bg_mask_data)
	
	# SNR
	anat_snr 		= snr(mean_gm, std_bg)


## Spatial Functional

The spatial functional measures will make use of the mean functional image and include EFC, FBER, and FWHM (code which has been described above) along with GSR (code which has been described below).

* **Entopy Focus Criterion [func_efc]:** Shannon’s entropy is used to summarize the principal directions distribution, higher energy indicating the distribution is more uniform (i.e., less noisy) [^2]. _Uses mean functional._
* **Foreground to Background Energy Ratio [func_fber]:** Mean energy of image values (i.e., mean of squares) within the head relative to outside the head. _Uses mean functional._
* **Smoothness of Voxels [func_fwhm]:** The full-width half maximum (FWHM) of the spatial distribution of the image intensity values. _Uses mean functional._
* **Ghost to Signal Ratio (GSR) [func_gsr]:** A measure of the mean signal in the ‘ghost’ image (signal present outside the brain due to acquisition in the phase encoding direction) relative to mean signal within the brain. _Uses mean functional._

### GSR

You should know the phase encoding direction to decide if you want to use `func_ghost_x` (RL/LR) or `func_ghost_y` (AP/PA).
	
	from qap import ghost_all
	func_ghost_x, func_ghost_y, _ = ghost_all(mean_func_data, func_mask_data)


## Temporal Functional

* **Standardized DVARS [func_dvars]:** The spatial standard deviation of the temporal derivative of the data, normalized by the temporal standard deviation and temporal autocorrelation [^5][^6]. _Uses functional time-series._
* **Outlier Detection [func_outlier]:** The mean fraction of outliers found in each volume using 3dTout command in AFNI (http://afni.nimh.nih.gov/afni) [^7]. _Uses functional time-series._
* **Median Distance Index [func_quality]:** The mean distance (1 – spearman’s rho) between each time-point's volume and the median volume using AFNI’s 3dTqual command (http://afni.nimh.nih.gov/afni) [^7]. _Uses functional time-series._
* **Mean Fractional Displacement - Jenkinson [func_mean_fd]:** A measure of subject head motion, which compares the motion between the current and previous volumes. This is calculated by summing the absolute value of displacement changes in the x, y and z directions and rotational changes about those three axes. The rotational changes are given distance values based on the changes across the surface of a 50mm radius sphere [^5][^8]. _Uses functional time-series._
* **Number of volumes with FD greater than 0.2mm [func_num_fd]:** _Uses functional time-series._
* **Percent of volumes with FD greater than 0.2mm [func_perc_fd]:** _Uses functional time-series._
* **Median temporal Signal to Noise ratio [tsnr]:** _Uses functional time-series and a mask._


### Standardized DVARS

The `output_all` option if `True` will spit out a three column matrix with standardized DVARS, regular DVARS, and a voxelwise standardization of DVARS.

	from qap import calc_dvars
	func_dvars	= calc_dvars(func_data, output_all=False)

### Outlier Detection

The `out_fraction` option if `True` will return the mean _fraction_ of time-points per voxel that are outliers, whereas `False` will return the mean _number_ of time-points per voxel that are outliers.

	from qap import mean_outlier_timepoints
	func_outlier	= mean_outlier_timepoints(func_file, func_mask_file, out_fraction=True)

### Median Distance Index

The `mask` allows you to specify a mask file. If set to "auto" mask will be calculated from data. If not set the whole dataset will be used.

	from qap import mean_quality_timepoints
	func_quality	= mean_quality_timepoints(func_file, mask="auto")

### FD

Here we describe computing `mean_fd`, `num_fd`, and `perc_fd`. This requires that you have the input coordinate transforms that is output by AFNI's `3dvolreg` during motion correction. The option `threshold` sets the threshold for determining the number and percent of volumes with FD greater than said threshold.

	from qap import summarize_fd
	mean_fd, num_fd, perc_fd = summarize_fd(motion_matrix_file, threshold=0.2)
	
### tSNR

tSNR is a ratio of signal (mean across time) and noise (standard deviation around the mean across time) averaged over all voxels withing the specified mask

    from qap import median_tsnr
	tsnr = median_tsnr(func_data, func_mask_data)


## Determining Outliers

If you run the above procedure on an array of subjects, then you can take 1.5x or 3x the inter-quartile range to determine subjects that are outliers.


## References

[^1]: Magnotta, V. A., & Friedman, L. (2006). Measurement of signal-to-noise and contrast-to-noise in the fBIRN multicenter imaging study. Journal of Digital Imaging, 19(2), 140-147.

[^2]: Farzinfar, M., Dietrich, C., Smith, R. G., Li, Y., Gupta, A., Liu, Z., & Styner, M. A. (2012, May). Entropy based DTI quality control via regional orientation distribution. In Biomedical Imaging (ISBI), 2012 9th IEEE International Symposium on (pp. 22-25). IEEE.

[^3]: Friedman, L., Stern, H., Brown, G. G., Mathalon, D. H., Turner, J., Glover, G. H., ... & Potkin, S. G. (2008). Test–retest and between‐site reliability in a multicenter fMRI study. Human brain mapping, 29(8), 958-972.

[^4]: Mortamet, B., Bernstein, M. A., Jack, C. R., Gunter, J. L., Ward, C., Britson, P. J., ... & Krueger, G. (2009). Automatic quality assessment in structural brain magnetic resonance imaging. Magnetic Resonance in Medicine, 62(2), 365-372.

[^5]: Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L. & Petersen, S. E. Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage 59, 2142-2154 (2012).

[^6]: Nichols, T. (2012, Oct 28). Standardizing DVARS. Retrieved from http://blogs.warwick.ac.uk/nichols/entry/standardizing_dvars.

[^7]: Cox, R.W.. AFNI: Software for analysis and visualization of functional magnetic resonance neuroimages. Computers and Biomedical Research, 29:162-173, 1996.

[^8]: Jenkinson, M., Bannister, P., Brady, M., & Smith, S. (2002). Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage, 17(2), 825-841.

[^9]: Giannelli, M., Diciotti, S., Tessa, C., & Mascalchi, M. (2010). Characterization of Nyquist ghost in EPI-fMRI acquisition sequences implemented on two clinical 1.5 T MR scanner systems: effect of readout bandwidth and echo spacing. Journal of Applied Clinical Medical Physics, 11(4).
