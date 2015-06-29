---
layout: page
title: PCP Quality Assessment Protocol
---

{ introduction here with more background info }

**Table of Contents**

* [System Requirements](#system-requirements)
* [Installing the QAP Package](#installing-the-qap-package)
* [A Description of QA Measures](#a-description-of-qa-measures)
* [Pipeline Configuration YAML Files](#pipeline-configuration-yaml-files)
* [Subject List YAML Files](#subject-list-yaml-files)
* [Running the QAP Pipelines](#running-the-qap-pipelines)
* [Running the QAP Pipelines on AWS Cloud Instances](#running-the-qap-pipelines-on-aws-amazon-cloud-instances)
* [References](#references)

## System Requirements

* Any *nix-based operating system capable of running QAP's dependencies will work.  

## Installing the QAP Package

Before you can install QAP, you must first install a number of key dependencies.

### Application Dependencies

QAP requires AFNI, C3D, and FSL to run. Links to installation instructions for AFNI and FSL are listed below:

* [AFNI Installation](#http://afni.nimh.nih.gov/pub/dist/HOWTO/howto/ht00_inst/html)
* [FSL Installation](#http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FslInstallation)

If you are using a Debian-based Linux distribution, you can use `apt-get` by first adding Neurodebian to the apt repository list and then installing the Neurodebian FSL and AFNI packages:

    wget -O- http://neuro.debian.net/lists/$(lsb_release -cs).us-nh.full | tee /etc/apt/sources.list.d/neurodebian.sources.list
    apt-key adv --recv-keys --keyserver pgp.mit.edu 2649A5A9
    apt-get update
    apt-get install -y fsl-5.0-complete afni

To install C3D, download the appropriate version for your platform from [Sourceforge](#http://sourceforge.net/projects/c3d/).  Unzip or mount the archive containing C3D and copy it to a memorable location.  Finish by adding the following line to your `.bashrc` file:

    export PATH=/path_to/C3D/bin:$PATH

### Python Dependencies

QAP requires Numpy, Scipy, Nipype, Nibabel and Nitime.  If you have `pip`, you may install these by typing in the command below:

    pip install numpy scipy nipype nibabel nitime

### QAP

Once the pre-requisites have been satisfied, you can install the QAP package itself:

    cd /tmp
    git clone https://github.com/preprocessed-connectomes-project/quality-assessment-protocol.git qap
    cd /qap
    python setup.py install

## A Description of QA Measures

There are three collections of measures that can be run using the QAP software package:

* Anatomical spatial measures
* Functional (4D timeseries) spatial measures
* Functional (4D timeseries) temporal measures

### Spatial Anatomical

* **Contrast to Noise Ratio (CNR) [anat_cnr]:** Calculated as the mean of the gray matter values minus the mean of the white matter values, divided by the standard deviation of the air values, higher values are better [^1].
* **Entropy Focus Criterion (EFC) [anat_efc]:** Uses the Shannon entropy of voxel intensities as an indication of ghosting and blurring induced by head motion, lower is better [^2].
* **Foreground to Background Energy Ratio (FBER) [anat_fber]:** Mean energy of image values (i.e., mean of squares) within the head relative to outside the head, higher values are better.
* **Smoothness of Voxels (FWHM) [anat_fwhm]:** The full-width half maximum (FWHM) of the spatial distribution of the image intensity values in units of voxels, lower values are better [^3].
* **Artifact Detection (Qi1) [anat_qi1]:** The proportion of voxels with intensity corrupted by artifacts normalized by the number of voxels in the background, lower values are better [^4].
* **Signal-to-Noise Ratio (SNR) [anat_snr]:** The mean of image values within gray matter divided by the standard deviation of the image values within air (i.e., outside the head), higher values are better [^1].

### Spatial Functional

* **Entropy Focus Criterion [func_efc]:** Uses the Shannon entropy of voxel intensities as an indication of ghosting and blurring induced by head motion, lower is better [^2]. _Uses mean functional._
* **Foreground to Background Energy Ratio [func_fber]:** Mean energy of image values (i.e., mean of squares) within the head relative to outside the head, higher values are better. _Uses mean functional._
* **Smoothness of Voxels [func_fwhm]:** The full-width half maximum (FWHM) of the spatial distribution of the image intensity values in units of voxels, lower values are better. _Uses mean functional._
* **Ghost to Signal Ratio (GSR) [func_gsr]:** A measure of the mean signal in the ‘ghost’ image (signal present outside the brain due to acquisition in the phase encoding direction) relative to mean signal within the brain, lower values are better. _Uses mean functional._

### Temporal Functional

* **Standardized DVARS [func_dvars]:** The spatial standard deviation of the temporal derivative of the data, normalized by the temporal standard deviation and temporal autocorrelation, lower values are better [^5][^6]. _Uses functional time-series._
* **Outlier Detection [func_outlier]:** The mean fraction of outliers found in each volume using the [3dTout](http://afni.nimh.nih.gov/pub/dist/doc/program_help/3dToutcount.html) command from [AFNI](http://afni.nimh.nih.gov/afni), lower values are better [^7]. _Uses functional time-series._
* **Median Distance Index [func_quality]:** The mean distance (1 – spearman’s rho) between each time-point's volume and the median volume using AFNI’s 3dTqual command (http://afni.nimh.nih.gov/afni), lower values are better [^7]. _Uses functional time-series._
* **Mean Fractional Displacement - Jenkinson [func_mean_fd]:** A measure of subject head motion, which compares the motion between the current and previous volumes. This is calculated by summing the absolute value of displacement changes in the x, y and z directions and rotational changes about those three axes. The rotational changes are given distance values based on the changes across the surface of a 80mm radius sphere, lower values are better [^8][^10]. _Uses functional time-series._
* **Number of volumes with FD greater than 0.2mm [func_num_fd]:** Lower values are better _Uses functional time-series._
* **Percent of volumes with FD greater than 0.2mm [func_perc_fd]:** Lower values are better _Uses functional time-series._

## Pipeline Configuration YAML Files

Certain pre-processed files derived from the raw data are required to calculate these measures. By default, the QAP software package will generate these pre-requisite files given the raw data you have (anatomical/structural scans for the anatomical measures, 4D anatomical+timeseries scans for the functional). A preprocessing pipeline will be constructed to create these files, and this pipeline can be customized with the pipeline configuration file you provide.

Some examples of customizable features include segmentation thresholds for anatomical preprocessing, and the option to run slice timing correction for functional preprocessing. Computer resource allocation can also be customized using the configuration files, such as dedicating multiple cores/threads to processing, etc.

Templates for these files are provided in the /configs folder in the QAP main repository directory. Below is a list of options which can be configured for each of the pipelines.

### General (both types)

* **num_cores_per_subject**: Number of cores (on a single machine) or slots on a node (cluster/grid) per subject (or per instance of the pipeline). Slots are cores on a cluster/grid node. Dedicating multiple nodes allows each subject's processing pipeline to run certain operations in parallel to save time. 
* **num_subjects_at_once**: Similar to *num_cores_per_subject*, except this determines how many pipelines to run at once.   
* **output_directory**: The directory to write output files to.
* **working_directory**: The directory to store intermediary processing files in.

### Anatomical pipelines

* **num_ants_threads**: Number of cores to dedicate to ANTS anatomical registration. More cores will result in a faster registration process.
* **template_brain_for_anat**: Template brain to be used during anatomical registration, as a reference.
* **csf_threshold**: A decimal value - only voxels with a CSF probability greater than this value will be classified as CSF.            
* **gm_threshold**: A decimal value - only voxels with a gray matter probability greater than this value will be classified as gray matter.
* **wm_threshold**: A decimal value - only voxels with a white matter probability greater than this value will be classified as white matter.

### Functional pipelines

* **start_idx**: This allows you to select an arbitrary range of volumes to include from your 4-D functional timeseries. Enter the number of the first timepoint you wish to include in the analysis. Enter *0* to include the first volume.          
* **stop_idx**: This allows you to select an arbitrary range of volumes to include from your 4-D functional timeseries. Enter the number of the last timepoint you wish to include in the analysis. Enter *End* to include the final volume. Enter *0* in start_idx and *End* in stop_idx to include the entire timeseries. 
* **slice_timing_correction**: Whether or not to run slice timing correction - *True* or *False*. Interpolates voxel timeseries so that sampling occurs at the same time.

Multiply *num_cores_per_subject*, *num_subjects_at_once*, and *num_ants_threads* for the maximum amount of cores that could potentially be used during the pipeline run. (Or just *num_cores_per_subject* and *num_subjects_at_once* for functional pipeline runs).

## Subject List YAML Files

### Providing Raw Data

The QAP pipelines take in subject list YAML (.yml) files as an input. The filepaths to your raw data are defined in these subject lists, and these .yml files can be easily generated using the *qap_raw_data_sublist_generator.py* script included in the QAP software package. After installing the QAP software package, this script can be run from any directory. This subject list generator script assumes a specific file structure for your input data:

	/data_folder/site_name/subject_id/session_id/scan_id/file.nii.gz

The script will go through your data folder structure and generate the subject list YAML file for you. Note that these subject lists can also be created or edited by hand if you wish, though this can be cumbersome for larger data sets. For reference, an example of the subject list format follows:

	'1019436':
	  session_1:
	    anatomical_scan:
	      anat_1: /test-data/site_1/1019436/session_1/anat_1/mprage.nii.gz
	'2014113':
	  session_1:
	    anatomical_scan:
	      anat_1: /test-data/site_1/2014113/session_1/anat_1/mprage.nii.gz
	'3154996':
	  session_1:
	    anatomical_scan:
	      anat_1: /test-data/site_1/3154996/session_1/anat_1/mprage.nii.gz

Note that *anatomical_scan* is the label for the type of resource (in this case, anatomical raw data for the anatomical spatial QAP measures), and *anat_1* is the name of the scan. There can be multiple scans, and the QAP pipeline will output a CSV spreadsheet with calculated values in rows, with each subject-session-scan combination having its own row of values.

### Providing Already Pre-Processed Data

Alternatively, if you have already preprocessed some or all of your raw data, you can provide these already-existing files as inputs directly to the QAP pipelines via your subject list manually. If these files were processed using the CPAC software package, there is a script named *qap_cpac_output_sublist_generator.py* which will create a subject list YAML file pointing to these already generated files. The QAP pipelines will then use these files and skip any pre-processing steps involved in creating them, saving time and allowing you to use your own method of processing your data.

Below is a list of intermediary files used in the steps leading to the final QAP measures calculations. If you already have some of these processed for your data, they can be included in the subject list with the label on the left. For example, if you've already deobliqued, reoriented and skull-stripped your anatomical scans, you would list them in your subject list YAML file like so:

	anatomical_brain:  /path/to/image.nii.gz

###Anatomical Spatial measures workflow resources

* **anatomical_reorient**: anatomical (structural) scan that has been deobliqued and reoriented to RPI (.nii/.nii.gz)
* **anatomical_brain**: deobliqued & reoriented anatomical which has been skull-stripped (.nii/.nii.gz)
* **ants_initial_xfm**: first of three warp matrix files output by ANTS linear registration (.mat)
* **ants_rigid_xfm**: second of three warp matrix files output by ANTS (.mat)
* **ants_affine_xfm**: third of three warp matrix files output by ANTS (.mat)
* **ants_linear_warped_image**:	the ANTS-warped anatomical scan (.nii/.nii.gz)
* **anatomical_csf_mask**: segmentation mask of the anatomical scan's CSF (.nii/.nii.gz)
* **anatomical_gm_mask**: segmentation mask of the anatomical scan's gray matter (.nii/.nii.gz)
* **anatomical_wm_mask**: segmentation mask of the anatomical scan's white matter (.nii/.nii.gz)
* **qap_head_mask**: a whole-skull binarized mask

in the QAP head mask workflow, we also mask the background immediately in front of the scan participant's mouth - we do this to exclude breathing-induced noise from the calculation for the FBER QAP measure

### Functional Spatial measures workflow resources

* **func_motion_correct**: motion-corrected 4-D functional timeseries (.nii/.nii.gz)
* **functional_brain_mask**: a binarized mask of the functional scan (.nii/.nii.gz)
* **mean_functional**: a 3-D file containing the mean of the functional 4-D timeseries (.nii/.nii.gz)

### Functional Temporal measures workflow resources

* **func_motion_correct**: motion-corrected 4-D functional timeseries (.nii/.nii.gz)
* **functional_brain_mask**: a binarized mask of the functional scan (.nii/.nii.gz)
* **coordinate_transformation**: the matrix transformation from base to input coordinates, produced during motion correction (.aff12.1D)

Note that these are complete lists- obviously, not all intermediary files are required if you choose to provide them. For example, if you provide the skull-stripped *anatomical_brain*, then *anatomical_reorient* would not be necessary, and the pipeline will skip all steps before skull-stripping. Alternatively, if you do not have the *functional_brain_mask* for either of the functional pipelines, providing the *func_motion_correct* file will allow the pipeline to create it for you. Having none of these will simply cause the pipeline to take in the original *functional_scan* and produce all of these files on its own.

## Running the QAP Pipelines

There is a launch script for each of these measures, with each one featuring a similar interface. The Python-friendly YAML file format is used for the input subject list and pipeline configuration files. You can use these scripts from the command line, from within iPython, or with AWS Cloud instances. After installing the QAP software package, these scripts can be run from any directory:

* qap_anatomical_spatial.py
* qap_functional_spatial.py
* qap_functional_temporal.py

For command-line runs:

    qap_anatomical_spatial.py --sublist {path to subject list YAML file} {path to pipeline configuration YAML file}

Executing either of the scripts with only the *-h* flag will produce a short help manual listing the command line arguments.

##Running the QAP Pipelines on AWS Amazon Cloud Instances

*** STILL IN PROGRESS ***

With access to the Amazon Cloud, the QAP measures can be calculated for a large volume of subjects quickly.

[AMI setup.. documentation necessary here? or link to somewhere else]

### install QAP on your AMI, and pre-requisites

### Generating Your S3 Subject Dictionary File

The QAP software package comes with a script called *qap_aws_s3_dict_generator.py*, which can be run from any directory once the package is installed. This script will create a YAML file containing the filepaths to your data stored on your AWS S3 bucket storage. You will need this dictionary YAML file to start an AWS Cloud run for QAP. This script takes in five input parameters:

* **scan_type**: *anat* or *func*, depending on which QAP measures you will be using the S3 subject dictionary for
* **bucket_name**: the name of your AWS S3 bucket
* **bucket_prefix**:  the filepath prefix to the top level of your raw data directory on S3 storage

For example, if your S3 storage is arranged like so:

	/data/project/raw_data/sub001/session_1/scan_1/file.nii.gz
	/data/project/raw_data/sub001/session_1/scan_2/file.nii.gz
	/data/project/raw_data/sub002/session_1/scan_1/file.nii.gz
                    
Then the bucket_prefix would be:

	/data/project/raw_data


* **creds_path**: the path to the file containing your AWS credentials
* **outfile_path**: the full filepath for the S3 subject YAML dictionary this script will create

Once this script is run, it will output the S3 dictionary YAML file, and it will give you the total number of subject-session-scans. Take note of this number, because you will need to list it in your SGE batch file (more below).

### Setting Up Your SGE File

	#! /bin/bash
	#$ -cwd
	#$ -S /bin/bash
	#$ -V
	#$ -t 1-{size of sub dictionary}
	#$ -q all.q
	#$ -pe mpi_smp 4
	#$ -e /home/ubuntu/qap_run_anat.err
	#$ -o /home/ubuntu/qap_run_anat.out
	source /etc/profile.d/cpac_env.sh
	ANAT_S3_DICT=/home/ubuntu/configs/anat_s3_dict.yml
	ANAT_SP_CONFIG_FILE=/home/ubuntu/configs/qap_config_anat_spatial.yml
	echo "Start - TASKID " $SGE_TASK_ID " : " $(date)
	# Run anatomical spatial qap
	qap_anatomical_spatial.py --subj_idx $SGE_TASK_ID --s3_dict_yml $ANAT_S3_DICT $ANAT_SP_CONFIG_FILE
	echo "End - TASKID " $SGE_TASK_ID " : " $(date)

## References

[^1]: Magnotta, V. A., & Friedman, L. (2006). Measurement of signal-to-noise and contrast-to-noise in the fBIRN multicenter imaging study. Journal of Digital Imaging, 19(2), 140-147.

[^2]: Atkinson D, Hill DL, Stoyle PN, Summers PE, Keevil SF (1997). Automatic correction of motion artifacts in magnetic resonance images using an entropy focus criterion. IEEE Trans Med Imaging. 16(6):903-10.

[^3]: Friedman, L., Stern, H., Brown, G. G., Mathalon, D. H., Turner, J., Glover, G. H., ... & Potkin, S. G. (2008). Test–retest and between‐site reliability in a multicenter fMRI study. Human brain mapping, 29(8), 958-972.

[^4]: Mortamet, B., Bernstein, M. A., Jack, C. R., Gunter, J. L., Ward, C., Britson, P. J., ... & Krueger, G. (2009). Automatic quality assessment in structural brain magnetic resonance imaging. Magnetic Resonance in Medicine, 62(2), 365-372.

[^5]: Power, J. D., Barnes, K. A., Snyder, A. Z., Schlaggar, B. L. & Petersen, S. E. (2012) Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage 59, 2142-2154.

[^6]: Nichols, T. (2012, Oct 28). Standardizing DVARS. Retrieved from http://blogs.warwick.ac.uk/nichols/entry/standardizing_dvars.

[^7]: Cox, R.W. (1996) AFNI: Software for analysis and visualization of functional magnetic resonance neuroimages. Computers and Biomedical Research, 29:162-173.

[^8]: Jenkinson, M., Bannister, P., Brady, M., & Smith, S. (2002). Improved optimization for the robust and accurate linear registration and motion correction of brain images. Neuroimage, 17(2), 825-841.

[^9]: Giannelli, M., Diciotti, S., Tessa, C., & Mascalchi, M. (2010). Characterization of Nyquist ghost in EPI-fMRI acquisition sequences implemented on two clinical 1.5 T MR scanner systems: effect of readout bandwidth and echo spacing. Journal of Applied Clinical Medical Physics, 11(4).

[^10]: Yan CG, Cheung B, Kelly C, Colcombe S, Craddock RC, Di Martino A, Li Q, Zuo XN, Castellanos FX, Milham MP (2013). A comprehensive assessment of regional variation in the impact of head micromovements on functional connectomics. Neuroimage. 76:183-201.
