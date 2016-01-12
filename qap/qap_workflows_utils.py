

def run_3dClipLevel(input_skull):

    import commands

    cmd = "3dClipLevel %s" % input_skull

    thresh_out = int(float((commands.getoutput(cmd).split("\n")[1])))


    return thresh_out



def slice_head_mask(infile, transform, standard):

    import os
    import sys

    import nibabel as nb
    import numpy as np
    import subprocess
    import pkg_resources as p

    # get file info
    infile_img = nb.load(infile)
    infile_header = infile_img.get_header()
    infile_affine = infile_img.get_affine()

    infile_dims = infile_header.get_data_shape()

    # these are stored in the files listed below, just here for reference
    inpoint_a = "78 -110 -72"
    inpoint_b = "-78 -110 -72"
    inpoint_c = "0 88 -72"  # nose, apparently

    # these each contain a set of coordinates for drawing the plane across
    # the image (to "slice" it)
    inpoint_files = [p.resource_filename("qap", "inpoint_a.txt"),
                     p.resource_filename("qap", "inpoint_b.txt"),
                     p.resource_filename("qap", "inpoint_c.txt")]

    # let's convert the coordinates into voxel coordinates

    coords = []

    for inpoint in inpoint_files:

        coord_cmd = "std2imgcoord -img %s -std %s -xfm %s -vox %s" \
                    % (infile, standard, transform, inpoint)

        coord_out = subprocess.check_output(coord_cmd, shell=True)

        if "Could not" in coord_out:
            raise Exception(coord_out)

        coords.append(coord_out)

    # get the converted coordinates into a list format, and also check to make
    # sure they are not "out of bounds"
    new_coords = []

    for coord in coords:

        co_nums = coord.split(" ")

        co_nums_newlist = []

        for num in co_nums:

            if num != "":
                co_nums_newlist.append(int(num.split(".")[0]))

        for ind in range(0, 3):

            if co_nums_newlist[ind] > infile_dims[ind]:
                co_nums_newlist[ind] = infile_dims[ind]

            elif co_nums_newlist[ind] < 1:
                co_nums_newlist[ind] = 1

        new_coords.append(co_nums_newlist)

    # get the vectors connecting the points

    u = []

    for a_pt, c_pt in zip(new_coords[0], new_coords[2]):

        u.append(int(a_pt - c_pt))

    v = []

    for b_pt, c_pt in zip(new_coords[1], new_coords[2]):

        v.append(int(b_pt - c_pt))

    u_vector = np.asarray(u)
    v_vector = np.asarray(v)

    # vector cross product
    n = np.cross(u, v)

    # normalize the vector
    n = n / np.linalg.norm(n, 2)

    constant = np.dot(n, np.asarray(new_coords[0]))

    # now determine the z-coordinate for each pair of x,y
    plane_dict = {}

    for yvox in range(0, infile_dims[1]):

        for xvox in range(0, infile_dims[0]):

            zvox = (constant - (n[0] * xvox + n[1] * yvox)) / n[2]

            zvox = np.floor(zvox)

            if zvox < 1:
                zvox = 1
            elif zvox > infile_dims[2]:
                zvox = infile_dims[2]

            plane_dict[(xvox, yvox)] = zvox

    # create the mask
    mask_array = np.zeros(infile_dims)

    for x in range(0, infile_dims[0]):

        for y in range(0, infile_dims[1]):

            for z in range(0, infile_dims[2]):

                if plane_dict[(x, y)] > z:
                    mask_array[x, y, z] = 1

    new_mask_img = nb.Nifti1Image(mask_array, infile_affine, infile_header)

    infile_filename = infile.split("/")[-1].split(".")[0]

    outfile_name = infile_filename + "_slice_mask.nii.gz"
    outfile_path = os.path.join(os.getcwd(), outfile_name)

    nb.save(new_mask_img, outfile_path)

    return outfile_path



def qap_anatomical_spatial(anatomical_reorient, head_mask_path,
                           anatomical_gm_mask, anatomical_wm_mask,
                           anatomical_csf_mask, subject_id, session_id,
                           scan_id, site_name=None, out_vox=True):

    import os
    import sys

    from qap.spatial_qc import summary_mask, snr, cnr, fber, efc, \
        artifacts, fwhm
    from qap.qap_utils import load_image, load_mask

    # Load the data
    anat_data = load_image(anatomical_reorient)
    fg_mask = load_mask(head_mask_path, anatomical_reorient)
    bg_mask = 1 - fg_mask

    gm_mask = load_mask(anatomical_gm_mask, anatomical_reorient)
    wm_mask = load_mask(anatomical_wm_mask, anatomical_reorient)
    csf_mask = load_mask(anatomical_csf_mask, anatomical_reorient)

    # Initialize QC
    qc = dict()

    qc['Participant'] = subject_id

    qc['Session'] = session_id

    qc['Series'] = scan_id

    if site_name:
        qc['Site'] = site_name

    # FBER
    qc['FBER'] = fber(anat_data, fg_mask)

    # EFC
    qc['EFC'] = efc(anat_data)

    # Artifact
    qc['Qi1'], _ = artifacts(anat_data, fg_mask, calculate_qi2=False)

    # Smoothness in voxels
    tmp = fwhm(anatomical_reorient, head_mask_path, out_vox=out_vox)
    qc['FWHM_x'], qc['FWHM_y'], qc['FWHM_z'], qc['FWHM'] = tmp

    # Summary Measures
    fg_mean, fg_std, fg_size = summary_mask(anat_data, fg_mask)
    bg_mean, bg_std, bg_size = summary_mask(anat_data, bg_mask)

    gm_mean, gm_std, gm_size = (None, None, None)
    wm_mean, wm_std, wm_size = (None, None, None)
    csf_mean, csf_std, csf_size = (None, None, None)
    qc['CNR'] = None
    qc['SNR'] = None

    # More Summary Measures
    gm_mean, gm_std, gm_size = summary_mask(anat_data, gm_mask)
    wm_mean, wm_std, wm_size = summary_mask(anat_data, wm_mask)
    csf_mean, csf_std, csf_size = summary_mask(anat_data, csf_mask)

    # SNR
    qc['SNR'] = snr(fg_mean, bg_std)

    # CNR
    qc['CNR'] = cnr(gm_mean, wm_mean, bg_std)


    return qc



def qap_functional_spatial(mean_epi, func_brain_mask, direction, subject_id,
                           session_id, scan_id, site_name=None, out_vox=True):

    import os
    import sys

    from qap.spatial_qc import summary_mask, snr, fber, efc, fwhm, \
        ghost_direction
    from qap.qap_utils import load_image, load_mask

    # Load the data
    anat_data = load_image(mean_epi)
    fg_mask = load_mask(func_brain_mask, mean_epi)
    bg_mask = 1 - fg_mask

    # Initialize QC
    qc = dict(Participant=subject_id, Session=session_id, Series=scan_id)

    if site_name:
        qc['Site'] = site_name

    # FBER
    qc['FBER'] = fber(anat_data, fg_mask)

    # EFC
    qc['EFC'] = efc(anat_data)

    # Smoothness in voxels
    tmp = fwhm(mean_epi, func_brain_mask, out_vox=out_vox)
    qc['FWHM_x'], qc['FWHM_y'], qc['FWHM_z'], qc['FWHM'] = tmp

    # Ghosting
    if (direction == "all"):
        qc['Ghost_x'] = ghost_direction(anat_data, fg_mask, "x")
        qc['Ghost_y'] = ghost_direction(anat_data, fg_mask, "y")
        qc['Ghost_z'] = ghost_direction(anat_data, fg_mask, "z")

    else:
        qc['Ghost_%s' % direction] = ghost_direction(anat_data, fg_mask,
                                                     direction)

    # Summary Measures
    fg_mean, fg_std, fg_size = summary_mask(anat_data, fg_mask)
    bg_mean, bg_std, bg_size = summary_mask(anat_data, bg_mask)

    qc['SNR'] = None

    # SNR
    qc['SNR'] = snr(fg_mean, bg_std)


    return qc



def qap_functional_temporal(
        func_timeseries, func_brain_mask, tsnr_volume, fd_file,
        subject_id, session_id, scan_id, site_name=None, 
        motion_threshold=1.0):

    import sys
    import nibabel as nb
    import numpy as np

    from qap.temporal_qc import mean_dvars_wrapper, mean_outlier_timepoints, \
        mean_quality_timepoints, global_correlation, \
        calculate_percent_outliers

    # DVARS
    mean_dvars = mean_dvars_wrapper(func_timeseries, func_brain_mask)

    # Mean FD (Jenkinson)
    fd = np.loadtxt(fd_file)

    meanfd_outliers, meanfd_iqr = calculate_percent_outliers(fd)

    # Calculate Outliers
    # Number and Percent of frames (time points) where
    # movement (FD) exceeded threshold
    #num_fd = np.float((fd > motion_threshold).sum())
    #percent_fd = (num_fd * 100) / (len(fd) + 1)

    # 3dTout
    mean_outlier, outlier_perc_out, outlier_IQR = \
        mean_outlier_timepoints(func_timeseries)

    # 3dTqual
    mean_quality, qual_perc_out, qual_IQR = \
        mean_quality_timepoints(func_timeseries)

    # GCOR
    gcor = global_correlation(func_timeseries, func_brain_mask)

    # Compile
    qc = {
        "Participant":   subject_id,
        "Session":   session_id,
        "Series":      scan_id,
        "Std. DVARS":     mean_dvars,
        "FD (Mean)":   fd.mean(),
        "FD percent outliers": meanfd_outliers,
        "FD IQR": meanfd_iqr,
        "Outliers Values (Mean)": mean_outlier,
        "Outliers Values percent outliers": outlier_perc_out,
        "Outliers Values IQR": outlier_IQR,
        "Quality":   mean_quality,
        "Quality percent outliers": qual_perc_out,
        "Quality IQR": qual_IQR,
        "GCOR":      gcor
    }

    if site_name:
        qc['Site'] = site_name


    return qc
