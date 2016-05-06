

def create_expr_string(clip_level_value):

    expr_string = "step(a-%s)" % clip_level_value

    return expr_string



def slice_head_mask(infile, transform, standard):

    import os
    import sys

    import nibabel as nb
    import numpy as np
    import numpy.linalg as npl
    import subprocess
    import pkg_resources as p

    # get file info
    infile_img = nb.load(infile)
    infile_header = infile_img.get_header()
    infile_affine = infile_img.get_affine()

    infile_dims = infile_header.get_data_shape()

    # these coordinates correspond to the points of defining the slice plane
    # on the MNI template
    inpoint_a = [78, -110, -72, 0]
    inpoint_b = [-78, -110, -72, 0]
    inpoint_c = [-1, 91, -29, 0]  # nose

    inpoint_coords = [inpoint_a, inpoint_b, inpoint_c]

    # get the affine output matrix of 3dallineate
    with open(transform,"r") as f:
        allineate_mat_list = f.readlines()

    # get the 3dAllineate output affine matrix into a list
    mat_list = filter(None,allineate_mat_list[1].rstrip("\n").split(" "))

    # put together the 4x4 matrix
    #   (encode the offset as the fourth row and include the 0,0,0,1
    #    dummy row, so that this will only require one matrix multiplication)
    row1 = [float(mat_list[0]), float(mat_list[1]), float(mat_list[2]), \
                float(mat_list[3])]
    row2 = [float(mat_list[4]), float(mat_list[5]), float(mat_list[6]), \
                float(mat_list[7])]
    row3 = [float(mat_list[8]), float(mat_list[9]), float(mat_list[10]), \
                float(mat_list[11])]
    row4 = [0,0,0,1]

    allineate_mat = np.asarray([row1,row2,row3,row4])

    coords_list = []

    for inpoint in inpoint_coords:      

        # using the transform, calculate what the three coordinates are in
        # native space, as it corresponds to the anatomical scan
        coord_out = \
            list(np.dot(np.linalg.inv(allineate_mat),inpoint))

        # remove the one resulting zero at the end
        coord_out = coord_out[:-1]

        # convert the coordinates from mm to voxels
        coord_out = \
            nb.affines.apply_affine(npl.inv(infile_affine),coord_out)

        coords_list.append(coord_out)


    # make sure converted coordinates are not "out of bounds"
    new_coords = []

    for coords in coords_list:

        co_nums_newlist = []

        for num in coords:

            if num != "":
                co_nums_newlist.append(int(num)) #.split(".")[0]))

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

    try:
        nb.save(new_mask_img, outfile_path)
    except:
        raise_smart_exception(locals())


    return outfile_path



def qap_anatomical_spatial(anatomical_reorient, qap_head_mask_path,
                           whole_head_mask_path, skull_mask_path,
                           anatomical_gm_mask, anatomical_wm_mask,
                           anatomical_csf_mask, subject_id, session_id,
                           scan_id, site_name=None, out_vox=True,
                           starter=None):

    import os
    import sys

    from qap.spatial_qc import summary_mask, snr, cnr, fber, efc, \
        artifacts, fwhm
    from qap.qap_utils import load_image, load_mask

    # Load the data
    anat_data = load_image(anatomical_reorient)

    fg_mask = load_mask(qap_head_mask_path, anatomical_reorient)
    bg_mask = 1 - fg_mask

    whole_head_mask = load_mask(whole_head_mask_path, anatomical_reorient)
    skull_mask = load_mask(skull_mask_path, anatomical_reorient)

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
    qc['FBER'] = fber(anat_data, skull_mask, bg_mask)

    # EFC
    qc['EFC'] = efc(anat_data)

    # Artifact
    qc['Qi1'], _ = artifacts(anat_data, fg_mask, calculate_qi2=False)

    # Smoothness in voxels
    tmp = fwhm(anatomical_reorient, whole_head_mask_path, out_vox=out_vox)
    qc['FWHM_x'], qc['FWHM_y'], qc['FWHM_z'], qc['FWHM'] = tmp

    # Summary Measures
    fg_mean, fg_std, fg_size = summary_mask(anat_data, whole_head_mask)
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
                           session_id, scan_id, site_name=None, out_vox=True,
                           starter=None):

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
    qc['FBER'] = fber(anat_data, fg_mask, bg_mask)

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
        func_timeseries, func_brain_mask, fd_file,
        subject_id, session_id, scan_id, site_name=None, 
        motion_threshold=1.0, starter=None):

    import sys
    import nibabel as nb
    import numpy as np

    from qap.temporal_qc import outlier_timepoints, quality_timepoints, \
                                global_correlation, calculate_percent_outliers
    from qap.dvars import calc_dvars

    # DVARS
    dvars = calc_dvars(func_timeseries, func_brain_mask)
    dvars_outliers, dvars_IQR = calculate_percent_outliers(dvars)

    mean_dvars = dvars.mean(0)
    mean_dvars = mean_dvars[0]

    # Mean FD (Jenkinson)
    fd = np.loadtxt(fd_file)
    meanfd_outliers, meanfd_IQR = calculate_percent_outliers(fd)

    # 3dTout
    outliers = outlier_timepoints(func_timeseries)
    # calculate the outliers of the outliers! AAHH!
    outlier_perc_out, outlier_IQR = calculate_percent_outliers(outliers)

    # 3dTqual
    quality = quality_timepoints(func_timeseries)
    quality_outliers, quality_IQR = calculate_percent_outliers(quality)

    # GCOR
    gcor = global_correlation(func_timeseries, func_brain_mask)

    # Compile
    qc = {
        "Participant": subject_id,
        "Session": session_id,
        "Series": scan_id,
        "Std. DVARS (Mean)": mean_dvars,
        "Std. DVARS (Std Dev)": np.std(dvars),
        "Std. DVARS (Median)": np.median(dvars),
        "Std. DVARs IQR": dvars_IQR,
        "Std. DVARS percent outliers": dvars_outliers,
        "RMSD (Mean)": np.mean(fd),
        "RMSD (Std Dev)": np.std(fd),
        "RMSD (Median)": np.median(fd),
        "RMSD IQR": meanfd_IQR,
        "RMSD percent outliers": meanfd_outliers,
        "Fraction of Outliers (Mean)": np.mean(outliers),
        "Fraction of Outliers (Std Dev)": np.std(outliers),
        "Fraction of Outliers (Median)": np.median(outliers),
        "Fraction of Outliers IQR": outlier_IQR,
        "Fraction of Outliers percent outliers": outlier_perc_out,
        "Quality (Mean)": np.mean(quality),
        "Quality (Std Dev)": np.std(quality),
        "Quality (Median)": np.median(quality),
        "Quality IQR": quality_IQR,
        "Quality percent outliers": quality_outliers,
        "GCOR": gcor
    }

    if site_name:
        qc['Site'] = site_name


    return qc
