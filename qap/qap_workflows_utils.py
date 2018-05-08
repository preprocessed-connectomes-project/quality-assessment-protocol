
"""Calculates the global signal time series of a functional file
   by averaging all voxels in each volume

    :type functional_file: nifti file
    :param functional_file: Functional file
    :rtype: float array
    :return: global signal time series 
    """
def global_signal_time_series(functional_file):
    import numpy as np
    import nibabel as nb

    func = nb.load(functional_file).get_data()
    time = func.shape[-1]
    output = [0]*time

    # that is a stupid way of calculating the mean of each volume,
    # maybe numpy can do it in one line
    for i in range(time):
        output[i] = func[:,:,:,i].mean()
    return output


def convert_allineate_xfm(mat_list):
    """Convert the affine transform output of AFNI's 3dAllineate into an
    equivalent 4x4 matrix.

    - Takes the 3x3 + offset format of the input matrix and turns it into a
      4x4 with the last row being 0,0,0,1.

    :type mat_list: list
    :param mat_list: A vector (list) of the flattened affine matrix.
    :rtype: NumPy array
    :return: A NumPy array of the converted affine matrix.
    """

    import numpy as np

    # put together the 4x4 matrix
    #   (encode the offset as the fourth row and include the 0,0,0,1
    #    dummy row, so that this will only require one matrix multiplication)
    row1 = [float(mat_list[0]), float(mat_list[1]), float(mat_list[2]),
                float(mat_list[3])]
    row2 = [float(mat_list[4]), float(mat_list[5]), float(mat_list[6]),
                float(mat_list[7])]
    row3 = [float(mat_list[8]), float(mat_list[9]), float(mat_list[10]),
                float(mat_list[11])]
    row4 = [0, 0, 0, 1]

    allineate_mat = np.asarray([row1, row2, row3, row4])

    return allineate_mat


def warp_coordinates(inpoint, allineate_mat, infile_affine, infile_dims):
    """Warp spatial coordinates using a 4x4 affine matrix.

    :type inpoint: list
    :param inpoint: A list of three numbers describing a coordinate in 3D
                    space.
    :type allineate_mat: NumPy array
    :param allineate_mat: A NumPy array describing the 4x4 affine transform.
    :type infile_affine: Nibabel affine
    :param infile_affine: A Nibabel NIFTI affine info object.
    :type infile_dims: tuple
    :param infile_dims: A tuple containing the dimensions of the NIFTI file
                        the coordinates pertain to.
    :rtype: list
    :return: A list of three warped coordinates.
    """

    import numpy as np
    import numpy.linalg as npl
    import nibabel as nb

    # using the transform, calculate what the three coordinates are in
    # native space, as it corresponds to the anatomical scan
    coord_out = \
        list(np.dot(np.linalg.inv(allineate_mat), inpoint))

    # remove the one resulting zero at the end
    coord_out = coord_out[:-1]

    # convert the coordinates from mm to voxels
    coord_out = \
        nb.affines.apply_affine(npl.inv(infile_affine), coord_out)

    # make sure converted coordinates are not "out of bounds"
    co_nums_newlist = []

    for num in coord_out:
        if num != "":
            co_nums_newlist.append(int(num))

    for ind in range(0, 3):

        if co_nums_newlist[ind] > infile_dims[ind]:
            co_nums_newlist[ind] = infile_dims[ind]

        elif co_nums_newlist[ind] < 1:
            co_nums_newlist[ind] = 1

    return co_nums_newlist


def calculate_plane_coords(coords_list, infile_dims):
    """Calculate the coordinates of the triangular plane spanning from roughly
    around the participant's nose, down to roughly below the rear of the neck.

    :type coords_list: list
    :param coords_list: A list of lists, describing the three coordinate
                        points of the triangular plane.
    :type infile_dims: list
    :param infile_dims: A list of the NIFTI file's dimensions.
    :rtype: dict
    :return: A dictionary mapping the z coordinates to the (x,y) coordinate
             pairs.
    """

    import numpy as np

    # get the vectors connecting the points
    u = []
    for a_pt, c_pt in zip(coords_list[0], coords_list[2]):
        u.append(int(a_pt - c_pt))

    v = []
    for b_pt, c_pt in zip(coords_list[1], coords_list[2]):
        v.append(int(b_pt - c_pt))

    # vector cross product
    n = np.cross(u, v)

    # normalize the vector
    n = n / np.linalg.norm(n, 2)
    constant = np.dot(n, np.asarray(coords_list[0]))

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

    return plane_dict


def scipy_fill(in_file):

    import os
    import nibabel as nb
    from scipy.ndimage.morphology import binary_fill_holes

    mask_img = nb.load(in_file)
    mask_data = mask_img.get_data()

    new_mask_data = binary_fill_holes(mask_data).astype('int8')

    new_mask_img = nb.Nifti1Image(new_mask_data, mask_img.affine)

    out_file = os.path.join(os.getcwd(), "qap_headmask_filled.nii.gz")
    new_mask_img.to_filename(out_file)

    return out_file


def create_slice_mask(plane_dict, infile_dims):
    """Create a binary array defining a mask covering the area below a given
    plane in the 3D image.

    :type plane_dict: dict
    :param plane_dict: A dictionary matching z voxel coordinates to
                       corresponding (x,y) coordinates.
    :type infile_dims: list
    :param infile_dims: A list of the NIFTI file's dimensions.
    :rtype: NumPy array
    :return: A NumPy array defining the binary mask of the slice.
    """

    import numpy as np

    mask_array = np.zeros(infile_dims)

    for x in range(0, infile_dims[0]):

        for y in range(0, infile_dims[1]):

            for z in range(0, infile_dims[2]):

                if plane_dict[(x, y)] > z:
                    mask_array[x, y, z] = 1

    return mask_array


def slice_head_mask(infile, transform):
    """Write out a binary mask NIFTI image defining a triangular area covering
    the region below the nose and mouth.

    :type infile: str
    :param infile: Filepath to the participant's anatomical scan.
    :type transform: str
    :param transform: Filepath to the text file containing the affine matrix
                      output of AFNI's 3dAllineate describing the warp from
                      the anatomical scan to a template.
    :rtype: str
    :return: Filepath to the new head mask NIFTI file.
    """

    import os
    import nibabel as nb

    from qap.script_utils import read_txt_file
    from qap.qap_workflows_utils import convert_allineate_xfm, \
                                        warp_coordinates, \
                                        calculate_plane_coords, \
                                        create_slice_mask
    from qap.qap_utils import read_nifti_image, write_nifti_image

    # get file info
    infile_img = read_nifti_image(infile)

    infile_header = infile_img.get_header()
    infile_affine = infile_img.get_affine()
    infile_dims = infile_header.get_data_shape()

    # get the affine output matrix of 3dallineate
    allineate_mat_list = read_txt_file(transform)

    # get the 3dAllineate output affine matrix into a list
    mat_list = filter(None, allineate_mat_list[1].rstrip("\n").split(" "))

    # convert 3dAllineate's affine transform
    allineate_mat = convert_allineate_xfm(mat_list)

    # get participant-specific coordinates
    #   these coordinates correspond to the points of defining the slice plane
    #   on the MNI template
    inpoint_a = [78, -110, -72, 0]
    inpoint_b = [-78, -110, -72, 0]
    inpoint_c = [-1, 91, -29, 0]  # nose
    inpoint_coords = [inpoint_a, inpoint_b, inpoint_c]

    coords_list = map(warp_coordinates, inpoint_coords, 
        [allineate_mat] * 3, [infile_affine] * 3, [infile_dims] * 3)

    # calculate normalized vector and get z coordinate for each x,y pair
    plane_dict = calculate_plane_coords(coords_list, infile_dims)

    # create the mask
    mask_array = create_slice_mask(plane_dict, infile_dims)

    # create new slice mask img file
    new_mask_img = nb.Nifti1Image(mask_array, infile_affine, infile_header)

    infile_filename = infile.split("/")[-1].split(".")[0]

    outfile_name = "_".join([infile_filename, "slice_mask.nii.gz"])
    outfile_path = os.path.join(os.getcwd(), outfile_name)

    write_nifti_image(new_mask_img, outfile_path)

    return outfile_path


def calc_temporal_std(voxel_ts):
    '''this can be used in a map in the below function later when you move to
    optimize it
    '''
    import numpy as np
    return np.std(voxel_ts)


def get_temporal_std_map(func_reorient, func_mask):
    """Create a map of the standard deviations of each voxel's timeseries.

    :param func_reorient: String filepath of the deobliqued, reoriented
                          functional timeseries NIFTI file.
    :param func_mask: String filepath of the functional brain mask NIFTI file.
    :return: String filepath of the temporal standard deviation map NIFTI
             file.
    """

    import os
    import numpy as np
    import nibabel as nb
    from qap.qap_utils import get_masked_data

    func_data = get_masked_data(func_reorient, func_mask)
    temporal_std_map = np.std(func_data, axis=-1)

    # write the image
    mask_img = nb.load(func_mask)
    tstd_img = nb.Nifti1Image(temporal_std_map, mask_img.affine)

    # this writes it into the node's folder in the working directory
    temporal_std_map_file = os.path.join(os.getcwd(), "tstd.nii.gz")
    tstd_img.to_filename(temporal_std_map_file)

    return temporal_std_map_file


def create_threshold_mask(data, threshold):
    """Create a binary mask of a dataset based on a given threshold value.

    :param data: Numpy array of a dataset.
    :param threshold: The value to threshold the dataset with in order to
                      define the binary mask.
    :return: Numpy array of the binary mask values for each voxel.
    """

    import numpy as np
    mask = np.zeros(data.shape)
    mask[data > threshold] = 1
    return mask


def calc_estimated_csf_nuisance(temporal_std_map):
    """Calculate the estimated CSF nuisance using a map of the temporal
    standard deviation.

    :param temporal_std_map: Numpy array of a map of the standard deviations
                             of each voxel's timeseries.
    :return: Numpy array of the map of estimated CSF nuisance for each voxel.
    """

    import numpy as np
    from qap.qap_utils import get_masked_data

    all_tstd = np.asarray(temporal_std_map.nonzero()).flatten()
    all_tstd_sorted = sorted(all_tstd)
    top_2 = 0.98 * len(all_tstd)

    top_2_std = all_tstd_sorted[int(top_2):]
    cutoff = top_2_std[0]

    estimated_nuisance_mask = create_threshold_mask(temporal_std_map, cutoff)

    nuisance_stds = get_masked_data(temporal_std_map, estimated_nuisance_mask)

    return nuisance_stds


def calc_compcor_nuisance(func_mean):
    """Calculate the estimated Compcor nuisance using the tSTD based
    determination of noise ROI.

    :param temporal_std_map: Numpy array of a map of the standard deviations
                             of each voxel's timeseries.
    :return: Numpy array of the map of estimated CSF nuisance for each voxel.
    """

    import numpy as np

    timepoints = func_mean.shape[-1]

    data = func_mean.reshape((-1, timepoints))
    X = np.ones((timepoints, 1))
    betas = np.linalg.pinv(X).dot(data.T)
    data_hat = X.dot(betas).T
    
    data = (data - data_hat).reshape(func_mean.shape).T
    std = np.std(data, axis=0)
    std[std == 0] = 1.
    std[np.isnan(std)] = 1.

    data = data / std
    u, _, _ = np.linalg.svd(data, full_matrices=False)

    return u.T


def sfs_voxel(voxel_ts, total_func_mean, voxel_ts_std, nuisance_mean_std):
    """Calculate the Signal Fluctuation Intensity (SFS) of one voxel's
    functional time series.

    - From "Signal Fluctuation Sensitivity: An Improved Metric for Optimizing
      Detection of Resting-State fMRI Networks", Daniel J. DeDora1,
      Sanja Nedic, Pratha Katti, Shafique Arnab, Lawrence L. Wald, Atsushi
      Takahashi, Koene R. A. Van Dijk, Helmut H. Strey and
      Lilianne R. Mujica-Parodi. More info here:
        http://journal.frontiersin.org/article/10.3389/fnins.2016.00180/full

    :param: Tuple of arguments, containing 1.) the timeseries of the current
            voxel, 2.) the mean value of all of the voxel timeseries, 3.) the
            standard deviation of the current voxel timeseries, and 4.) the
            mean standard deviation of the estimated CSF nuisance.
    :return: Numpy array of the signal fluctuation intensity timecourse for
             the voxel timeseries provided.
    """

    sfs_vox = \
        (voxel_ts/total_func_mean) * (voxel_ts_std/nuisance_mean_std)

    return sfs_vox


def sfs_timeseries(func_mean, func_mask, temporal_std_file):
    """Average the SFS timecourses of each voxel into one SFS timeseries.

    :param func_mean: String filepath to the mean functional NIFTI file.
    :param func_mask: String filepath to the functional brain mask NIFTI file.
    :param temporal_std_file: String filepath to the temporal standard
                              deviation map NIFTI file.
    :return: The string filepaths of the SFS map file, and the estimated
             nuisance map file.
    """

    import os
    import numpy as np
    import nibabel as nb
    from qap.qap_workflows_utils import calc_compcor_nuisance, sfs_voxel

    func_mean_img = nb.load(func_mean)
    func_mean_data = func_mean_img.get_data()
    tstd_img = nb.load(temporal_std_file)
    temporal_std_map = tstd_img.get_data()

    total_func_mean = np.mean(func_mean_data)
    nuisance_stds = calc_compcor_nuisance(func_mean_data)
    nuisance_mean_std = np.mean(np.asarray(nuisance_stds.nonzero()).flatten())

    sfs_voxels = np.asarray([   
        sfs_voxel(voxel_ts, total_func_mean, voxel_std, nuisance_mean_std)
        for voxel_ts, voxel_std in zip(func_mean_data, temporal_std_map)
    ])

    # write the images
    mask_img = nb.load(func_mask)
    est_n_img = nb.Nifti1Image(nuisance_stds, mask_img.affine)
    # this writes it into the node's folder in the working directory
    est_nuisance_file = os.path.join(os.getcwd(), "estimated_nuisance.nii.gz")
    est_n_img.to_filename(est_nuisance_file)

    sfs_img = nb.Nifti1Image(sfs_voxels, mask_img.affine)
    # this writes it into the node's folder in the working directory
    sfs_file = os.path.join(os.getcwd(), "SFS.nii.gz")
    sfs_img.to_filename(sfs_file)

    return sfs_file, est_nuisance_file


def create_header_dict_entry(in_file, subject, session, scan, type):
    """Gather the header information from a NIFTI file and arrange it into a
    Python dictionary.

    :type in_file: str
    :param in_file: Filepath to the NIFTI raw data scan.
    :type subject: str
    :param subject: The participant ID.
    :type session: str
    :param session: The session ID.
    :type scan: str
    :param scan: The scan ID.
    :type type: str
    :param type: The data type ("anatomical" or "functional").
    :rtype: dict
    :return: A dictionary with the header information of the file, keyed by
             the participant's ID data.
    """

    if not os.path.isfile(in_file):
        err = "Filepath doesn't exist!\nFilepath: %s" % in_file
        raise Exception(err)

    try:
        img = nb.load(in_file)
        img_header = img.header
    except:
        err = "You may not have an up-to-date installation of the Python " \
              "Nibabel package.\nYour Nibabel version: %s" % \
              str(nb.__version__)
        raise Exception(err)

    subkey = "%s_header_info" % type
    id_string = "%s %s %s" % (subject, session, scan)
    qap_dict = {id_string: {subkey: {}}}

    info_labels = ["descrip", "db_name", "bitpix", "slice_start",
                   "scl_slope", "scl_inter", "slice_end", "slice_duration",
                   "toffset", "quatern_b", "quatern_c", "quatern_d",
                   "qoffset_x", "qoffset_y", "qoffset_z", "aux_file",
                   "intent_name", "slice_code", "data_type", "qform_code",
                   "srow_x", "srow_y", "srow_z", "sform_code"]

    for info_label in info_labels:
        try:
            qap_dict[id_string][subkey][info_label] = \
                str(img_header[info_label])
        except KeyError:
            print ("\n\n{0} field not in NIFTI header of {1}"
                   "\n\n".format(info_label, in_file))
            qap_dict[id_string][subkey][info_label] = ""
            pass

    try:
        pixdim = img_header['pixdim']
        qap_dict[id_string][subkey]["pix_dimx"] = str(pixdim[1])
        qap_dict[id_string][subkey]["pix_dimy"] = str(pixdim[2])
        qap_dict[id_string][subkey]["pix_dimz"] = str(pixdim[3])
        qap_dict[id_string][subkey]["tr"] = str(pixdim[4])
    except KeyError:
        print("\n\npix_dim/TR fields not in NIFTI header of {0}"
              "\n\n".format(in_file))
        pass

    try:
        qap_dict[id_string][subkey]["extensions"] = \
            len(img.header.extensions.get_codes())
    except:
        print("\n\nExtensions not in NIFTI header of {0}\n\n".format(in_file))
        pass

    return qap_dict


def qap_anatomical_spatial(anatomical_reorient, qap_head_mask_path,
                           qap_bg_head_mask_path, whole_head_mask_path,
                           skull_mask_path, anatomical_gm_mask,
                           anatomical_wm_mask, anatomical_csf_mask,
                           fav_artifacts, subject_id, session_id, scan_id,
                           run_name, site_name=None, exclude_zeroes=False,
                           out_vox=True, session_output_dir=None,
                           starter=None):
    """Calculate the anatomical spatial QAP measures for an anatomical scan.

    - The exclude_zeroes flag is useful for when a large amount of zero
      values have been artificially injected into the image, for example,
      when removing the faces and ears in scans for privacy compliance
      reasons; these zeroes can artificially skew the spatial quality metric
      results and make it seem that there is far less noise or artifacts in
      the image than there really is.
    - The inclusion of the starter node allows several QAP measure pipelines
      which are not dependent on one another to be executed as one pipeline.
      This allows the MultiProc Nipype plugin to efficiently manage
      resources when parallelizing.

    :type anatomical_reorient: str
    :param anatomical_reorient: Filepath to the reoriented anatomical scan.
    :type qap_head_mask_path: str
    :param qap_head_mask_path: Filepath to mask of the head, plus the slice
                               covering the region below the nose and in
                               front of the mouth.
    :type whole_head_mask_path: str
    :param whole_head_mask_path: Filepath to mask of the entire head only (no
                                 slice in front of head).
    :type skull_mask_path: str
    :param skull_mask_path: Filepath to the mask of the upper portion of the
                            head only (the whole head mask subtracted by the
                            slice mask).
    :type anatomical_gm_mask: str
    :param anatomical_gm_mask: Filepath to the binary mask of the gray matter.
    :type anatomical_wm_mask: str
    :param anatomical_wm_mask: Filepath to the binary mask of the white
                               matter.
    :type anatomical_csf_mask: str
    :param anatomical_csf_mask: Filepath to the binary mask of the CSF.
    :type subject_id: str
    :param subject_id: The participant ID.
    :type session_id: str
    :param session_id: The session ID.
    :type scan_id: str
    :param scan_id: The scan ID.
    :type run_name: str
    :param run_name: The name of the pipeline run.
    :type site_name: str
    :param site_name: (default: None) The name of the site where the scan was
                      acquired.
    :type exclude_zeroes: bool
    :param exclude_zeroes: (default: False) Whether or not to exclude the
                           pure zero values when defining the background mask.
    :type out_vox: bool
    :param out_vox: (default: True) For FWHM measure: output the FWHM as
                    number of voxels (otherwise as mm).
    :type starter: str
    :param starter: (default: None) If this function is being pulled into a
                    Nipype pipeline, this is the dummy input for the function
                    node.
    :return: A dictionary mapping out the QAP measure values for the current
             participant.
    """

    import os
    from time import strftime
    import qap
    from qap.spatial_qc import summary_mask, snr, cnr, fber, efc, fav, fwhm, \
        cortical_contrast, skew_and_kurt
    from qap.qap_utils import load_image, load_mask, \
                              create_anatomical_background_mask

    # Load the data
    anat_data, anat_aff, anat_hdr = \
        load_image(anatomical_reorient, return_affine=True)

    fg_mask = load_mask(qap_head_mask_path, anatomical_reorient)

    # bg_mask is the inversion of the "qap_head_mask"
    bg_mask = create_anatomical_background_mask(anat_data, fg_mask,
                                                exclude_zeroes)

    whole_head_mask = load_mask(whole_head_mask_path, anatomical_reorient)
    skull_mask = load_mask(skull_mask_path, anatomical_reorient)

    gm_mask = load_mask(anatomical_gm_mask, anatomical_reorient)
    wm_mask = load_mask(anatomical_wm_mask, anatomical_reorient)
    csf_mask = load_mask(anatomical_csf_mask, anatomical_reorient)

    # FBER
    fber_out = fber(anat_data, skull_mask, bg_mask)

    # EFC
    efc_out = efc(anat_data)

    # Artifact
    qi1, _ = fav(fav_artifacts, anatomical_reorient, qap_bg_head_mask_path,
                 qap_head_mask_path)

    # Smoothness in voxels
    tmp = fwhm(anatomical_reorient, whole_head_mask_path, out_vox=out_vox)
    fwhm_x, fwhm_y, fwhm_z, fwhm_out = tmp

    # Summary Measures
    fg, fg_mean, fg_std, fg_size = summary_mask(anat_data, whole_head_mask)
    bg, bg_mean, bg_std, bg_size = summary_mask(anat_data, bg_mask)

    # More Summary Measures
    gm, gm_mean, gm_std, gm_size = summary_mask(anat_data, gm_mask)
    wm, wm_mean, wm_std, wm_size = summary_mask(anat_data, wm_mask)

    # SNR
    snr_out = snr(fg_mean, bg_std)

    # CNR
    cnr_out = cnr(gm_mean, wm_mean, bg_std)

    # Cortical contrast
    cort_out = cortical_contrast(gm_mean, wm_mean)

    # Skewness and Kurtosis
    gm_skew, gm_kurt = skew_and_kurt(gm.flatten())
    wm_skew, wm_kurt = skew_and_kurt(wm.flatten())
    bg_skew, bg_kurt = skew_and_kurt(bg.flatten())

    id_string = "%s %s %s" % (subject_id, session_id, scan_id)
    qap_version = qap.__version__
    qap = {
        id_string: {
            "QAP_Version": "QAP version %s" % qap_version,
            "QAP_pipeline_id": run_name,
            "Time": strftime("%Y-%m-%d %H:%M:%S"),
            "Participant": str(subject_id),
            "Session": str(session_id),
            "Series": str(scan_id),
            "anatomical_spatial":
            {
                "FBER": fber_out,
                "EFC": efc_out,
                "Qi1": qi1,
                "FWHM_x": fwhm_x,
                "FWHM_y": fwhm_y,
                "FWHM_z": fwhm_z,
                "FWHM": fwhm_out,
                "CNR": cnr_out,
                "SNR": snr_out,
                "Cortical Contrast": cort_out,
                "Gray Matter Skewness": gm_skew,
                "Gray Matter Kurtosis": gm_kurt,
                "White Matter Skewness": wm_skew,
                "White Matter Kurtosis": wm_kurt,
                "Background Skewness": bg_skew,
                "Background Kurtosis": bg_kurt
            }
        }
    }

    if site_name:
        qap[id_string]['Site'] = str(site_name)

    if exclude_zeroes:
        qap[id_string]['_zeros_excluded'] = "True"

    for key in qap[id_string]["anatomical_spatial"].keys():
        qap[id_string]["anatomical_spatial"][key] = \
            str(qap[id_string]["anatomical_spatial"][key])

    # prospective filepaths
    if session_output_dir:

        qap[id_string]["filepaths"] = {}

        anat_file = os.path.join(session_output_dir, run_name, site_name,
                                 subject_id, session_id, "anat",
                                 "_".join([subject_id, session_id, scan_id,
                                           "anat-reorient.nii.gz"]))

        if os.path.exists(anat_file):
            qap[id_string]["filepaths"]["scan filepath"] = anat_file

        qi_file = os.path.join(session_output_dir, run_name, site_name,
                               subject_id, session_id, "anat",
                               "_".join([subject_id, session_id, scan_id,
                                         "anat-fav-artifacts-background"
                                         ".nii.gz"]))
        if os.path.exists(qi_file):
            qap[id_string]["filepaths"]["FAV artifacts background"] = qi_file

    return qap


def qap_functional_spatial(mean_epi, func_brain_mask, direction, subject_id,
                           session_id, scan_id, run_name, site_name=None,
                           out_vox=True, starter=None):
    """ Calculate the functional spatial QAP measures for a functional scan.

    - The inclusion of the starter node allows several QAP measure pipelines
      which are not dependent on one another to be executed as one pipeline.
      This allows the MultiProc Nipype plugin to efficiently manage
      resources when parallelizing.

    :type mean_epi: str
    :param mean_epi: Filepath to the mean of the functional timeseries image
                     (should be 3D).
    :type func_brain_mask: str
    :param func_brain_mask: Filepath to the binary mask defining the brain
                            within the functional image.
    :type direction: str
    :param direction: For ghost-to-signal ratio; the phase-encoding direction
                      of the image - this is often "y".
    :type subject_id: str
    :param subject_id: The participant ID.
    :type session_id: str
    :param session_id: The session ID.
    :type scan_id: str
    :param scan_id: The scan ID.
    :type run_name: str
    :param run_name: The pipeline run name.
    :type site_name: str
    :param site_name: (default: None) The name of the site where the scan was
                      acquired.
    :type out_vox: bool
    :param out_vox: (default: True) For FWHM measure: output the FWHM as
                    number of voxels (otherwise as mm).
    :type starter: str
    :param starter: (default: None) If this function is being pulled into a
                    Nipype pipeline, this is the dummy input for the function
                    node.
    :rtype: dict
    :return: A dictionary mapping out the QAP measure values for the current
             participant.
    """

    from time import strftime

    import qap
    from qap.spatial_qc import summary_mask, snr, fber, efc, fwhm, \
        ghost_direction
    from qap.qap_utils import load_image, load_mask

    # Load the data
    anat_data = load_image(mean_epi)
    fg_mask = load_mask(func_brain_mask, mean_epi)
    bg_mask = 1 - fg_mask

    # FBER
    fber_out = fber(anat_data, fg_mask, bg_mask)

    # EFC
    efc_out = efc(anat_data)
    
    # Smoothness in voxels
    tmp = fwhm(mean_epi, func_brain_mask, out_vox=out_vox)
    fwhm_x, fwhm_y, fwhm_z, fwhm_out = tmp

    # Summary Measures
    fg, fg_mean, fg_std, fg_size = summary_mask(anat_data, fg_mask)
    bg, bg_mean, bg_std, bg_size = summary_mask(anat_data, bg_mask)

    # SNR
    snr_out = snr(fg_mean, bg_std)

    id_string = "%s %s %s" % (subject_id, session_id, scan_id)
    qap = {
        id_string:
        {
           "QAP_Version": "QAP version %s" % qap.__version__,
           "QAP_pipeline_id": run_name,
           "Time": strftime("%Y-%m-%d %H:%M:%S"),
           "Participant": str(subject_id),
           "Session": str(session_id),
           "Series": str(scan_id),
           "functional_spatial":
           {
              "FBER": fber_out,
              "EFC": efc_out,
              "FWHM": fwhm_out,
              "FWHM_x": fwhm_x,
              "FWHM_y": fwhm_y,
              "FWHM_z": fwhm_z,
              "SNR": snr_out
           }
        }
    }

    # Ghosting
    if direction == "all":
        qap[id_string]["functional_spatial"]['Ghost_x'] = \
            ghost_direction(anat_data, fg_mask, "x")
        qap[id_string]["functional_spatial"]['Ghost_y'] = \
            ghost_direction(anat_data, fg_mask, "y")
        qap[id_string]["functional_spatial"]['Ghost_z'] = \
            ghost_direction(anat_data, fg_mask, "z")
    else:
        qap[id_string]["functional_spatial"]['Ghost_%s' % direction] = \
            ghost_direction(anat_data, fg_mask, direction)

    if site_name:
        qap[id_string]['Site'] = str(site_name)

    for key in qap[id_string]["functional_spatial"].keys():
        qap[id_string]["functional_spatial"][key] = \
            str(qap[id_string]["functional_spatial"][key])

    return qap


def qap_functional_temporal(
        func_timeseries, func_mean, func_brain_mask, bg_func_brain_mask,
        fd_file, sfs, subject_id, session_id, scan_id, run_name,
        site_name=None, session_output_dir=None, starter=None):
    """ Calculate the functional temporal QAP measures for a functional scan.

    - The inclusion of the starter node allows several QAP measure pipelines
      which are not dependent on one another to be executed as one pipeline.
      This allows the MultiProc Nipype plugin to efficiently manage resources
      when parallelizing.

    :type func_timeseries: str
    :param func_timeseries: Filepath to the 4D functional timeseries.
    :type func_mean: str
    :param func_mean: Filepath to the averaged functional timeseries.
    :type func_brain_mask: str
    :param func_brain_mask: Filepath to the binary mask defining the brain
                            within the functional image.
    :type bg_func_brain_mask: str
    :param bg_func_brain_mask: Filepath to the inversion of the functional
                               brain mask.
    :type fd_file: str
    :param fd_file: File containing the RMSD values (calculated previously).
    :type subject_id: str
    :param subject_id: The participant ID.
    :type session_id: str
    :param session_id: The session ID.
    :type scan_id: str
    :param scan_id: The scan ID.
    :type run_name: str
    :param run_name: The pipeline run name.
    :type site_name: str
    :param site_name: (default: None) The name of the site where the scan was
                      acquired.
    :type starter: str
    :param starter: (default: None) If this function is being pulled into a
                    Nipype pipeline, this is the dummy input for the function
                    node.
    :rtype: dict
    :return: A dictionary mapping out the QAP measure values for the current
             participant.
    :rtype: dict
    :return: A dictionary containing the output measure timeseries for the
             quality assurance output file for visualization and reporting.
    """

    import os
    import numpy as np
    import nibabel as nb
    from time import strftime

    import qap
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
    outliers = outlier_timepoints(func_timeseries, mask_file=func_brain_mask)
    # calculate the outliers of the outliers! AAHH!
    outlier_perc_out, outlier_IQR = calculate_percent_outliers(outliers)

    # 3dTout (outside of brain)
    oob_outliers = outlier_timepoints(func_timeseries,
                                      mask_file=bg_func_brain_mask)
    oob_outlier_perc_out, oob_outlier_IQR = \
        calculate_percent_outliers(oob_outliers)

    # 3dTqual
    quality = quality_timepoints(func_timeseries)
    quality_outliers, quality_IQR = calculate_percent_outliers(quality)

    # Signal Fluctuation Sensitivity (SFS)
    sfs_data = nb.load(sfs).get_data()
    sfs_inbrain = sfs_data.nonzero()

    # GCOR
    gcor = global_correlation(func_timeseries, func_brain_mask)

    # Compile
    id_string = "%s %s %s" % (subject_id, session_id, scan_id)
    qap_version = qap.__version__
    qap = {
        id_string:
        {
          "QAP_Version": "QAP version %s" % qap_version,
          "QAP_pipeline_id": run_name,
          "Time": strftime("%Y-%m-%d %H:%M:%S"),
          "Participant": str(subject_id),
          "Session": str(session_id),
          "Series": str(scan_id),
          "functional_temporal":
          {
             "Std DVARS (Mean)": mean_dvars,
             "Std DVARS (Std Dev)": np.std(dvars),
             "Std DVARS (Median)": np.median(dvars),
             "Std DVARs IQR": dvars_IQR,
             "Std DVARS percent outliers": dvars_outliers,
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
             "Fraction of OOB Outliers (Mean)": np.mean(oob_outliers),
             "Fraction of OOB Outliers (Std Dev)": np.std(oob_outliers),
             "Fraction of OOB Outliers (Median)": np.median(oob_outliers),
             "Fraction of OOB Outliers IQR": oob_outlier_IQR,
             "Fraction of OOB Outliers percent outliers": oob_outlier_perc_out,
             "Quality (Mean)": np.mean(quality),
             "Quality (Std Dev)": np.std(quality),
             "Quality (Median)": np.median(quality),
             "Quality IQR": quality_IQR,
             "Quality percent outliers": quality_outliers,
             "Signal Fluctuation Sensitivity (Mean)": np.mean(sfs_inbrain),
             "GCOR": gcor
          }
        }
    }

    if site_name:
        qap[id_string]['Site'] = str(site_name)

    for key in qap[id_string]["functional_temporal"].keys():
        qap[id_string]["functional_temporal"][key] = \
            str(qap[id_string]["functional_temporal"][key])

    # prospective filepaths
    if session_output_dir:

        qap[id_string]["filepaths"] = {}

        func_file = os.path.join(session_output_dir, run_name, site_name,
                                 subject_id, session_id, "func",
                                 "_".join([subject_id, session_id, scan_id,
                                           "func-mean.nii.gz"]))
        if os.path.exists(func_file):
            qap[id_string]["filepaths"]["scan filepath"] = func_file

        tstd_file = os.path.join(session_output_dir, run_name, site_name,
                                 subject_id, session_id, "func",
                                 "_".join([subject_id, session_id, scan_id,
                                           "func-temporal-std-map.nii.gz"]))
        if os.path.exists(tstd_file):
            qap[id_string]["filepaths"]["temporal STD"] = tstd_file

        estn_file = os.path.join(session_output_dir, run_name, site_name,
                                 subject_id, session_id, "func",
                                 "_".join([subject_id, session_id, scan_id,
                                           "func-estimated-nuisance.nii.gz"]))
        if os.path.exists(estn_file):
            qap[id_string]["filepaths"]["estimated nuisance"] = estn_file

        sfs_file = os.path.join(session_output_dir, run_name, site_name,
                                subject_id, session_id, "func",
                                "_".join([subject_id, session_id, scan_id,
                                          "func-SFS.nii.gz"]))
        if os.path.exists(sfs_file):
            qap[id_string]["filepaths"]["SFS filepath"] = sfs_file

    return qap