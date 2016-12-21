

def create_expr_string(clip_level_value):
    """Create the expression arg string to run AFNI 3dcalc via Nipype.

    Keyword Arguments:
      clip_level_value -- [integer] the integer of the clipping threshold

    Returns:
      expr_string -- [string] the string intended for the Nipype AFNI 3dcalc 
                     "expr" arg inputs
    """

    expr_string = "step(a-%s)" % clip_level_value

    return expr_string


def read_nifti_image(nifti_infile):
    """Read a NIFTI file into Nibabel-format image data.

    Keyword arguments:
      nifti_infile -- [string] the filepath of the NIFTI image to read in

    Returns:
      nifti_img -- [Nibabel image] image data in NIFTI format
    """

    import nibabel as nb
    from qap.workflow_utils import raise_smart_exception

    try:
        nifti_img = nb.load(nifti_infile)
    except:
        err = "\n\n[!] Could not load the NIFTI image using Nibabel:\n" \
              "%s\n\n" % nifti_infile
        raise_smart_exception(locals(), err)

    return nifti_img


def write_nifti_image(nifti_img, file_path):
    """Write image data in Nibabel format into a NIFTI file.

    Keyword arguments:
      nifti_img -- [Nibabel image] the image data Nibabel object to write out
      file_path -- [string] the filepath of the NIFTI image to create

    Returns:
      N/A
    """

    import nibabel as nb
    from qap.workflow_utils import raise_smart_exception

    try:
        nb.save(nifti_img, file_path)
    except:
        err = "\n\n[!] Could not save the NIFTI image using Nibabel:\n" \
              "%s\n\n" % file_path
        raise_smart_exception(locals(), err)


def read_json(json_filename):
    """Read the contents of a JSON file.

    Keyword arguments:
      json_filename -- [string] the path to the JSON file

    Returns:
      json_dict -- [Python dictionary] dictionary containing the info from the
                   JSON file
    """

    import os
    import json
    from qap.workflow_utils import raise_smart_exception

    if not os.path.exists(json_filename):
        err = "\n\n[!] The JSON file provided does not exist.\nFilepath: " \
              "%s\n\n" % json_filename
        raise_smart_exception(locals(),err)

    with open(json_filename, "r") as f:
        json_dict = json.load(f)

    return json_dict


def write_json(output_dict, json_file):
    """Either update or write a dictionary to a JSON file.

    Keyword arguments:
      output_dict -- [Python dictionary] the dictionary to write or append to 
                     the JSON file
      json_file -- [string] the filepath of the JSON file to write or update

    Returns:
      json_file -- [string] filepath of the JSON file written to
    """

    import os
    import json
    from lockfile import FileLock

    from qap.qap_workflows_utils import read_json

    write = True

    if os.path.exists(json_file):
        current_dict = read_json(json_file)
        if current_dict == output_dict:
            # nothing to update
            write = False
        else:
            for key in output_dict.keys():
                try:
                    current_dict[key].update(output_dict[key])
                except KeyError:
                    current_dict[key] = output_dict[key]
    else:
        current_dict = output_dict

    if write:
        lock = FileLock(json_file)
        lock.acquire()
        with open(json_file, "wt") as f:
            json.dump(current_dict, f, indent=2, sort_keys=True)
        lock.release()

    if os.path.exists(json_file):
        return json_file


def convert_allineate_xfm(mat_list):
    """Convert the affine transform output of AFNI's 3dAllineate into an
    equivalent 4x4 matrix.

    Keyword arguments:
      mat_list -- [Python list] a vector (list) of the flattened affine matrix

    Returns:
      allineate_mat -- [Numpy array] a NumPy array of the converted affine 
                       matrix

    Notes:
      - Takes the 3x3 + offset format of the input matrix and turns it into a
        4x4 with the last row being 0,0,0,1
    """

    import numpy as np

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

    return allineate_mat


def warp_coordinates(inpoint, allineate_mat, infile_affine, infile_dims):
    """Warp spatial coordinates using a 4x4 affine matrix.

    Keyword arguments:
      inpoint -- [Python list] a list of three numbers describing a coordinate 
                 in 3D space
      allineate_mat -- [Numpy array] a NumPy array describing the 4x4 affine 
                       transform
      infile_affine -- [Nibabel affine] a NiBabel NIFTI affine info object
      infile_dims -- [Python tuple] a tuple containing the dimensions of the 
                     NIFTI file the coordinates pertain to

    Returns:
      co_nums_newlist -- [Python list] a list of three warped coordinates
    """

    import numpy as np
    import numpy.linalg as npl
    import nibabel as nb

    # using the transform, calculate what the three coordinates are in
    # native space, as it corresponds to the anatomical scan
    coord_out = \
        list(np.dot(np.linalg.inv(allineate_mat),inpoint))

    # remove the one resulting zero at the end
    coord_out = coord_out[:-1]

    # convert the coordinates from mm to voxels
    coord_out = \
        nb.affines.apply_affine(npl.inv(infile_affine),coord_out)

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

    Keyword arguments:
      coords_list -- [Python list] a list of lists, describing the three 
                     coordinate points of the triangular plane
      infile_dims -- [Python list] a list of the NIFTI file's dimensions

    Returns:
      plane_dict -- [Python dictionary] a dictionary mapping the z coordinates
                    to the (x,y) coordinate pairs
    """

    import numpy as np

    # get the vectors connecting the points
    u = []

    for a_pt, c_pt in zip(coords_list[0], coords_list[2]):
        u.append(int(a_pt - c_pt))

    v = []

    for b_pt, c_pt in zip(coords_list[1], coords_list[2]):
        v.append(int(b_pt - c_pt))

    u_vector = np.asarray(u)
    v_vector = np.asarray(v)

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


def create_slice_mask(plane_dict, infile_dims):
    """Create a binary array defining a mask covering the area below a given
    plane in the 3D image.

    Keyword arguments:
      plane_dict -- [Python dictionary] a dictionary matching z voxel 
                    coordinates to corresponding (x, y) coordinates
      infile_dims -- [Python list] a list of the NIFTI file's dimensions

    Returns:
      mask_array -- [Numpy array] a NumPy array defining the binary mask of 
                    the slice
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

    Keyword arguments:
      infile -- [string] filepath to the participant's anatomical scan
      transform -- [string] filepath to the text file containing the affine 
                   matrix output of AFNI's 3dAllineate describing the warp 
                   from the anatomical scan to a template

    Returns:
      outfile_path -- [string] filepath to the new head mask NIFTI file
    """

    import os
    import sys

    import nibabel as nb
    import numpy as np
    import numpy.linalg as npl
    import subprocess
    import pkg_resources as p

    from qap.script_utils import read_txt_file
    from qap.qap_workflows_utils import read_nifti_image, \
                                        convert_allineate_xfm, \
                                        warp_coordinates, \
                                        calculate_plane_coords, \
                                        create_slice_mask, \
                                        write_nifti_image

    # get file info
    infile_img = read_nifti_image(infile)

    infile_header = infile_img.get_header()
    infile_affine = infile_img.get_affine()
    infile_dims = infile_header.get_data_shape()

    # get the affine output matrix of 3dallineate
    allineate_mat_list = read_txt_file(transform)

    # get the 3dAllineate output affine matrix into a list
    mat_list = filter(None,allineate_mat_list[1].rstrip("\n").split(" "))

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


def create_header_dict_entry(in_file, subject, session, scan, type):
    """Gather the header information from a NIFTI file and arrange it into a
    Python dictionary.

    Keyword arguments:
      in_file -- [string] filepath to the NIFTI raw data scan
      subject -- [string] the participant ID
      session -- [string] the session ID
      scan -- [string] the scan ID
      type -- [string] the data type ("anatomical" or "functional")

    Returns:
      qap_dict -- [Python dictionary] a dictionary with the header information
                  of the file assigned to the participant's data
    """

    import nibabel as nb
    from qap.workflow_utils import raise_smart_exception

    try:
        img = nb.load(in_file)
        img_header = img.header
    except:
        err = "You may not have an up-to-date installation of the Python " \
              "Nibabel package.\nYour Nibabel version: %s" % \
              str(nb.__version__)
        raise_smart_exception(locals(),err)

    subkey = "%s_header_info" % type
    id_string = "%s %s %s" % (subject, session, scan)
    qap_dict = {id_string: {subkey: {}}}

    info_labels = ["descrip", "db_name", "bitpix", "slice_start", \
                   "scl_slope", "scl_inter", "slice_end", "slice_duration", \
                   "toffset", "quatern_b", "quatern_c", "quatern_d", \
                   "qoffset_x", "qoffset_y", "qoffset_z", "srow_x", "srow_y",\
                   "srow_z", "aux_file", "intent_name", "slice_code", \
                   "data_type", "qform_code", "sform_code"]

    for info_label in info_labels:
        try:
            qap_dict[id_string][subkey][info_label] = str(img_header[info_label])
        except:
            print "\n\n%s field not in NIFTI header of %s\n\n" % \
                  (info_label, in_file)
            qap_dict[id_string][subkey][info_label] = ""
            pass

    try:
        pixdim = img_header['pixdim']
        qap_dict[id_string][subkey]["pix_dimx"] = str(pixdim[1])
        qap_dict[id_string][subkey]["pix_dimy"] = str(pixdim[2])
        qap_dict[id_string][subkey]["pix_dimz"] = str(pixdim[3])
        qap_dict[id_string][subkey]["tr"] = str(pixdim[4])
    except:
        print "\n\npix_dim/TR fields not in NIFTI header of %s\n\n" % in_file
        pass

    try:
        qap_dict[id_string][subkey]["extensions"] = \
            len(img.header.extensions.get_codes())
    except:
        print "\n\nExtensions not in NIFTI header of %s\n\n" % in_file
        pass

    return qap_dict


def qap_anatomical_spatial(anatomical_reorient, qap_head_mask_path,
                           whole_head_mask_path, skull_mask_path,
                           anatomical_gm_mask, anatomical_wm_mask,
                           anatomical_csf_mask, subject_id, session_id,
                           scan_id, site_name=None, exclude_zeroes=False,
                           out_vox=True, starter=None):
    """ Calculate the anatomical spatial QAP measures for an anatomical scan.

    Keyword arguments:
      anatomical_reorient -- [string] filepath to the reoriented anatomical 
                             scan
      qap_head_mask_path -- [string] filepath to mask of the head, plus the 
                            slice covering the region below the nose and in 
                            front of the mouth
      whole_head_mask_path -- [string] filepath to mask of the entire head 
                              only (no slice in front of the mouth)
      skull_mask_path -- [string] filepath to the mask of the upper portion of
                         the head only (the whole head mask subtracted by the 
                         slice mask)
      anatomical_gm_mask -- [string] filepath to the binary mask of the gray 
                            matter
      anatomical_wm_mask -- [string] filepath to the binary mask of the white 
                            matter
      anatomical_csf_mask -- [string] filepath to the binary mask of the CSF
      subject_id -- [string] the participant ID
      session_id -- [string] the session ID
      scan_id -- [string] the scan ID
      site_name -- [string] (default: None) the name of the site where the 
                   scan was acquired
      exclude_zeroes -- [boolean] (default: False) whether or not to exclude 
                        the pure zero values when defining the background mask
      out_vox -- [boolean] (default: True) for FWHM measure: output the FWHM 
                 as # of voxels (otherwise as mm)
      starter -- [string] (default: None) if this function is being pulled 
                 into a Nipype pipeline, this is the dummy input for the 
                 function node

    Returns:
      qc -- [Python dictionary] a dictionary mapping out the QAP measure 
            values for the current participant

    Notes:
      - The exclude_zeroes flag is useful for when a large amount of zero
        values have been artificially injected into the image, for example,
        when removing the faces and ears in scans for privacy compliance
        reasons; these zeroes can artificially skew the spatial quality metric
        results and make it seem that there is far less noise or artifacts in
        the image than there really is
      - The inclusion of the starter node allows several QAP measure pipelines
        which are not dependent on one another to be executed as one pipeline.
        This allows the MultiProc Nipype plugin to efficiently manage
        resources when parallelizing.
    """

    import os
    import sys
    from time import strftime

    import qap
    from qap.spatial_qc import summary_mask, snr, cnr, fber, efc, \
        artifacts, fwhm, cortical_contrast
    from qap.qap_utils import load_image, load_mask, \
                              create_anatomical_background_mask

    # Load the data
    anat_data = load_image(anatomical_reorient)

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
    qi1, _ = artifacts(anat_data, fg_mask, bg_mask, calculate_qi2=False)

    # Smoothness in voxels
    tmp = fwhm(anatomical_reorient, whole_head_mask_path, out_vox=out_vox)
    fwhm_x, fwhm_y, fwhm_z, fwhm_out = tmp

    # Summary Measures
    fg_mean, fg_std, fg_size = summary_mask(anat_data, whole_head_mask)
    bg_mean, bg_std, bg_size = summary_mask(anat_data, bg_mask)

    gm_mean, gm_std, gm_size = (None, None, None)
    wm_mean, wm_std, wm_size = (None, None, None)
    csf_mean, csf_std, csf_size = (None, None, None)

    # More Summary Measures
    gm_mean, gm_std, gm_size = summary_mask(anat_data, gm_mask)
    wm_mean, wm_std, wm_size = summary_mask(anat_data, wm_mask)
    csf_mean, csf_std, csf_size = summary_mask(anat_data, csf_mask)

    # SNR
    snr_out = snr(fg_mean, bg_std)

    # CNR
    cnr_out = cnr(gm_mean, wm_mean, bg_std)

    # Cortical contrast
    cort_out = cortical_contrast(gm_mean, wm_mean)

    id_string = "%s %s %s" % (subject_id, session_id, scan_id)
    qc = {
            id_string:
            {
               "QAP_pipeline_id": "QAP version %s" % qap.__version__,
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
                  "Cortical Contrast": cort_out
               }
            }
    }

    if site_name:
        qc[id_string]['Site'] = str(site_name)

    if exclude_zeroes:
        qc[id_string]['_zeros_excluded'] = "True"

    for key in qc[id_string]["anatomical_spatial"].keys():
        qc[id_string]["anatomical_spatial"][key] = \
            str(qc[id_string]["anatomical_spatial"][key])

    return qc


def qap_functional_spatial(mean_epi, func_brain_mask, direction, subject_id,
                           session_id, scan_id, site_name=None, out_vox=True,
                           starter=None):
    """ Calculate the functional spatial QAP measures for a functional scan.

    Keyword arguments:
      mean_epi -- [string] filepath to the mean of the functional timeseries 
                  image (should be 3D)
      func_brain_mask -- [string] filepath to the binary mask defining the 
                         brain within the functional image
      direction -- [string] for ghost-to-signal ratio; the phase encoding 
                   direction of the image - this is often "y"
      subject_id -- [string] the participant ID
      session_id -- [string] the session ID
      scan_id -- [string] the scan ID
      site_name -- [string] (default: None) the name of the site where the 
                   scan was acquired
      out_vox -- [boolean] (default: True) for FWHM measure: output the FWHM 
                 as # of voxels (otherwise as mm)
      starter -- [string] (default: None) if this function is being pulled 
                 into a Nipype pipeline, this is the dummy input for the
                 function node

    Returns:
      qc -- [Python dictionary] a dictionary mapping out the QAP measure 
            values for the current participant

    Notes:
      - The inclusion of the starter node allows several QAP measure pipelines
        which are not dependent on one another to be executed as one pipeline.
        This allows the MultiProc Nipype plugin to efficiently manage
        resources when parallelizing.
    """

    import os
    import sys
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
    fg_mean, fg_std, fg_size = summary_mask(anat_data, fg_mask)
    bg_mean, bg_std, bg_size = summary_mask(anat_data, bg_mask)

    # SNR
    snr_out = snr(fg_mean, bg_std)

    id_string = "%s %s %s" % (subject_id, session_id, scan_id)
    qc = {
            id_string:
            {
               "QAP_pipeline_id": "QAP version %s" % qap.__version__,
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
    if (direction == "all"):
        qc[id_string]["functional_spatial"]['Ghost_x'] = \
            ghost_direction(anat_data, fg_mask, "x")
        qc[id_string]["functional_spatial"]['Ghost_y'] = \
            ghost_direction(anat_data, fg_mask, "y")
        qc[id_string]["functional_spatial"]['Ghost_z'] = \
            ghost_direction(anat_data, fg_mask, "z")
    else:
        qc[id_string]["functional_spatial"]['Ghost_%s' % direction] = \
            ghost_direction(anat_data, fg_mask, direction)

    if site_name:
        qc[id_string]['Site'] = str(site_name)

    for key in qc[id_string]["functional_spatial"].keys():
        qc[id_string]["functional_spatial"][key] = \
            str(qc[id_string]["functional_spatial"][key])

    return qc


def qap_functional_temporal(
        func_timeseries, func_brain_mask, bg_func_brain_mask, fd_file,
        subject_id, session_id, scan_id, site_name=None, starter=None):
    """ Calculate the functional temporal QAP measures for a functional scan.

    Keyword arguments:
      func_timeseries -- [string] filepath to the 4D functional timeseries
      func_brain_mask -- [string] filepath to the binary mask defining the 
                         brain within the functional image
      bg_func_brain_mask -- [string] filepath to the inversion of the 
                            functional brain mask
      fd_file -- [string] file containing the RMSD values (calculated 
                 previously)
      subject_id -- [string] the participant ID
      session_id -- [string] the session ID
      scan_id -- [string] the scan ID
      site_name -- [string] (default: None) the name of the site where the 
                   scan was acquired
      starter -- [string] (default: None) if this function is being pulled 
                 into a Nipype pipeline, this is the dummy input for the 
                 function node

    Returns:
      qc -- [Python dictionary] a dictionary mapping out the QAP measure 
            values for the current participant

    Notes:
      - The inclusion of the starter node allows several QAP measure pipelines
        which are not dependent on one another to be executed as one pipeline.
        This allows the MultiProc Nipype plugin to efficiently manage
        resources when parallelizing.
    """

    import sys
    import nibabel as nb
    import numpy as np
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

    # GCOR
    gcor = global_correlation(func_timeseries, func_brain_mask)

    # Compile
    id_string = "%s %s %s" % (subject_id, session_id, scan_id)
    qc = {
            id_string:
            {
              "QAP_pipeline_id": "QAP version %s" % qap.__version__,
              "Time": strftime("%Y-%m-%d %H:%M:%S"),
              "Participant": str(subject_id),
              "Session": str(session_id),
              "Series": str(scan_id),
              "functional_temporal":
              {
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
                 "GCOR": gcor
              }
            }
    }

    if site_name:
        qc[id_string]['Site'] = str(site_name)

    for key in qc[id_string]["functional_temporal"].keys():
        qc[id_string]["functional_temporal"][key] = \
            str(qc[id_string]["functional_temporal"][key])

    return qc
