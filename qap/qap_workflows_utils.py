

def select_thresh(input_skull):

    import os
    import commands

    avg_in = "3dmaskave %s" % input_skull

    avg_out = commands.getoutput(avg_in)


    avg = int(float(avg_out.split("\n")[-1].split(" ")[0]))
    max_limit = int(avg * 3)


    # get the voxel intensity bins
    cmd_in = "3dHist -input %s -nbin 10 -max %d -showhist" % \
             (input_skull, max_limit)

    cmd_out = commands.getoutput(cmd_in)

    os.system("rm HistOut.niml.hist")

    bins = {}

    for line in cmd_out.split("\n"):

        if "*" in line and not line.startswith("*"):

            vox_bin = line.replace(" ","").split(":")[0]

            voxel_value = int(float(vox_bin.split(",")[0]))

            bins[int(vox_bin.split(",")[1])] = voxel_value


    thresh_out = bins[min(bins.keys())]

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
    inpoint_c = "0 88 -72" # nose, apparently


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

        for ind in range(0,3):

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
    n = np.cross(u,v)

    # normalize the vector
    n = n / np.linalg.norm(n,2)
      

    constant = np.dot(n, np.asarray(new_coords[0]))


    # now determine the z-coordinate for each pair of x,y
    plane_dict = {}

    for yvox in range(0,infile_dims[1]):

        for xvox in range(0, infile_dims[0]):

            zvox = (constant - (n[0]*xvox + n[1]*yvox))/n[2]

            zvox = np.floor(zvox)

            if zvox < 1:
                zvox = 1
            elif zvox > infile_dims[2]:
                zvox = infile_dims[2]

            plane_dict[(xvox,yvox)] = zvox



    # create the mask
    mask_array = np.zeros(infile_dims)

    for x in range(0, infile_dims[0]):

        for y in range(0, infile_dims[1]):

            for z in range(0, infile_dims[2]):

                if plane_dict[(x,y)] > z:
                    mask_array[x,y,z] = 1


    new_mask_img = nb.Nifti1Image(mask_array, infile_affine, infile_header)

    infile_filename = infile.split("/")[-1].split(".")[0]

    outfile_name = infile_filename + "_slice_mask.nii.gz"
    outfile_path = os.path.join(os.getcwd(), outfile_name)

    nb.save(new_mask_img, outfile_path)


    return outfile_path



def qap_anatomical_spatial(anatomical_reorient, head_mask_path, \
                               anatomical_gm_mask, anatomical_wm_mask, \
                               anatomical_csf_mask, subject_id, session_id, \
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
    qc              = dict()

    qc['subject'] = subject_id

    qc['session'] = session_id

    qc['scan'] = scan_id

    if site_name:
        qc['site'] = site_name

   
    # FBER
    qc['fber'] = fber(anat_data, fg_mask)
    
    # EFC
    qc['efc'] = efc(anat_data)
    
    # Artifact
    qc['qi1'], _    = artifacts(anat_data, fg_mask, calculate_qi2=False)
    
    # Smoothness in voxels
    tmp             = fwhm(anatomical_reorient, head_mask_path, out_vox=out_vox)
    qc['fwhm_x'], qc['fwhm_y'], qc['fwhm_z'], qc['fwhm'] = tmp
    
    
    # Summary Measures
    qc['fg_mean'], qc['fg_std'], qc['fg_size']      = summary_mask(anat_data, fg_mask)
    qc['bg_mean'], qc['bg_std'], qc['bg_size']      = summary_mask(anat_data, bg_mask)
    
    qc['gm_mean'], qc['gm_std'], qc['gm_size']      = (None, None, None)
    qc['wm_mean'], qc['wm_std'], qc['wm_size']      = (None, None, None)
    qc['csf_mean'], qc['csf_std'], qc['csf_size']   = (None, None, None)
    qc['cnr']   = None
    qc['snr']   = None
    
    # More Summary Measures
    qc['gm_mean'], qc['gm_std'], qc['gm_size']      = summary_mask(anat_data, gm_mask)
    qc['wm_mean'], qc['wm_std'], qc['wm_size']      = summary_mask(anat_data, wm_mask)
    qc['csf_mean'], qc['csf_std'], qc['csf_size']   = summary_mask(anat_data, csf_mask)

    # SNR
    qc['snr']       = snr(qc['fg_mean'], qc['bg_std'])

    # CNR
    qc['cnr']   = cnr(qc['gm_mean'], qc['wm_mean'], qc['bg_std'])
    

    return qc



def qap_functional_spatial(mean_epi, func_brain_mask, direction, subject_id, \
                               session_id, scan_id, site_name=None, \
                               out_vox=True):

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
    qc              = dict()

    qc['subject'] = subject_id

    qc['session'] = session_id

    qc['scan'] = scan_id

    if site_name:
        qc['site'] = site_name

   
    # FBER
    qc['fber'] = fber(anat_data, fg_mask)
    
    # EFC
    qc['efc'] = efc(anat_data)
    
    
    # Smoothness in voxels
    tmp             = fwhm(mean_epi, func_brain_mask, out_vox=out_vox)
    qc['fwhm_x'], qc['fwhm_y'], qc['fwhm_z'], qc['fwhm'] = tmp
    
    # Ghosting
    if (direction == "all"):
        qc['ghost_x'] = ghost_direction(anat_data, fg_mask, "x")
        qc['ghost_y'] = ghost_direction(anat_data, fg_mask, "y")
        qc['ghost_z'] = ghost_direction(anat_data, fg_mask, "z")

    else:
        qc['ghost_%s' % direction] = ghost_direction(anat_data, fg_mask, \
                                         direction)

    
    # Summary Measures
    qc['fg_mean'], qc['fg_std'], qc['fg_size'] = summary_mask(anat_data, fg_mask)
    qc['bg_mean'], qc['bg_std'], qc['bg_size'] = summary_mask(anat_data, bg_mask)
    

    qc['snr']   = None
    

    # SNR
    qc['snr']       = snr(qc['fg_mean'], qc['bg_std'])
    

    return qc



def qap_functional_temporal(func_motion_correct, func_brain_mask, \
                                coord_xfm_matrix, subject_id, session_id, \
                                scan_id, site_name=None, \
                                motion_threshold=1.0):

    import sys

    from qap.temporal_qc import mean_dvars_wrapper, summarize_fd, \
                                mean_outlier_timepoints, \
                                mean_quality_timepoints, global_correlation

    # DVARS
    mean_dvars  = mean_dvars_wrapper(func_motion_correct, func_brain_mask)

    # Mean FD (Jenkinson)
    (mean_fd, num_fd, percent_fd) = summarize_fd(coord_xfm_matrix, \
                                                 threshold=motion_threshold)

    # 3dTout
    mean_outlier = mean_outlier_timepoints(func_motion_correct, \
                                               func_brain_mask)

    # 3dTqual
    mean_quality = mean_quality_timepoints(func_motion_correct)

    # new thing
    gcor = global_correlation(func_motion_correct, func_brain_mask)

    # Compile
    qc = {
        "subject":  subject_id,
        "session":  session_id,
        "scan":     scan_id,
        "dvars":    mean_dvars, 
        "mean_fd":  mean_fd, 
        'num_fd':   num_fd, 
        'perc_fd':  percent_fd, 
        "outlier":  mean_outlier,
        "quality":  mean_quality,
        "gcor": gcor
    }

    if site_name:
        qc['site'] = site_name
    

    return qc
    
    
    
def write_to_csv(sub_qap_dict):  #, outfile):

    import os
    import csv

    fields = sub_qap_dict.keys()

    # put these at the forefront of the list of header items, to make the
    # output CSV's more readable

    fields = sorted(fields)

    if "subject" in fields:
        fields.remove("subject")
        fields.insert(0, "subject")

    if "session" in fields:
        fields.remove("session")
        fields.insert(1, "session")

    if "scan" in fields:
        fields.remove("scan")
        fields.insert(2, "scan")

    if "site" in fields:
        fields.remove("site")
        fields.insert(3, "site")


    outfile = os.path.join(os.getcwd(), "qap_measures.csv")


    with open(outfile, "wt") as out_f:

        csv_writer = csv.DictWriter(out_f, fields)

        csv_writer.writeheader()

        csv_writer.writerow(sub_qap_dict)


    return outfile



