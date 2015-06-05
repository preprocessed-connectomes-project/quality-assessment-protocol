
import os
import sys


base_QAP_dir = "/tdata/qap_project/QAP"

base_test_dir = os.path.join(base_QAP_dir, "qc_test")

test_sub_dir = os.path.join(base_test_dir, "pipeline_folder", "2014113")



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

        if "*" in line:

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
    import commands
    
    import pkg_resources as p

    # get the directory this script is in (not the current working one)
    script_dir = os.path.dirname(os.path.realpath('__file__'))

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
    inpoint_files = [p.resource_filename("pcp_qap", "qclib/inpoint_a.txt"),
                     p.resource_filename("pcp_qap", "qclib/inpoint_b.txt"),
                     p.resource_filename("pcp_qap", "qclib/inpoint_c.txt")]


    # convert the ANTS anat->standard affine to FSL format, so we can use
    # std2imgcoord below

    try:

        c3d_cmd = "c3d_affine_tool -itk %s -ref %s -src %s -ras2fsl -o " \
                  "qc_fsl_affine_xfm.mat" % (transform, standard, infile)

        #print "\nExecuting:\n%s\n" % c3d_cmd

        os.system(c3d_cmd)

    except:

        err = "\n[!] The C3D Affine Tool failed to convert the ANTS " \
              "affine transform file to FSL format.\n\nCommand: %s" \
              % c3d_cmd
        raise Exception(err)



    # let's convert the coordinates into voxel coordinates

    coords = []

    for inpoint in inpoint_files:

        coord_cmd = "std2imgcoord -img %s -std %s -xfm " \
                    "qc_fsl_affine_xfm.mat -vox %s" \
                    % (infile, standard, inpoint)

        try:

            coord_out = commands.getoutput(coord_cmd)

        except:

            raise Exception


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
    
    

def test_select_thresh():

    input_skull = os.path.join(test_sub_dir, "anatomical_reorient",
                               "mprage_resample.nii.gz")

    thresh_out = select_thresh(input_skull)


    assert thresh_out == 147



def test_slice_head_mask():

    import nibabel as nb
    import numpy as np

    infile = os.path.join(test_sub_dir, "anatomical_reorient",
                          "mprage_resample.nii.gz")

    transform = os.path.join(test_sub_dir, "ants_affine_xfm",
                             "transform2Affine.mat")

    standard = os.path.join(base_test_dir, "MNI152_T1_2mm.nii.gz")


    slice_mask_path = slice_head_mask(infile, transform, standard)

    # get file info
    slice_mask_img = nb.load(slice_mask_path)

    slice_mask_data = slice_mask_img.get_data()

    # test file
    test_mask_img = nb.load(os.path.join(base_test_dir,
                            "test_slice_mask.nii.gz"))

    test_mask_data = test_mask_img.get_data()
    

    os.system("rm %s" % slice_mask_path)


    assert slice_mask_data.all() == test_mask_data.all()
  


def run_all_tests():

    test_select_thresh()
    test_slice_head_mask()


