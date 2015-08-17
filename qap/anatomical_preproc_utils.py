
def ants_lin_reg(anatomical_brain, reference_brain):

    import os

    cmd = "antsRegistration " \
          "--collapse-output-transforms 0 " \
          "--dimensionality 3 " \
          "--initial-moving-transform [ %s, %s, 0 ] " \
          "--interpolation Linear " \
          "--output [ transform, transform_Warped.nii.gz ] " \
          "--transform Rigid[ 0.1 ] " \
          "--metric MI[ %s, %s, 1, 32, Regular, 0.25 ] " \
          "--convergence [ 1000x500x250x100, 1e-08, 10 ] " \
          "--smoothing-sigmas 3.0x2.0x1.0x0.0 " \
          "--shrink-factors 8x4x2x1 " \
          "--use-histogram-matching 1 " \
          "--transform Affine[ 0.1 ] " \
          "--metric MI[ %s, %s, 1, 32, Regular, 0.25 ] " \
          "--convergence [ 1000x500x250x100, 1e-08, 10 ] " \
          "--smoothing-sigmas 3.0x2.0x1.0x0.0 " \
          "--shrink-factors 8x4x2x1 " \
          "--use-histogram-matching 1 " % \
          (reference_brain, anatomical_brain, reference_brain, \
           anatomical_brain, reference_brain, anatomical_brain)


    os.system(cmd)


    warp_list = []

    files = [f for f in os.listdir('.') if os.path.isfile(f)]

    for f in files:

        if ("transform" in f) and ("Warped" not in f):
            warp_list.append(os.getcwd() + "/" + f)

        if "Warped" in f:
            warped_image = os.getcwd() + "/" + f


    return warp_list, warped_image



def separate_warps_list(warp_list, selection):

    for warp in warp_list:

        if selection in warp:

            selected_warp = warp
     
    return selected_warp



def pick_seg_type(probability_maps, seg_type):

    """
    Returns the selected probability map from the list of segmented
    probability maps

    Parameters
    ----------

    probability_maps : list (string)
        List of Probability Maps

    Returns
    -------

    file : string
        Path to segment_prob_0.nii.gz is returned

    """

    import os
    import sys


    if(isinstance(probability_maps, list)):

        if(len(probability_maps) == 1):

            probability_maps = probability_maps[0]

        for filename in probability_maps:

            if seg_type == "csf":

                if filename.endswith("_0.nii.gz"):
                    return filename

            elif seg_type == "gm":

                if filename.endswith("_1.nii.gz"):
                    return filename

            elif seg_type == "wm":

                if filename.endswith("_2.nii.gz"):
                    return filename


    return None



