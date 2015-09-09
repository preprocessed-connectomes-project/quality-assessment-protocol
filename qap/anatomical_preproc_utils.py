
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



