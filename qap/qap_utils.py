
def load_image(image_file):
    """Load a raw scan image from a NIFTI file and check it.

    Keyword arguments:
      image_file -- path to the image, usually a structural or functional scan

    Returns:
      dat -- image data in Nibabel format
    """

    import nibabel as nib
    import numpy as np

    from workflow_utils import raise_smart_exception

    try:
        img = nib.load(image_file)
    except:
        raise_smart_exception(locals())

    dat = img.get_data()

    # Ensure that data is cast as at least 32-bit
    if np.issubdtype(dat.dtype, float):
        dat = dat.astype('float32')
        # Check for negative values
        if (dat < 0).any():
            print "found negative values, setting to zero (see file: %s)" \
                  % image_file
            dat[dat<0] = 0

    elif np.issubdtype(dat.dtype, int):
        dat = dat.astype('int32')

    elif np.issubdtype(dat.dtype, np.uint8):
        dat = dat.astype(np.uint8)

    else:
        msg = "Error: Unknown datatype %s" % dat.dtype
        raise_smart_exception(locals(),msg)

    return dat


def load_mask(mask_file, ref_file):
    """Load a mask from a NIFTI file and check the shape and dimensions.

    Keyword arguments:
      mask_file -- binarized mask file
      ref_file -- anatomical file the mask is meant for

    Returns:
      mask_dat -- the mask data in Nibabel format
    """

    import nibabel as nib
    import numpy as np

    from workflow_utils import raise_smart_exception

    try:
        mask_img = nib.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask_dat = mask_img.get_data()
    ref_img = nib.load(ref_file)

    # Check that the specified mask is binary.
    mask_vals   = np.unique(mask_dat)
    if (mask_vals.size != 2) or not (mask_vals == [0, 1]).all():
        err = "Error: Mask is not binary, has %i unique val(s) of %s " \
              "(see file %s)" % (mask_vals.size, mask_vals, mask_file)
        raise_smart_exception(locals(),err)

    # Verify that the mask and anatomical images have the same dimensions.
    if ref_img.shape != mask_img.shape:
        err = "Error: Mask and anatomical image are different dimensions " \
              "for %s" % mask_file
        raise_smart_exception(locals(),err)

    # Verify that the mask and anatomical images are in the same space
    # (have the samme affine matrix)
    if (mask_img.get_affine() == ref_img.get_affine()).all == False:
        err = "Error: Mask and anatomical image are not in the same space " \
              "for %s vs %s" % (mask_file, ref_file)
        raise_smart_exception(locals(),err)

    return mask_dat


def create_anatomical_background_mask(anatomical_data, fg_mask_data, 
    exclude_zeroes=False):
    """Create a mask of the area outside the head in an anatomical scan by
    inverting a provided foreground mask.

    Keyword arguments:
      anatomical_data -- a NumPy array of the raw anatomical data
      fg_mask_data -- a NumPy array of the binary foreground mask data
      exclude_zeroes -- (default: False) flag to exclude pure zero values when
                        creating the background mask

    Returns:
      bg_mask_data -- background mask data in Nibabel format
    """

    import numpy as np
    from workflow_utils import raise_smart_exception

    # invert the foreground mask
    try:
        bg_mask_data = 1 - fg_mask_data
    except Exception as e:
        err = "\n\n[!] Input data must be a NumPy array object, and not a " \
              "list.\n\nError details: %s\n\n" % e
        raise_smart_exception(locals(),err)

    if exclude_zeroes:
        # modify the mask to exclude zeroes in the background of the
        # anatomical image, as these are often introduced artificially and can
        # skew the QAP metric results
        bool_anat_data = anatomical_data > 0
        bg_mask_data = bg_mask_data * bool_anat_data

    return bg_mask_data

    
def json_to_csv(json_file, csv_output_dir=None):
    """Extract the data from the JSON output file and write it to a CSV file.

    Keyword arguments:
      json_file -- filepath to the JSON file to be written to a CSV
      csv_output_dir -- (default: None) path to the directory to write the CSV
                        file into

    Returns:
      csv_file -- the CSV file path
    """

    import os
    import json
    import pandas as pd
    from lockfile import FileLock

    from qap.qap_workflows_utils import read_json
    from qap.workflow_utils import raise_smart_exception

    qap_types = ["anatomical_spatial",
                 "functional_spatial",
                 "functional_temporal"]

    output_dict = {}

    json_dict = read_json(json_file)

    for sub_sess_scan in json_dict.keys():
        # flatten the JSON dict
        sub_json_dict = json_dict[sub_sess_scan]
        header_dict = {}
        qap_dict = {}

        try:
            header_dict = sub_json_dict["anatomical_header_info"]
        except KeyError:
            pass

        try:
            header_dict = sub_json_dict["functional_header_info"]
        except KeyError:
            pass

        for qap_type in qap_types:
            try:
                qap_dict = sub_json_dict[qap_type]
            except KeyError:
                continue

            for key in sub_json_dict.keys():
                if "anatomical" not in key and "functional" not in key:
                    qap_dict[key] = sub_json_dict[key]

            qap_dict.update(header_dict)

            try:
                output_dict[qap_type].append(qap_dict)
            except KeyError:
                output_dict[qap_type] = [qap_dict]

    for qap_type in output_dict.keys():

        json_df = pd.DataFrame(output_dict[qap_type])
        if not csv_output_dir:
            csv_output_dir = json_file.replace(json_file.split("/")[-1],"")
        csv_file = os.path.join(csv_output_dir, "qap_%s.csv" % qap_type)

        lock = FileLock(csv_file)
        lock.acquire()
        try:
            json_df.to_csv(csv_file)
        except:
            lock.release()
            err = "Could not write CSV file!\nCSV file: %s" % csv_file
            raise_smart_exception(locals(),err)
        lock.release()

    return csv_file
