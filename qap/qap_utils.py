
def create_expr_string(clip_level_value):
    """Create the expression arg string to run AFNI 3dcalc via Nipype.

    :type clip_level_value: int
    :param clip_level_value: The integer of the clipping threshold.
    :rtype: str
    :return The string intended for the Nipype AFNI 3dcalc "expr" arg inputs.
    """

    expr_string = "step(a-%s)" % clip_level_value

    return expr_string


def read_nifti_image(nifti_infile):
    """Read a NIFTI file into Nibabel-format image data.

    :type nifti_infile: str
    :param nifti_infile: The filepath of the NIFTI image to read in.
    :rtype: Nibabel image
    :return: Image data in Nibabel format.
    """

    import nibabel as nb
    from qap.qap_utils import raise_smart_exception

    try:
        nifti_img = nb.load(nifti_infile)
    except:
        err = "\n\n[!] Could not load the NIFTI image using Nibabel:\n" \
              "%s\n\n" % nifti_infile
        raise_smart_exception(locals(), err)

    return nifti_img


def write_nifti_image(nifti_img, file_path):
    """Write image data in Nibabel format into a NIFTI file.

    :type nifti_img: Nibabel image
    :param nifti_img: The image data Nibabel object to write out.
    :type file_path: str
    :param file_path: The filepath of the NIFTI image to create.
    """

    import nibabel as nb
    from qap.qap_utils import raise_smart_exception

    try:
        nb.save(nifti_img, file_path)
    except:
        err = "\n\n[!] Could not save the NIFTI image using Nibabel:\n" \
              "%s\n\n" % file_path
        raise_smart_exception(locals(), err)


def read_json(json_filename):
    """Read the contents of a JSON file.

    :type json_filename: str
    :param json_filename: The path to the JSON file.
    :rtype: dict
    :return: Dictionary containing the info from the JSON file.
    """

    import os
    import json
    from qap.qap_utils import raise_smart_exception

    if not os.path.exists(json_filename):
        err = "\n\n[!] The JSON file provided does not exist.\nFilepath: " \
              "%s\n\n" % json_filename
        raise_smart_exception(locals(),err)

    with open(json_filename, "r") as f:
        json_dict = json.load(f)

    return json_dict


def write_json(output_dict, json_file):
    """Either update or write a dictionary to a JSON file.

    :type output_dict: dict
    :param output_dict: The dictionary to write or append to the JSON file.
    :type json_file: str
    :param json_file: The filepath of the JSON file to write to or update.
    :rtype: str
    :return: Filepath of the JSON file written to.
    """

    import os
    import json
    from lockfile import FileLock

    from qap.qap_utils import read_json

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


def load_image(image_file):
    """Load a raw scan image from a NIFTI file and check it.

    :type image_file: str
    :param image_file: Path to the image, usually a structural or functional
                       scan.
    :rtype: Nibabel data
    :return: Image data in Nibabel format.
    """

    import nibabel as nib
    import numpy as np
    from qap.qap_utils import raise_smart_exception

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

    :type mask_file: str
    :param mask_file: Filepath to the binarized mask file.
    :type ref_file: str
    :param ref_file: Filepath to the anatomical file the mask is meant for.
    :rtype: Nibabel data
    :return: The mask data in Nibabel format.
    """

    import nibabel as nib
    import numpy as np

    from qap.qap_utils import raise_smart_exception

    try:
        mask_img = nib.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask_dat = mask_img.get_data()
    ref_img = nib.load(ref_file)

    # Check that the specified mask is binary.
    mask_vals = np.unique(mask_dat)
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
    # (have the same affine matrix)
    if (mask_img.get_affine() == ref_img.get_affine()).all == False:
        err = "Error: Mask and anatomical image are not in the same space " \
              "for %s vs %s" % (mask_file, ref_file)
        raise_smart_exception(locals(),err)

    return mask_dat


def create_anatomical_background_mask(anatomical_data, fg_mask_data, 
    exclude_zeroes=False):
    """Create a mask of the area outside the head in an anatomical scan by
    inverting a provided foreground mask.

    :type anatomical_data: NumPy array
    :param anatomical_data: An array of the raw anatomical data.
    :type fg_mask_data: NumPy array
    :param fg_mask_data: An array of binary foreground mask data.
    :type exclude_zeroes: bool
    :param exclude_zeroes: (default: False) Flag to exclude pure zero values
                           when creating the background mask.
    :rtype: Nibabel data
    :return bg_mask_data: Background mask data in Nibabel format.
    """

    from qap.qap_utils import raise_smart_exception

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


def raise_smart_exception(local_vars, msg=None):
    """Raise an exception with more information about the traceback, and
    enforce inclusion of the locals().

    :type local_vars: dict
    :param local_vars: Input for locals().
    :type msg: str
    :param msg: (default: None) The custom error message for the exception in
                question.
    """

    import traceback
    e = "\n\nLocal variables:\n%s\n\n%s\n\n" \
        % (str(local_vars), str(traceback.format_exc()))
    if msg:
        e = "%s\n\n%s\n\n" % (e, str(msg))
    raise Exception(e)


def check_input_resources(resource_pool, resource_name):
    """Check to make sure a specific resource/file is present in the
    resource pool.

    :type resource_pool: dict
    :param resource_pool: The resource pool of resources (which includes
                          output files of sub-workflows, and connection
                          pointers for Nipype nodes/workflows).
    :type resource_name: str
    :param resource_name: The name of the output/intermediary file to check
                          for within the resource pool.
    """

    import os

    if resource_name not in resource_pool.keys():
        err = "Resource pool: %s\n\n[!] The resource '%s' is missing " \
              "from the resource pool, and it is needed in one of the steps " \
              "of the pipeline. Please make sure it is specified " \
              "properly." % (resource_pool, resource_name)

        raise_smart_exception(locals(), err)

    else:
        if len(resource_pool[resource_name]) > 2:
            if not os.path.isfile(resource_pool[resource_name]):
                err = "[!] The path provided for the resource '%s' " \
                      "does not exist!\nPath provided: %s" % \
                      (resource_name, resource_pool[resource_name])

                raise_smart_exception(locals(), err)


def check_config_settings(config, parameter):
    """Check to make sure a configuration setting/parameter is present in the
    pipeline configuration dictionary.

    :type config: dict
    :param config: A dictionary keying configuration options to their chosen
                   selections.
    :type parameter: str
    :param parameter: The key of the configuration parameter to be checked.
    """

    if parameter not in config.keys():
        err = "[!] The parameter '%s' is missing from your pipeline " \
              "configuration .YML file. Please make sure this is specified " \
              "properly." % parameter
        raise_smart_exception(locals(), err)


def generate_nipype_workflow_graphs(workflow, out_dir=None):
    """Generate the Nipype workflow dependency graphs given the workflow
    object.

    :type workflow: Nipype workflow object
    :param workflow: The connected workflow object.
    :type out_dir: str
    :param out_dir: (default: None) The directory where to write the
                    dependency graph .dot and .png files to.
    """

    if not out_dir:
        pass

    """
    workflow.write_graph(
        dotfilename=op.join(config["output_directory"], \
                            "".join([run_name, ".dot"])),
        simple_form=False)
    workflow.write_graph(
        graph2use="orig",
        dotfilename=op.join(config["output_directory"], \
                            "".join([run_name, ".dot"])),
        simple_form=False)
    workflow.write_graph(
        graph2use="hierarchical",
        dotfilename=op.join(config["output_directory"], \
                            "".join([run_name, ".dot"])),
        simple_form=False)
    """