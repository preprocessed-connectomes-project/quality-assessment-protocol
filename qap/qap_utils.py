import json
import numpy
import nibabel
import traceback
import os

# encoder to help with exporting numpy data types and exceptions to json
class NumpyEncoder(json.JSONEncoder):

    def default(self, obj):

        if isinstance(obj, numpy.ndarray):
            return obj.tolist()
        if isinstance(obj, numpy.float32):
            return float(obj)
        if isinstance(obj, numpy.int64):
            return int(obj)
        if isinstance(obj, numpy.uint64):
            return int(obj)
        if isinstance(obj, Exception):
            return str(obj)

        return json.JSONEncoder.default(self, obj)


def raise_smart_exception(local_vars, msg=None):
    """Raise an exception with more information about the traceback, and
    enforce inclusion of the locals().

    :type local_vars: dict
    :param local_vars: Input for locals().
    :type msg: str
    :param msg: (default: None) The custom error message for the exception in
                question.
    """

    ll = "\n".join([":".join([k, str(v)]) for k, v in local_vars.items()])

    e = "\n\nLocal variables:\n{0}\n\n{1}\n\n".format(ll, str(traceback.format_exc()))
    if msg:
        e = "{0}\n\n{1}\n\n".format(e, str(msg))
    raise Exception(e)


def create_expr_string(clip_level_value):
    """Create the expression arg string to run AFNI 3dcalc via Nipype.

    :type clip_level_value: float
    :param clip_level_value: The desired clipping threshold
    :rtype: str
    :return The string intended for the Nipype AFNI 3dcalc "expr" arg inputs.
    """

    expr_string = "step(a-{0})".format(clip_level_value)

    return expr_string


def read_nifti_image(nifti_infile):
    """Read a NIFTI file into Nibabel-format image data.

    :type nifti_infile: str
    :param nifti_infile: The filepath of the NIFTI image to read in.
    :rtype: Nibabel image
    :return: Image data in Nibabel format.
    """
    nifti_image = None

    try:
        nifti_image = nibabel.load(nifti_infile)
    except: # TODO: figure out the right way to handle exceptions to avoid PEP8 problems
        err = "\n\n[!] Could not load the NIFTI image using Nibabel:\n{0}\n\n".format(nifti_infile)
        raise_smart_exception(locals(), err)

    return nifti_image


def write_nifti_image(nifti_img, file_path):
    """Write image data in Nibabel format into a NIFTI file.

    :type nifti_img: Nibabel image
    :param nifti_img: The image data Nibabel object to write out.
    :type file_path: str
    :param file_path: The filepath of the NIFTI image to create.
    """

    try:
        nibabel.save(nifti_img, file_path)
    except:
        err = "\n\n[!] Could not save the NIFTI image using Nibabel:\n{0}\n\n".format(file_path)
        raise_smart_exception(locals(), err)


def read_json(json_filename):
    """Read the contents of a JSON file.

    :type json_filename: str
    :param json_filename: The path to the JSON file.
    :rtype: dict
    :return: Dictionary containing the info from the JSON file.
    """

    json_dict = {}

    try:
        with open(json_filename, "r") as f:
            json_dict = json.load(f)
    except:
        err = "\n\n[!] Could not load JSON file {0}\n\n".format(json_filename)
        raise_smart_exception(locals(), err)

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
    from qap.qap_utils import NumpyEncoder

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
                except AttributeError:
                    current_dict[key] = output_dict[key]
    else:
        current_dict = output_dict

    if write:
        with open(json_file, "wt") as f:
            json.dump(current_dict, f, indent=2,
                      sort_keys=True, cls=NumpyEncoder)

    if os.path.exists(json_file):
        return json_file


def load_image(image_file, return_img=False):
    """Load a raw scan image from a NIFTI file and check it.

    :type image_file: str
    :param image_file: Path to the image, usually a structural or functional
                       scan.
    :type return_img: bool
    :param return_img: flag to indicate whether a nibabel.NiftiImage object containing the header information should
        be returned in addition to the voxel data
    :rtype: tuple, or numpy array
    :return: Depending on the value of return_img, returns either just the voxel data as a numpy array, or a tuple
        that contains both the numpy array and a nibabel.NiftiImage object containing the header information
    """

    nifti_image = None

    try:
        nifti_image = nibabel.load(image_file)
    except:
        raise_smart_exception(locals())

    voxel_data = nifti_image.get_data()

    # Ensure that data is cast as at least 32-bit
    if numpy.issubdtype(voxel_data.dtype, float):
        dat = voxel_data.astype('float32')
        # Check for negative values
        if (voxel_data < 0).any():
            print("found negative values, setting to zero (see file: %s)".format(image_file))
            dat[dat < 0] = 0

    elif numpy.issubdtype(voxel_data.dtype, int):
        dat = voxel_data.astype('int32')

    elif numpy.issubdtype(voxel_data.dtype, numpy.uint8):
        dat = voxel_data.astype(numpy.uint8)

    else:
        msg = "Error: Unknown datatype %s" % voxel_data.dtype
        raise_smart_exception(locals(), msg)

    if return_img:
        return voxel_data, nifti_image
    else:
        return voxel_data


def load_mask(mask_file, ref_file):
    """Load a mask from a NIFTI file and check the shape and dimensions.

    :type mask_file: str
    :param mask_file: Filepath to the binarized mask file.
    :type ref_file: str
    :param ref_file: Filepath to the anatomical file the mask is meant for.
    :rtype: Nibabel data
    :return: The mask data in Nibabel format.
    """

    nifti_image = None

    try:
        nifti_image = nibabel.load(mask_file)
    except:
        raise_smart_exception(locals())

    mask_data = nifti_image.get_data()
    ref_img = nibabel.load(ref_file)

    # Check that the specified mask is binary.
    mask_values = numpy.unique(mask_data)
    if (mask_values.size != 2) or not (mask_values == [0, 1]).all():
        err = "Error: Mask is not binary, has {0} unique val(s) of {1} (see file {2})".format(mask_values.size,
                                                                                              mask_values,
                                                                                              mask_file)
        raise_smart_exception(locals(), err)

    # Verify that the mask and reference images have the same dimensions.
    if ref_img.shape != nifti_image.shape:
        err = "Error: Mask and anatomical image are different dimensions for {0}".format(mask_file)
        raise_smart_exception(locals(), err)

    # Verify that the mask and reference images are in the same space
    # (have the same affine matrix)
    if not numpy.alltrue(nifti_image.get_affine() == ref_img.get_affine()):
        err = "Error: Mask and anatomical image are not in the same space for {0} vs {1}".format(mask_file, ref_file)
        raise_smart_exception(locals(), err)

    return mask_data


def create_anatomical_background_mask(anatomical_data, fg_mask_data, exclude_zeroes=False):
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

    # invert the foreground mask
    bg_mask_data = None

    try:
        bg_mask_data = 1 - fg_mask_data
    except Exception as e:
        err = "\n\n[!] Error calculating background mask. {0}\n\n".format(e)
        raise_smart_exception(locals(), err)

    if exclude_zeroes:
        # modify the mask to exclude zeroes in the background of the
        # anatomical image, as these are often introduced artificially and can
        # skew the QAP metric results
        bg_mask_data[numpy.isclose(anatomical_data, 0)] = 0

    return bg_mask_data


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

    if resource_name not in resource_pool.keys():
        err = "[!] Mandatory resource '{0}' is missing from resource pool ({1})".format(resource_pool, resource_name)
        raise_smart_exception(locals(), err)

    else:
        if len(resource_pool[resource_name]) > 2:
            if not os.path.isfile(resource_pool[resource_name]):
                err = "[!] Could not find resource '{0}' on the path provided ({1})!".format(resource_name,
                                                                                             resource_pool[
                                                                                                 resource_name])

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
        err = "[!] The parameter '%s' is missing from your pipeline configuration .YML file.".format(parameter)
        raise_smart_exception(locals(), err)


def generate_nipype_workflow_graphs(workflow, out_dir=None):
    """Generate the Nipype workflow dependency graphs given the workflow
    object.

    TODO: Why are the contents of this function comment out?

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
