import os
import yaml
import json
import copy

import re


def bids_decode_filename(file_path, dbg=False):
    """
    Decode BIDS filename into a dictionary.

    :param file_path: path containing the filename to be decoded
    :param dbg: indicates whether or not debugging messages should be printed
    :return: dictionary containing information decoded from the filename
    """

    bids_dict = {}

    filename = os.path.basename(file_path)

    # first lets make sure that we know how to handle the file
    if 'nii' not in filename.lower() and 'json' not in filename.lower():
        raise IOError("File {0} does not appear to be a nifti or json file".format(filename))

    if dbg:
        print("Parsing BIDS information from {0}".format(file_path))

    site_index = -1

    # first figure out if there is a site directory level, this isn't
    # specified in BIDS currently, but hopefully will be in the future
    file_path_values = os.path.dirname(file_path).split('/')

    # if this is a derivative, then the site directory is just above the "derivative" directory
    if "derivatives" in file_path:
        site_index = file_path_values.index(b"derivatives") - 1

        if site_index + 2 > 0 and file_path_values[site_index + 2] and \
                "sub-" in file_path_values[site_index + 3]:
            # go ahead and include the pipeline
            bids_dict["pipeline"] = file_path_values[site_index + 2]

        else:
            raise ValueError("Could not find pipeline level for derivative {0}. Does the path conform to the "
                             "BIDS standard?".format(file_path))
    else:
        reg_result = re.search("(sub-[0-9]+)", file_path)
        if reg_result:
            site_index = file_path_values.index(reg_result.group(0)) - 1

    # site the site name
    if site_index != -1 and file_path_values[site_index]:
        bids_dict["site"] = file_path_values[site_index]
    else:
        bids_dict["site"] = "none"

    # remove filename extensions
    filename = filename.split(".")[0]

    # convert the filename string into a dictionary to pull out the other
    # key value pairs
    singleton_count = 0
    for key_val_pair in filename.split(b'_'):
        if "-" in key_val_pair:
            chunks = key_val_pair.split(b'-')
            bids_dict[chunks[0]] = "-".join([str(c) for c in chunks[1:]])
        else:
            if singleton_count == 0:
                bids_dict["scantype"] = key_val_pair
            elif singleton_count == 1:
                bids_dict["derivative"] = key_val_pair
            else:
                raise ValueError("Received {0} singletons in filename {1}, expecting 2 at most.".format(singleton_count,
                                                                                                        filename))
            singleton_count += 1

    if not bids_dict["scantype"]:
        raise ValueError("Filename ({0}) does not appear to contain scan type. "
                         "Does it conform to the BIDS format?".format(filename))

    if "derivative" in bids_dict and "derivatives" not in file_path:
        raise ValueError("{0} appears to be a derivative, but the path does not contain a 'derivatives' directory. "
                         "Does it conform to the BIDS format?".format(filename))

    if "derivatives" in file_path and "derivative" not in bids_dict:
        raise ValueError("Derivative ({0}) does not appear to contain a derivative tag. "
                         "Does it conform to the BIDS format?".format(filename))

    if "derivatives" in file_path and "pipeline" not in bids_dict:
        raise ValueError("Derivative ({0}) does not appear to contain a pipeline tag. "
                         "Does it conform to the BIDS format?".format(filename))

    if 'bold' in bids_dict["scantype"] and 'task' not in bids_dict:
        raise ValueError("Filename ({0}) is a BOLD file, but doesn't contain a task. "
                         "Does it conform to the BIDS format?".format(filename))

    bids_dict["site"] = re.sub('[\s\-_]+', '', bids_dict["site"])

    return bids_dict


def bids_retrieve_parameters(bids_parameters_dict, bids_dict, dbg=False):
    """

    Retrieve the BIDS parameters from bids_parameters_dict for BIDS file
    corresponding to bids_dict. If an exact match for bids_dict is not found
    the nearest match is returned, corresponding to the BIDS inheritance
    principle.

    :param bids_parameters_dict: BIDS configuration dictionary, this is a
      multi-level dictionary that maps the components of a bids filename
      (i.e. sub, ses, acq, run) to a dictionary that contains the BIDS
      parameters (RepetitionTime, EchoTime, etc). This information is
      extracted from sidecar json files using the principle of inheritance
      using the bids_parse_configs function
    :param bids_dict: Dictionary built from the name of a file in the BIDS
      format. This is built using the bids_decode_filename by splitting on
      "-" and "_" delimiters
    :param dbg: boolean flag that indicates whether or not debug statements
      should be printed, defaults to "False"
    :return: returns a dictionary that contains the BIDS parameters
    """

    params = {}

    t_dict = bids_parameters_dict  # pointer to current dictionary

    # try to populate the configuration using information
    # already in the list
    for level in ['scantype', 'site', 'sub', 'ses', 'task', 'acq',
                  'rec', 'run']:
        if level in bids_dict:
            key = "-".join([level, bids_dict[level]])
        else:
            key = "-".join([level, "none"])

        if dbg:
            print("Key: {0}".format(key))

        # if the key doesn't exist in the config dictionary, check to see if
        # the generic key exists and return that
        if key in t_dict:
            t_dict = t_dict[key]
        else:
            if dbg:
                print("Could find {0}, so going with {1}".format(key, "-".join([level, "none"])))
            key = "-".join([level, "none"])
            if key in t_dict:
                t_dict = t_dict[key]

    # if we have an image parameter dictionary at this level, use it to
    # initialize our configuration we look for "RepetitionTime", because
    # according to the spec it is a mandatory parameter for JSON
    # side car files

    if dbg:
        print(t_dict)

    for key in t_dict.keys():
        if u'RepetitionTime' in key:
            params = copy.deepcopy(t_dict)
            break

    return params


def bids_parse_sidecar(bids_dict, dbg=False):
    """
    Uses the BIDS principle of inheritance to build a data structure that
    maps parameters in side car .json files to components in the names of
    corresponding nifti files.

    :param bids_dict: dictionary that maps paths of sidecar json files
       (the key) to a dictionary containing the contents of the files (the values)
    :param dbg: boolean flag that indicates whether or not debug statements
       should be printed
    :return: a dictionary that maps parameters to components from BIDS file names
       such as sub, sess, run, acq, and scan type
    """

    # we are going to build a large-scale data structure, consisting of many
    # levels of dictionaries to hold the data.
    bids_parameters_dict = {}

    # initialize 'default' entries, this essentially is a pointer traversal
    # of the dictionary
    t_dict = bids_parameters_dict
    for level in ['scantype', 'site', 'sub', 'ses', 'task',
                  'acq', 'rec', 'run']:
        key = '-'.join([level, 'none'])
        t_dict[key] = {}
        t_dict = t_dict[key]

    if dbg:
        print(bids_parameters_dict)

    # get the paths to the json yaml files in config_dict, the paths contain
    # the information needed to map the parameters from the json files (the values
    # of the config_dict) to corresponding nifti files. We sort the list
    # by the number of path components, so that we can iterate from the outer
    # most path to inner-most, which will help us address the BIDS inheritance
    # principle
    parameter_file_path_list = sorted(bids_dict.keys(), key=lambda p: len(p.split('/')))

    if dbg:
        print("parameter_file_paths {}".format(", ".join(parameter_file_path_list)))
        print(bids_dict)

    for parameter_file_path in parameter_file_path_list:

        if dbg:
            print("processing {0}".format(parameter_file_path))

        # decode the filepath into its various components as defined by  BIDS
        filename_bids_dict = bids_decode_filename(parameter_file_path)

        # handling inheritance is a complete pain, we will try to handle it by
        # build the key from the bottom up, starting with the most
        # parsimonious possible, incorporating configuration information that
        # exists at each level

        # first lets try to find any parameters that already apply at this
        # level using the information in the json's file path
        t_params = bids_retrieve_parameters(bids_parameters_dict, filename_bids_dict)

        # now populate the parameters
        bids_config = {}
        if t_params:
            bids_config.update(t_params)

        # add in the information from this config file
        t_config = bids_parameters_dict[parameter_file_path]
        while isinstance(t_config, list):
            print("{0} contents are in a list?".format(parameter_file_path))
            t_config = t_config[0]

        if dbg:
            print("updating {0} with {1}".format(bids_config, t_config))

        bids_config.update(t_config)

        # now put the configuration in the data structure, by first iterating
        # to the location of the key, and then inserting it. When a key isn't
        # defined we use the "none" value. A "none" indicates that the
        # corresponding parameters apply to all possible settings of that key
        # e.g. run-1, run-2, ... will all map to run-none if no json files
        # explicitly define values for those runs
        t_dict = bids_parameters_dict  # pointer to current dictionary
        for level in ['scantype', 'site', 'sub', 'ses', 'task', 'acq',
                      'rec', 'run']:
            if level in filename_bids_dict:
                key = "-".join([level, filename_bids_dict[level]])
            else:
                key = "-".join([level, "none"])

            if key not in t_dict:
                t_dict[key] = {}

            t_dict = t_dict[key]

        t_dict.update(bids_config)

    return bids_parameters_dict


def bids_generate_qap_data_configuration(bids_dir, paths_list, configuration_dictionary=None, credentials_path="",
                                         debug=False):
    """
    Generates a QAP formatted subject list from information contained in a
    BIDS formatted set of data.

    :param bids_dir:  base directory that contains all of the data, this could be
       a directory that contains data for a multiple BIDS data sets, in which
       case the intervening directories will be interpreted as site names
    :param paths_list: lists of all nifti files found in bids_dir, these paths
       are relative to bids_dir
    :param configuration_dictionary: a dictionary containing BIDS parameter
       information
    :param credentials_path: if using S3 bucket, this path credentials needed to
       access the bucket, if accessing anonymous bucket, this can be set
       to None
    :param debug: boolean indicating whether or not the debug statements should
       be printed
    :return: a data configuration dictionary suitable for use by QAP to
    """

    if debug:
        print("gen_bids_sublist called with:")
        print("  bids_dir: {0}".format(bids_dir))
        print("  # paths: {0}".format(str(len(paths_list))))
        print("  configuration_dictionary: {0}".format("missing" if not configuration_dictionary else "found"))
        print("  credentials_path: {0}".format(credentials_path))

    # if configuration information is not desired, config_dict will be empty,
    # otherwise parse the information in the sidecar json files into a dict
    # we can use to extract data for our nifti files
    bids_configuration_dictionary = {}
    if configuration_dictionary:
        bids_configuration_dictionary = bids_parse_sidecar(configuration_dictionary, debug)

    data_configuration = {}

    for file_path in paths_list:
        file_path = file_path.rstrip()
        file_name = os.path.basename(file_path)

        if file_name.endswith(".nii") or file_name.endswith(".nii.gz") or file_name.endswith(".json"):

            # 1. decode file path
            file_bids_dict = bids_decode_filename(file_path)

            if "ses" not in file_bids_dict:
                file_bids_dict["ses"] = "1"

            if "sub" not in file_bids_dict:
                raise IOError('"sub" not found in {}, perhaps it is not in BIDS format?'.format(file_path))

            # 2. construct site, participant, and session keys
            participant = "-".join(["sub", file_bids_dict["sub"]])
            session = "-".join(["ses", file_bids_dict["ses"]])
            site = "-".join(["site", file_bids_dict["site"]])

            # 3. construct scan_name to be used as index in data configuration
            scan_name = ''

            if "bold" in file_bids_dict["scantype"]:

                if 'task' not in file_bids_dict:
                    raise IOError('"task" not found in {}, perhaps it is not correct BIDS format for a'
                                  ' functional file?'.format(file_path))

                scan_name = "-".join(["task", file_bids_dict["task"]])

            for bids_key in ["run", "acq"]:
                if bids_key in file_bids_dict:
                    if scan_name:
                        scan_name = "_".join([scan_name, "-".join([bids_key, file_bids_dict[bids_key]])])
                    else:
                        scan_name = "-".join([bids_key, file_bids_dict[bids_key]])

            if scan_name:
                scan_name = "_".join([scan_name, file_bids_dict["scantype"]])
            else:
                scan_name = file_bids_dict["scantype"]

            # 4. Create a resource_name corresponding to this file
            resource_name = ''

            if "derivative" not in file_bids_dict:
                if "T1w" in file_bids_dict["scantype"]:
                    resource_name = 'anatomical_scan'
                elif "bold" in file_bids_dict["scantype"]:
                    resource_name = 'functional_scan'
                else:
                    if debug:
                        print("Do not know how to map {0} to a resource, skipping".format(file_path))
                        continue
            else:
                for bids_key in ["pipeline", "variant", "space", "res", "atlas", "roi"]:
                    if bids_key in file_bids_dict:
                        if resource_name:
                            resource_name = "_".join([resource_name, "-".join([bids_key, file_bids_dict[bids_key]])])
                        else:
                            resource_name = "-".join([bids_key, file_bids_dict[bids_key]])

                if resource_name:
                    resource_name = "_".join([resource_name, file_bids_dict["derivative"]])
                else:
                    resource_name = file_bids_dict["derivative"]

            # 5. create a temp resource pool with the info for this file
            resource_dict = {resource_name: file_path}

            # if not a derivative, populate the resource pool with information from bids sidecar files, if they exist
            if "derivative" not in file_bids_dict and bids_configuration_dictionary:
                file_params = bids_retrieve_parameters(bids_configuration_dictionary, file_bids_dict)
                if not file_params:
                    print("Did not receive any parameters for {0}, is this a problem?".format(file_path))
                else:
                    resource_dict.update(file_params)

            # 6. Add the resource pool into the data configuration
            if site not in data_configuration:
                data_configuration[site] = {}

            if participant not in data_configuration[site]:
                data_configuration[site][participant] = {}

            if session not in data_configuration[site][participant]:
                data_configuration[site][participant][session] = {}

            if scan_name not in data_configuration[site][participant][session]:
                data_configuration[site][participant][session][scan_name] = {}

            data_configuration[site][participant][session][scan_name].update(resource_dict)

    return data_configuration


def collect_bids_files_configs(bids_dir, aws_input_credentials=None,
                               raise_error=True):
    """

    :param bids_dir:
    :param aws_input_credentials:
    :return: raw_file_paths = list of paths to raw files
             derivative_file_paths = list of paths to derivatives
             parameter_dict = dictionary containing parameters for the various files
    """

    raw_file_paths = []
    derivative_file_paths = []
    parameter_dict = {}

    if bids_dir.lower().startswith("s3://"):
        # s3 paths begin with s3://bucket/
        bucket_name = bids_dir.split('/')[2]
        s3_prefix = '/'.join(bids_dir.split('/')[:3])
        prefix = bids_dir.replace(s3_prefix, '').lstrip('/')

        if aws_input_credentials and not os.path.isfile(aws_input_credentials):
            raise IOError("Could not find aws input credentials ({0})".format(aws_input_credentials))

        from indi_aws import fetch_creds
        bucket = fetch_creds.return_bucket(aws_input_credentials, bucket_name)

        print("gathering files from S3 bucket {0} for {1}".format(bucket, prefix))

        for s3_obj in bucket.objects.filter(Prefix=prefix):
            # we only know how to handle T1w and BOLD files, for now
            if 'T1w' in str(s3_obj.key) or 'bold' in str(s3_obj.key):
                if str(s3_obj.key).endswith("json") and 'derivatives' not in str(s3_obj.key):
                    try:
                        parameter_dict[s3_obj.key.replace(prefix, "").lstrip('/')] \
                            = json.loads(s3_obj.get()["Body"].read())
                    except Exception as e:
                        print("Error retrieving {0} ({1})".format(s3_obj.key.replace(prefix, ""), e.message))
                        raise
                        
                elif str(s3_obj.key).endswith("json") or str(s3_obj.key).endswith("nii") or \
                        str(s3_obj.key).endswith("nii.gz"):

                    # keep track of raw data and derivatives separately, in case it matters to the user
                    file_path = os.path.join(s3_prefix, str(s3_obj.key))
                    if "derivatives" in file_path:
                        derivative_file_paths.append(file_path)
                    else:
                        raw_file_paths.append(file_path)

    else:
        for root, dirs, files in os.walk(bids_dir, topdown=False):
            if files:
                for f in files:
                    file_path = os.path.join(root, f)
                    if 'T1w' in f or 'bold' in f:
                        if f.endswith('nii') or f.endswith("nii.gz") or (f.endswith('json') and
                                                                         'derivatives' in file_path):
                            if "derivatives" in file_path:
                                derivative_file_paths.append(file_path)
                            else:
                                raw_file_paths.append(file_path)
                        elif f.endswith('json'):
                            file_contents = json.load(open(file_path, 'r'))
                            while isinstance(file_contents, list):
                                print("file ({}) contents are in a list?".format(file_path))
                                file_contents = file_contents[0]
                            parameter_dict.update({os.path.join(root.replace(bids_dir, '').lstrip('/'), f):
                                                   file_contents})

                        # don't know what other files are, so just drop then

    if not raw_file_paths and not parameter_dict and not derivative_file_paths:
        if raise_error:
            raise IOError("Didn't find any files in {}. Please verify that "
                          "the path is typed correctly, that you have read "
                          "access to the directory, and that it is not "
                          "empty.".format(bids_dir))

    return raw_file_paths, derivative_file_paths, parameter_dict


def write_data_configuration(data_configuration_filename, data_configuration_dictionary):

    if not data_configuration_filename.endswith('yml'):
        data_configuration_filename += '.yml'

    with open(data_configuration_filename, 'w') as ofd:
        yaml.dump(data_configuration_dictionary, ofd, encoding='utf-8')


def test_bids_decode(bids_dir, dbg=True):
  
    (bids_files, derivative_files, config) = collect_bids_files_configs(bids_dir, [])
    if dbg:
        print("Found %d config files for %d image files and %d derivatives" % (len(config), len(bids_files),
                                                                               len(derivative_files)))
    bids_files = bids_files + derivative_files

    for bids_file in bids_files:
        print("{0} :: {1}".format(bids_file, bids_decode_filename(bids_file)))

    return

def test_gen_bids_sublist_qap(bids_dir, test_yml, aws_input_credentials=None, debug=False, cfg=False):

    (bids_files, derivative_files, config) = collect_bids_files_configs(bids_dir,
                                                                        aws_input_credentials=aws_input_credentials)
    if debug:
        num_json = 0
        for df in derivative_files:
            if df.endswith("json"):
                num_json += 1
        print("Found {0} config files for {1} image files and {2} derivatives ({3} json derivatives)".format(
            len(config), len(bids_files), len(derivative_files), num_json))

    bids_files = bids_files + derivative_files

    if cfg:
        data_configuration = bids_generate_qap_data_configuration(bids_dir, bids_files, config,
                                                                  credentials_path=aws_input_credentials, debug=debug)
    else:
        data_configuration = bids_generate_qap_data_configuration(bids_dir, bids_files,
                                                                  credentials_path=aws_input_credentials, debug=debug)

    with open(test_yml, "w") as ofd:
        yaml.dump(data_configuration, ofd, encoding='utf-8')

    assert data_configuration

if __name__ == '__main__':

    test_gen_bids_sublist_qap("/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/bids_dir",
                              "/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/qap_test.yml",
                              debug=True)
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/data/ADHD200/RawDataBIDS/",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/test/rs_subject_list.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)

    # test_bids_decode("/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/bids_dir")

    #
    # test_gen_bids_sublist(
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/data/ADHD200/RawDataBIDS/",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/test/rs_subject_list.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
    # test_gen_bids_sublist(
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/data/ADHD200/RawDataBIDS/Peking_3",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/test/rs_subject_list_pk3.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
    test_gen_bids_sublist_qap(
        "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/",
        "/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/rs_subject_list_s3.yml",
        aws_input_credentials="/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
        debug=False)
    #
    # test_gen_bids_sublist(
    #     "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/Peking_3",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC/test/"
    #        "rs_subject_list_pk3_s3.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
    # test_gen_bids_sublist(
    #     "s3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC/test/"
    #        "rs_subject_list_corr_bmb3_s3.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
