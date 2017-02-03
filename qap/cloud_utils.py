# cloud_utils.py
#
# Contributing authors: Daniel Clark, Steve Giavasis, 2015

def dl_subj_from_s3(subj_idx, cfg_file, s3_dict_yaml):
    """Download a single participant's data from the Amazon S3 bucket.

    :type subj_idx: int
    :param subj_idx: The participant's index within the list of S3 bucket
                     filepaths.
    :type cfg_file: str
    :param cfg_file: Filepath to the pipeline configuration file containing
                     S3 bucket and AWS credentials information.
    :type s3_dict_yaml: str
    :param s3_dict_yaml: Filepath to the YAML file containing the AWS S3
                         filepaths.
    :rtype: dict
    :return: A dictionary with one entry mapping the participant's ID info
             to local filepath of the newly-downloaded data file.
    """

    # Import packages
    from indi_aws import fetch_creds, aws_utils
    import yaml

    # Load config file
    with open(cfg_file,'r') as f:
        cfg_dict = yaml.load(f)

    # Init variables
    bucket_prefix = cfg_dict["bucket_prefix"]
    local_prefix = cfg_dict["local_prefix"]
    bucket_name = cfg_dict["bucket_name"]
    creds_path = cfg_dict["creds_path"]
    
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    s3_list = []
    s3_dict = {}

    # pull in S3 dict yaml
    with open(s3_dict_yaml,'r') as f:
        s3_dict = yaml.load(f)

    if len(s3_dict) == 0:
        err = "\n[!] Filepaths have not been successfully gathered from " \
              "the filepath YAML dictionary!\n"
        raise Exception(err)

    # Get list of subject keys for indexing
    sd_keys = s3_dict.keys()
    sd_keys.sort()

    # Grab subject dictionary of interest
    subj_key = sd_keys[subj_idx-1]
    sub_dict = s3_dict[subj_key]

    # Download subject data to local prefix
    s3_dl = []
    for s3_key, s3_path in sub_dict.items():
        s3_dl.append(s3_path)
        sub_dict[s3_key] = s3_path.replace(bucket_prefix, local_prefix)

    local_dl = [s.replace(bucket_prefix, local_prefix) for s in s3_dl]

    aws_utils.s3_download(bucket, (s3_dl,local_dl))

    sub_dict = {subj_key : sub_dict}

    # Return single subject dictionary
    return sub_dict


def download_single_s3_path(s3_path, cfg_dict):
    """Download a single file from an AWS s3 bucket.

    :type s3_path: str
    :param s3_path: An "s3://" pre-pended path to a file stored on an
                    Amazon AWS s3 bucket.
    :type cfg_dict: dictionary
    :param cfg_dict: A dictionary containing the pipeline setup
                     parameters.
    :rtype: str
    :return: The local filepath of the downloaded s3 file.
    """

    import os
    from indi_aws import fetch_creds, aws_utils
    from workflow_utils import raise_smart_exception

    # Init variables
    working_dir = cfg_dict["working_directory"]
    creds_path = cfg_dict["creds_path"]

    if "s3://" in s3_path:
        s3_prefix = s3_path.replace("s3://","")
    else:
        err = "[!] S3 filepaths must be pre-pended with the 's3://' prefix."
        raise_smart_exception(locals(),err)

    bucket_name = s3_prefix.split("/")[0]
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    data_dir = s3_path.split(bucket_name + "/")[1]
    local_dl = os.path.join(working_dir, data_dir)

    if os.path.isfile(local_dl):
        print "\nS3 bucket file already downloaded! Skipping download."
        print "S3 file: %s" % s3_path
        print "Local file already exists: %s\n" % local_dl
    else:
        aws_utils.s3_download(bucket, ([data_dir], [local_dl]))

    return local_dl


def upl_qap_output(cfg_file):
    """Upload a pipeline output file to an AWS S3 bucket.

    :type cfg_file: str
    :param cfg_file: Filepath to the pipeline configuration file containing
                     S3 bucket and AWS credentials information.
    """

    # Import packages
    from indi_aws import aws_utils, fetch_creds
    import os
    import yaml

    # Load config file
    with open(cfg_file["pipeline_config_yaml"],'r') as f:
        cfg_dict = yaml.load(f)

    # Init variables
    bucket_name = cfg_dict["bucket_name"]
    bucket_out_prefix = cfg_dict["bucket_prefix"]
    creds_path = cfg_dict["creds_path"]
    
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)
        
    output_dir = cfg_dict['output_directory']

    # And upload data
    upl_files = []
    for root, dirs, files in os.walk(output_dir):
        if files:
            upl_files.extend([os.path.join(root, fil) for fil in files])

    # Using INDI AWS utils
    s3_upl_files = [ufile.replace(output_dir, bucket_out_prefix) \
                   for ufile in upl_files]

    aws_utils.s3_upload(bucket, (upl_files, s3_upl_files))

