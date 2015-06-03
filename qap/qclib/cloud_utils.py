# cloud_utils.py
#
# Contributing authors: Daniel Clark, Steve Giavasis, 2015

'''
'''

# Function to build, select, and download subject data
def dl_subj_from_s3(subj_idx, img_type, creds_path):
    '''
    '''

    # Import packages
    from CPAC.AWS import fetch_creds, aws_utils

    # Init variables
    bucket_prefix = 'data/Projects/ABIDE_Initiative/RawData'
    local_prefix ='/mnt/inputs'
    bucket = fetch_creds.return_bucket(creds_path, 'fcp-indi')
    s3_list = []
    s3_dict = {}

    # Filter for anat/rest
    if img_type == 'anat':
        subkey_type = 'anatomical_scan'
    elif img_type == 'rest':
        subkey_type = 'functional_scan'

    # Build S3-subjects to download
    for bk in bucket.list(prefix=bucket_prefix):
        s3_list.append(str(bk.name))

    print 's3 list:'
    print s3_list

    # Build dictionary of filepaths
    for sfile in s3_list:
        ssplit = sfile.split('/')
        key = '_'.join(ssplit[5:7])
        sub_key = ssplit[-2]
        if img_type in sub_key:
            sub_key = subkey_type
        else:
            continue
        #elif 'rest' in sub_key:
        #    sub_key = 'functional_scan'
        if not s3_dict.has_key(key):
            s3_dict[key] = {}
            s3_dict[key][sub_key] = sfile
        else:
            s3_dict[key][sub_key] = sfile

    if len(s3_dict) == 0:
        err = "\n[!] Filepaths have not been successfully gathered from the S3 bucket!\n"
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

    aws_utils.s3_download(bucket, s3_dl, local_prefix=local_prefix, bucket_prefix=bucket_prefix)

    sub_dict = {subj_key : sub_dict}

    # Return single subject dictionary
    return sub_dict


def upl_qap_output(cfg_file, creds_path):
    '''
    '''

    # Import packages
    from CPAC.AWS import aws_utils, fetch_creds
    import os
    import yaml

    # Init variables
    bucket = fetch_creds.return_bucket(creds_path, 'fcp-indi')

    # Load config file
    cfg_dict = yaml.load(open(cfg_file, 'r'))
    output_dir = cfg_dict['output_directory']

    # And upload data
    upl_files = []
    for root, dirs, files in os.walk(output_dir):
        if files:
            upl_files.extend([os.path.join(root, fil) for fil in files])

    # Using CPAC AWS utils
    s3_upl_files = [ufile.replace(output_dir, 'data/Projects/ABIDE_Initiative/Outputs/qap') \
                   for ufile in upl_files]
    aws_utils.s3_upload(bucket, upl_files, s3_upl_files)

