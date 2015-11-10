# cloud_utils.py
#
# Contributing authors: Daniel Clark, Steve Giavasis, 2015

'''
'''

def pull_S3_sublist(yaml_outpath, img_type, cfg_file):

    # function example use:
    #
    # yamlpath = os.path.join(os.getcwd(), "s3dict.yml")
    #
    # # Build entire filepath dictionary from S3
    # s3_dict_yml = pull_S3_sublist(yamlpath, 'anat', args.config)

    import os
    from CPAC.AWS import fetch_creds
    import yaml

    s3_list = []
    s3_dict = {}

    # Load config file
    with open(cfg_file,'r') as f:
        cfg_dict = yaml.load(f)

    bucket_name = cfg_dict["bucket_name"]
    bucket_prefix = cfg_dict["bucket_prefix"]
    creds_path = cfg_dict["creds_path"]

    bucket = fetch_creds.return_bucket(creds_path, bucket_name)


    # Filter for anat/rest
    if img_type == 'anat':
        subkey_type = 'anatomical_scan'
    elif img_type == 'rest':
        subkey_type = 'functional_scan'


    # Build S3-subjects to download
    for bk in bucket.list(prefix=bucket_prefix):
        s3_list.append(str(bk.name))

    # Build dictionary of filepaths
    for sfile in s3_list:

        ssplit = sfile.split('/')

        sub_id = ssplit[-4]

        session_id = ssplit[-3]

        scan_id = ssplit[-2]

        if img_type in scan_id:
             
            # this ONLY handles raw data inputs, not CPAC-generated outputs!
            if not s3_dict.has_key((sub_id, session_id, scan_id)):

                resource_dict = {}
                resource_dict[subkey_type] = sfile

                s3_dict[(sub_id, session_id, scan_id)] = {}
                s3_dict[(sub_id, session_id, scan_id)].update(resource_dict)

            else:

                s3_dict[(sub_id, session_id, scan_id)].update(resource_dict)         

        else:

            continue
    
            
    if len(s3_dict) == 0:
        err = "\n[!] Filepaths have not been successfully gathered from " \
              "the S3 bucket!\n"
        raise Exception(err)
            
           
    # write yaml file
    with open(yaml_outpath,"wt") as f:
        f.write(yaml.dump(s3_dict))
        
    
    if os.path.isfile(yaml_outpath):
        return yaml_outpath
    else:
        err = "\n[!] Filepaths from the S3 bucket have not been " \
              "successfully saved to the YAML file!\nOutput filepath: %s\n" \
              % yaml_outpath
        raise Exception(err)



# Function to build, select, and download subject data
def dl_subj_from_s3(subj_idx, cfg_file, s3_dict_yaml):
    '''
    '''

    # Import packages
    from CPAC.AWS import fetch_creds, aws_utils
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

    aws_utils.s3_download(bucket, s3_dl, local_prefix=local_prefix, \
                              bucket_prefix=bucket_prefix)

    sub_dict = {subj_key : sub_dict}

    # Return single subject dictionary
    return sub_dict



def upl_qap_output(cfg_file):
    '''
    '''

    # Import packages
    from CPAC.AWS import aws_utils, fetch_creds
    import os
    import yaml

    # Load config file
    with open(cfg_file,'r') as f:
        cfg_dict = yaml.load(f)

    # Init variables
    bucket_name = cfg_dict["bucket_name"]
    bucket_out_prefix = cfg_dict["bucket_out_prefix"]
    creds_path = cfg_dict["creds_path"]
    
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)
        
    output_dir = cfg_dict['output_directory']

    # And upload data
    upl_files = []
    for root, dirs, files in os.walk(output_dir):
        if files:
            upl_files.extend([os.path.join(root, fil) for fil in files])

    # Using CPAC AWS utils
    s3_upl_files = [ufile.replace(output_dir, bucket_out_prefix) \
                   for ufile in upl_files]
    aws_utils.s3_upload(bucket, upl_files, s3_upl_files)



