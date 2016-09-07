
def gather_raw_data(site_folder, include_sites=False, subject_inclusion=None):

    import os
    
    sub_dict = {}
    inclusion_list = []
    
    # create subject inclusion list
    if subject_inclusion:
        with open(subject_inclusion, "r") as f:
            inclusion_list = f.readlines()
        # remove any /n's
        inclusion_list = map(lambda s: s.strip(), inclusion_list)      
    
    for root, folders, files in os.walk(site_folder):
    
        for filename in files:
        
            fullpath = os.path.join(root, filename)
            
            if ".nii" in fullpath:
            
                # /path_to_site_folder/subject_id/session_id/scan_id/..
                second_half = fullpath.split(site_folder)[1]
                second_half_list = second_half.split("/")

                try:
                    second_half_list.remove("")
                except:
                    pass

                if include_sites:
                    try:
                        site_id = second_half_list[-5]
                        subject_id = second_half_list[-4]
                        session_id = second_half_list[-3]
                        scan_id = second_half_list[-2]
                    except:
                        err = "\n\n[!] Could not parse the data directory " \
                              "structure for this file - is it in the " \
                              "correct format?\nFile path:\n%s\n\nIt should "\
                              "be something like this:\n/site_folder/subject"\
                              "_id/session_id/scan_id/file.nii.gz\n\n" \
                              % second_half
                        print err
                else:
                    try:
                        subject_id = second_half_list[-4]
                        session_id = second_half_list[-3]
                        scan_id = second_half_list[-2]
                    except:                    
                        err = "\n\n[!] Could not parse the data directory " \
                              "structure for this file - is it in the " \
                              "correct format?\nFile path:\n%s\n\nIt should "\
                              "be something like this:\n/site_folder/subject"\
                              "_id/session_id/scan_id/file.nii.gz\n\n" \
                              % second_half
                        print err

                if subject_inclusion == None:
                    inclusion_list.append(subject_id)
                               
                #sub_info = (subject_id, session_id, scan_id)
                if ("anat" in scan_id) or ("anat" in filename) or \
                    ("mprage" in filename):
                    resource = "anatomical_scan"
                    
                if ("rest" in scan_id) or ("rest" in filename) or \
                    ("func" in scan_id) or ("func" in filename):
                    resource = "functional_scan"
                
                if (resource == "anatomical_scan") and \
                        (subject_id in inclusion_list):

                    if subject_id not in sub_dict.keys():
                        sub_dict[subject_id] = {}
                    
                    if session_id not in sub_dict[subject_id].keys():
                        sub_dict[subject_id][session_id] = {}
                    
                    if resource not in sub_dict[subject_id][session_id].keys():
                        sub_dict[subject_id][session_id][resource] = {}

                        if include_sites:
                            sub_dict[subject_id][session_id]["site_name"] = \
                                site_id
                                                       
                    if scan_id not in sub_dict[subject_id][session_id][resource].keys():
                        sub_dict[subject_id][session_id][resource][scan_id] = fullpath
                        
                if (resource == "functional_scan") and \
                         (subject_id in inclusion_list):
                
                    if subject_id not in sub_dict.keys():
                        sub_dict[subject_id] = {}
                    
                    if session_id not in sub_dict[subject_id].keys():
                        sub_dict[subject_id][session_id] = {}
                    
                    if resource not in sub_dict[subject_id][session_id].keys():
                        sub_dict[subject_id][session_id][resource] = {}
                        
                        if include_sites:
                            sub_dict[subject_id][session_id]["site_name"] = \
                                site_id
                                    
                    if scan_id not in sub_dict[subject_id][session_id][resource].keys():
                        sub_dict[subject_id][session_id][resource][scan_id] = fullpath
 
    return sub_dict


def pull_S3_sublist(bucket_name, bucket_prefix, creds_path):

    import os
    from indi_aws import fetch_creds

    s3_list = []
    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    # Build S3-subjects to download
    for bk in bucket.objects.filter(Prefix=bucket_prefix):
        s3_list.append(str(bk.key))

    return s3_list


def create_subdict_from_s3_list(s3_list, bucket_prefix, session_list=None,
    series_list=None, BIDS=False):

    s3_dict = {}

    # Read in series_list, if it is provided
    if session_list:
        try:
            session_list = os.path.abspath(session_list)
            with open(session_list,"r") as f:
                sessions = f.readlines()
            sessions = [i.rstrip("\n") for i in sessions]
        except:
            err = "\n\nCould not successfully read the session list.\n%s" \
                  % session_list
            raise Exception(err)

    # Read in series_list, if it is provided
    if series_list:
        try:
            series_list = os.path.abspath(series_list)
            with open(series_list,"r") as f:
                series = f.readlines()
            series = [i.rstrip("\n") for i in series]
        except:
            err = "\n\nCould not successfully read the series list.\n%s" \
                  % series_list
            raise Exception(err)

    # Build dictionary of filepaths
    for sfile in s3_list:

        include = False

        if BIDS:
            ssplit = sfile.split("/")
            folder = ssplit[-5]
            sub_id = ssplit[-4]
            session_id = ssplit[-3]
            scan_type = ssplit[-2]
            filename = ssplit[-1]

            if ".nii" not in filename:
                continue

            scan_id = None

            if sub_id not in filename:
                err = "\n\n[!] The dataset provided is not in BIDS format!" \
                      "\n\n"
                raise Exception(err)
            elif sub_id in filename:
                scan_id = filename.replace(sub_id,"")
            if session_id in scan_id:
                scan_id = scan_id.replace(session_id,"")
            if ".nii" in scan_id:
                scan_id = scan_id.replace(".nii","")
            if ".gz" in scan_id:
                scan_id = scan_id.replace(".gz","")
            if "task-" in scan_id:
                scan_id = scan_id.replace("task-","")
            if "__" in scan_id:
                scan_id = scan_id.replace("__","")          

            if scan_type == "anat":
                # this requirement is subject to change with the BIDS spec!!
                if "T1w.nii" in filename:
                    include = True
            if scan_type == "func":
                include = True

            if ("sub-" not in sub_id) or (scan_id == None):
                err = "\n\n[!] This is not a BIDS-formatted dataset!\n\n"
                raise Exception(err)

            # this is to get the folder containing the sub-ID folders so that
            # we can avoid sub-folders within the directory
            bucket_prefix_split = bucket_prefix.split("/")
            if bucket_prefix_split[-1] != "":
                containing_folder = bucket_prefix_split[-1]
            else:
                containing_folder = bucket_prefix_split[-2]

            if containing_folder != folder:
                include = False

            if session_list:
                if session_id not in sessions:
                    include = False

            if series_list:
                for series_id in series:
                    if series_id in scan_id:
                        break
                else:
                    include = False 

        else:
            ssplit = sfile.split('/')
            sub_id = ssplit[-4]
            session_id = ssplit[-3]
            scan_id = ssplit[-2]
            filename = ssplit[-1]

            if ("anat" in scan_id) or ("anat" in filename) or \
                ("mprage" in filename):
                scan_type = "anat"
                include = True
            if ("func" in scan_id) or ("rest" in scan_id) or \
                ("func" in filename) or ("rest" in filename):
                scan_type = "func"
                include = True

            if session_list:
                if session_id not in sessions:
                    include = False

            if series_list:
                if scan_id not in series:
                    include = False

        if (include == True) and ("nii" in filename):
        
            if scan_type == "anat":
                subkey_type = "anatomical_scan"
            elif scan_type == "func":
                subkey_type = "functional_scan"

            resource_dict = {}
            resource_dict[subkey_type] = sfile

            # this ONLY handles raw data inputs, not CPAC-generated outputs!
            if not s3_dict.has_key((sub_id, session_id, scan_id)):

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

    dict_len = len(s3_dict)
    print "Total number of subject-session-scans: %d\n" % dict_len

    return s3_dict


def write_inputs_dict_to_yaml_file(input_dict, yaml_outpath):

    import os
    import yaml         

    yaml_outpath = os.path.abspath(yaml_outpath)
    if (".yml" not in yaml_outpath) and (".yaml" not in yaml_outpath):
        yaml_outpath += ".yml"

    # write yaml file
    try:
        with open(yaml_outpath,"wt") as f:
            f.write(yaml.dump(input_dict))
    except:
        err = "\n\n[!] Error writing YAML file output.\n1. Do you have " \
              "write permissions for the output path provided?\n2. Did you " \
              "provide a full path for the output path? Example: /home/data" \
              "/sublist.yml\n\nOutput path provided: %s\n\n" % yaml_outpath
        raise Exception(err)
        
    if os.path.isfile(yaml_outpath):
        print "\nInputs dictionary file successfully created: %s\n" \
              % yaml_outpath
    else:
        err = "\n[!] Filepaths from the have not been successfully " \
              "saved to the YAML file!\nOutput filepath: %s\n" \
              % yaml_outpath
        raise Exception(err)