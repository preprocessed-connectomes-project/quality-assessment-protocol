
def read_txt_file(txt_file):
    with open(txt_file,"r") as f:
        strings = f.read().splitlines()
    return strings


def gather_filepath_list(site_folder):

    # gathers all of the NIFTI files under the provided directory

    import os

    filepath_list = []
    for root, folders, files in os.walk(site_folder):
        for filename in files:
            fullpath = os.path.join(root, filename)
            if ".nii" in fullpath:
                filepath_list.append(fullpath)

    return filepath_list


def csv_to_pandas_df(csv_file):
    import pandas as pd
    data = pd.read_csv(csv_file)

    return data


def parse_raw_data_list(filepath_list, site_folder, include_sites=False, 
    subject_inclusion=None):
    
    # for script 'qap_raw_data_sublist_generator.py'

    sub_dict = {}
    inclusion_list = []
    
    # create subject inclusion list
    if subject_inclusion:
        inclusion_list = read_txt_file(subject_inclusion)    
    
    for fullpath in filepath_list:
            
        # /path_to_site_folder/subject_id/session_id/scan_id/..
        second_half = fullpath.split(site_folder)[1]
        second_half_list = second_half.split("/")
        filename = second_half_list[-1]

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
                      "be something like this:\n/data_folder/subject"\
                      "_id/session_id/scan_id/file.nii.gz\n\n" \
                      % second_half
                print err
            if subject_id not in second_half_list[0]:
                err = "\n\n[!] You are creating a participant list without " \
                      "site information included, but the site folder you " \
                      "provided is not the directory level containing the " \
                      "participants' folders - if there are multiple site " \
                      "or data folders in the site_folder directory, you " \
                      "might end up accidentally gathering data from other " \
                      "sites or collections.\n\nPlease either provide a " \
                      "different site folder input, or include sites by " \
                      "passing the --sites flag.\n\n"
                raise Exception(err)

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


def populate_custom_data_dict(data_dict, filepath, part_id, session_id, series_id,
    site_id, resource, default_series_label=None):

    # used in 'gather_custom_raw_data'

    new_dict = data_dict

    if not session_id:
        session_id = "session_1"
    if not series_id:
        if default_series_label:
            series_id = default_series_label
        else:
            series_id = "scan_1"

    if part_id not in new_dict.keys():
        new_dict[part_id] = {}
    if session_id not in new_dict[part_id].keys():
        new_dict[part_id][session_id] = {}
    if resource not in new_dict[part_id][session_id].keys():
        new_dict[part_id][session_id][resource] = {}
    if site_id not in new_dict[part_id][session_id].keys():
        new_dict[part_id][session_id]["site_name"] = site_id
    if series_id not in new_dict[part_id][session_id][resource].keys():
        new_dict[part_id][session_id][resource][series_id] = filepath    

    data_dict.update(new_dict)

    return data_dict


def gather_custom_raw_data(filepath_list, base_folder, directory_format, 
    anatomical_keywords=None, functional_keywords=None):

    # for script 'qap_flexible_sublist_generator.py'

    import os
    from qap.script_utils import populate_custom_data_dict

    if "{participant}" not in directory_format:
        pass

    data_dict = {}

    format_list = [x for x in directory_format.split("/") if x != ""]

    indices = {}
    operators = ["{site}", "{participant}", "{session}", "{series}"]

    for op in operators:
        if op in format_list:
            indices[op] = format_list.index(op)

    if anatomical_keywords:
        anatomical_keywords = [x for x in anatomical_keywords.split(" ") if x != ""]

    if functional_keywords:
        functional_keywords = [x for x in functional_keywords.split(" ") if x != ""]

    for filepath in filepath_list:

        second_half_list = [x for x in filepath.split(base_folder)[1].split("/") if x != ""]
        filename = second_half_list[-1]

        site_id = None
        part_id = None
        session_id = None
        series_id = None

        if "{site}" in indices.keys():
            site_id = second_half_list[indices["{site}"]]
        if "{participant}" in indices.keys():
            part_id = second_half_list[indices["{participant}"]]
        if "{session}" in indices.keys():
            session_id = second_half_list[indices["{session}"]]
        if "{series}" in indices.keys():
            series_id = second_half_list[indices["{series}"]]

        if anatomical_keywords:
            for word in anatomical_keywords:
                if (word in filename) or (word in session_id):

                    data_dict = populate_custom_data_dict(data_dict,
                                                   filepath, part_id,
                                                   session_id,
                                                   series_id,
                                                   site_id,
                                                   "anatomical_scan",
                                                   "anat_1")

        if functional_keywords:
            for word in functional_keywords:
                if (word in filename) or (word in session_id):

                    data_dict = populate_custom_data_dict(data_dict,
                                                   filepath, part_id,
                                                   session_id,
                                                   series_id,
                                                   site_id,
                                                   "functional_scan",
                                                   "func_1")

    if len(data_dict) == 0:
        err = "\n\n[!] No data files were found given the inputs you " \
              "provided! Double-check your setup.\n\n"
        raise Exception(err)

    return data_dict


def pull_s3_sublist(bucket_name, bucket_prefix, creds_path):

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

    # for script 'qap_aws_s3_dict_generator.py'

    s3_dict = {}

    # Read in series_list, if it is provided
    if session_list:
        try:
            session_list = os.path.abspath(session_list)
            sessions = read_txt_file(session_list)
        except:
            err = "\n\nCould not successfully read the session list.\n%s" \
                  % session_list
            raise Exception(err)

    # Read in series_list, if it is provided
    if series_list:
        try:
            series_list = os.path.abspath(series_list)
            series = read_txt_file(series_list)
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

            if ("sub-" in filename) and ("_ses-" in filename):
                err = "\n\n[!] You selected to create a data list from a " \
                      "a data directory that isn't in BIDS format, but the " \
                      "data is in BIDS format!\n\nPlease use the --BIDS " \
                      "flag for BIDS datasets.\n\n"
                raise Exception(err)

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


def create_CPAC_outputs_dict(cpac_outdir, qap_type, session_format):

    # for script 'qap_cpac_output_sublist_generator.py'

    ''' WARNING: OUT OF DATE!!! '''

    import os
    import glob

    if qap_type == "anat":

        outputs = ["anatomical_reorient", "anatomical_csf_mask", \
                   "anatomical_gm_mask", "anatomical_wm_mask", \
                   "anatomical_to_mni_linear_xfm"]

    elif qap_type == "func":

        outputs = ["mean_functional", "functional_brain_mask", \
                   "raw_functional", "coordinate_transformation"]

    outputs_dict = {}

    for sub_dir in os.listdir(cpac_outdir):

        if not os.path.isdir(os.path.join(cpac_outdir, sub_dir)):
            continue

        sessions = []

        # if the folder structure is sub_id/session_id/scan_id/...
        if session_format == 1:
            for session in os.listdir(os.path.join(cpac_outdir, sub_dir)):
                if os.path.isdir(os.path.join(cpac_outdir, sub_dir, session)):
                    sessions.append(session)

        # if there is no session label in the folder structure
        if session_format == 2:
            # since there is no session, let's assign one
            sessions = ["session_1"]

        # if the session is embedded in the subject ID
        if session_format == 3:
            subject_session = sub_dir

            if "_" not in sub_dir:
                err = "\n\n[!] You said your CPAC output directory had the " \
                      "session IDs embedded in the subject IDs, but it " \
                      "doesn't seem that way for subject ID %s!\n\nIs it " \
                      " formatted like this?   ../pipeline_output/subject_" \
                      "session/output/..\n\nIf not, you're using the wrong " \
                      "option for session_format! Use the -h flag to see " \
                      "the documentation.\n\n%s not being included in the " \
                      "subject list.\n\n" % (sub_dir, sub_dir)
                print err
                continue

            session_id = sub_dir.split("_",1)[1]
            sub_dir = sub_dir.split("_",1)[0]
            sessions = [session_id]

        for session in sessions:

            for resource in outputs:

                resource_path = ""

                if session_format == 1:
                    resource_folder = os.path.join(cpac_outdir, sub_dir, \
                                                       session, resource)
                elif session_format == 2:
                    resource_folder = os.path.join(cpac_outdir, sub_dir, \
                                                       resource)

                elif session_format == 3:
                    resource_folder = os.path.join(cpac_outdir, \
                                                       subject_session, \
                                                       resource)

                # if this current resource/output does not exist for this
                # subject, go to next resource in list
                if not os.path.isdir(resource_folder):
                    continue

                if qap_type == "anat":

                    ''' until CPAC writes multiple anat scans in the '''
                    ''' output folder structure '''
                    scans = ["anat_1"]

                if qap_type == "func":
    
                    scans = []

                    for item in os.listdir(resource_folder):
                        if os.path.isdir(os.path.join(resource_folder, item)):
                            item = item.replace("_scan_","")
                            item = item.replace("_rest","")
                            scans.append(item)

                for scan in scans:

                    if qap_type == "anat":

                        if "mask" in resource:
                            resource_paths = glob.glob(os.path.join(resource_folder, "*", "*"))
                        else:
                            resource_paths = glob.glob(os.path.join(resource_folder, "*"))

                        if len(resource_paths) == 1:
                            resource_path = resource_paths[0]
                        else:
                            print "\nMultiple files for %s for subject %s!!" \
                                  % (resource, sub_dir)
                            print "Check the directory: %s" \
                                      % resource_folder
                            print "%s for %s has not been included in the " \
                                  "subject list.\n" % (resource, sub_dir)
                            continue

                    if qap_type == "func":

                        fullscan = "_scan_" + scan + "_rest"

                        resource_paths = glob.glob(os.path.join(resource_folder, fullscan, "*"))

                        if len(resource_paths) == 1:
                            resource_path = resource_paths[0]
                        else:
                            print "\nMultiple files for %s for subject %s!!" \
                                  % (resource, sub_dir)
                            print "Check the directory: %s" \
                                      % resource_folder
                            print "%s for %s has not been included in the " \
                                  "subject list.\n" % (resource, sub_dir)
                            continue

                    ''' put a catch here for multiple files '''

                    if sub_dir not in outputs_dict.keys():
                        outputs_dict[sub_dir] = {}

                    if session not in outputs_dict[sub_dir].keys():
                        outputs_dict[sub_dir][session] = {}

                    if resource not in outputs_dict[sub_dir][session].keys():
                        outputs_dict[sub_dir][session][resource] = {}

                    if scan not in outputs_dict[sub_dir][session][resource].keys():
                        outputs_dict[sub_dir][session][resource][scan] = resource_path

    # make up for QAP - CPAC resource naming discrepancy
    for subid in outputs_dict.keys():

        for session in outputs_dict[subid].keys():

            for resource in outputs_dict[subid][session].keys():

                if resource == "motion_correct":

                    filepath = outputs_dict[subid][session]["motion_correct"]

                    outputs_dict[subid][session]["func_motion_correct"] = \
                        filepath

                    del outputs_dict[subid][session]["motion_correct"]

                if resource == "anatomical_to_mni_linear_xfm":

                    filepath = outputs_dict[subid][session]["anatomical_to_mni_linear_xfm"]

                    outputs_dict[subid][session]["flirt_affine_xfm"] = \
                        filepath

                    del outputs_dict[subid][session]["anatomical_to_mni_linear_xfm"]

    return outputs_dict


def qap_csv_correlations(data_old, data_new, replacements=None):

    # for script 'qap_test_correlations.py'

    import numpy as np
    import pandas as pd
    import scipy.stats

    metric_list = ["EFC","SNR","FBER","CNR","FWHM","Qi1","Cortical Contrast",
        "Ghost_x", "Ghost_y", "Ghost_z", "GCOR", "RMSD (Mean)", 
        "Quality (Mean)", "Fraction of Outliers (Mean)", "Std. DVARS (Mean)", 
        "Fraction of OOB Outliers (Mean)"]

    # update datasets if necessary
    if replacements:
        replace_dict = {}
        for word_couple in replacements:
            if "," not in word_couple:
                err = "\n\n[!] In the replacements text file, the old " \
                      "substring and its replacement must be separated " \
                      "by a comma.\n\nLine: %s\n\n" % word_couple
                raise Exception(err)
            word = word_couple.split(",")[0]
            new = word_couple.split(",")[1]
            replace_dict[word] = new

        data_old.rename(columns=replace_dict, inplace=True)
        data_new.rename(columns=replace_dict, inplace=True)

    # remove nulls
    data_old = data_old[pd.notnull(data_old["Participant"])]
    data_new = data_new[pd.notnull(data_new["Participant"])]

    for metric in metric_list:
        if metric in data_old:
            data_old = data_old[pd.notnull(data_old[metric])]
        if metric in data_new:
            data_new = data_new[pd.notnull(data_new[metric])]

    # make sure participant IDs are strings (if they are all digits, can be
    # mistakenly read in as ints or floats)
    if data_old["Participant"].dtype != str:
        data_old["Participant"] = data_old["Participant"].astype(int).astype(str)
    #if data_new[partic_label_new].dtype != str:
    #    data_new[partic_label_new] = data_new[partic_label_new].astype(int).astype(str)

    # make sure both DFs match
    data_merged = pd.merge(data_old, data_new, 
        on=["Participant","Session","Series"], how="inner", 
        suffixes=("_OLD","_NEW"))

    # correlate the numbers!
    correlations_dict = {}
    for metric in metric_list:
        metric_old = metric + "_OLD"
        metric_new = metric + "_NEW"
        if (metric_old in data_merged) and (metric_new in data_merged):
            metric_old_val = data_merged[metric_old]#.flatten()
            metric_new_val = data_merged[metric_new]#.flatten()
            correlations_dict[metric] = scipy.stats.pearsonr(metric_old_val, metric_new_val)

    return correlations_dict


def write_inputs_dict_to_yaml_file(input_dict, yaml_outpath):

    import os
    import yaml         

    yaml_outpath = os.path.abspath(yaml_outpath)
    if (".yml" not in yaml_outpath) and (".yaml" not in yaml_outpath):
        yaml_outpath += ".yml"

    # write yaml file
    try:
        with open(yaml_outpath,"wt") as f:
            f.write(yaml.dump(input_dict, default_flow_style=True))
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