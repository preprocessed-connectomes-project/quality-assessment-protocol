
def gather_raw_data(site_folder, yaml_outpath, scan_type, \
                        include_sites=False, subject_inclusion=None):

    import os
    import yaml
    
    sub_dict = {}
    inclusion_list = []
    

    # create subject inclusion list
    if subject_inclusion != None:
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
                              "structure. Is it in the correct format?\n\n" \
                              "Your directory structure:\n%s\n\nIt should " \
                              "be something like this:\n/site_folder/subject"\
                              "_id/session_id/scan_id/file.nii.gz\n\n" \
                              % second_half
                        raise Exception(err)

                else:

                    try:

                        subject_id = second_half_list[-4]
                        session_id = second_half_list[-3]
                        scan_id = second_half_list[-2]

                    except:
                        
                        err = "\n\n[!] Could not parse the data directory " \
                              "structure. Is it in the correct format?\n\n" \
                              "Your directory structure:\n%s\n\nIt should " \
                              "be something like this:\n/subject_id" \
                              "/session_id/scan_id/file.nii.gz\n\n" \
                              % second_half
                        raise Exception(err)


                if subject_inclusion == None:
                    inclusion_list.append(subject_id)
                

                # assign default scan names if the file structure of the data
                # does not have a directory level for scan
                if ".nii" in scan_id:
                    if scan_type == "anat":
                        scan_id = "anat_1"
                    elif scan_type == "func":
                        scan_id = "rest_1"
                
                #sub_info = (subject_id, session_id, scan_id)
                
                if ("anat" in scan_id) or ("anat" in filename) or \
                    ("mprage" in filename):
                    
                    resource = "anatomical_scan"
                    
                if ("rest" in scan_id) or ("rest" in filename) or \
                    ("func" in scan_id) or ("func" in filename):
                    
                    resource = "functional_scan"

                
                if (scan_type == "anat") and \
                    (resource == "anatomical_scan") and \
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
                    
                        
                if (scan_type == "func") and \
                     (resource == "functional_scan") and \
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
 
 
                
    with open(yaml_outpath,"wt") as f:
        f.write(yaml.dump(sub_dict))



def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("site_folder", type=str, \
                            help="full path to the directory holding the " \
                                 "raw data, organized by site/subject/" \
                                 "session/scan/file.nii.gz")

    parser.add_argument("outfile_path", type=str, \
                            help="full path for the generated subject list")
                            
    parser.add_argument("scan_type", type=str, \
                            help="'anat' or 'func', depending on which QAP " \
                                 "measures you will be using the subject " \
                                 "list for")
                                 
    parser.add_argument('--sites', action='store_true', \
                            help="include this flag if you wish to include " \
                                 "site information in your subject list - " \
                                 "data must be organized as /site_name/" \
                                 "subject_id/session_id/scan_id/..")

    parser.add_argument("--include", type=str, \
                            help="text file containing subject IDs of the " \
                                 "subjects you want to include - leave " \
                                 "this out if you want to run all of them")

    args = parser.parse_args()

    # run it!
    gather_raw_data(args.site_folder, args.outfile_path, args.scan_type, \
                        args.sites, args.include)



if __name__ == "__main__":
    main()
    
    
