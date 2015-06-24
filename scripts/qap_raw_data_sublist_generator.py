
def gather_raw_data(site_folder, yaml_outpath, scan_type, subject_inclusion):

    import os
    import yaml
    
    sub_dict = {}
    
    # create subject inclusion list
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

                subject_id = second_half_list[0]
                
                session_id = second_half_list[1]
                
                scan_id = second_half_list[2]
                
                
                #sub_info = (subject_id, session_id, scan_id)
                
                if ("anat" in scan_id) or ("anat" in filename) or \
                    ("mprage" in filename):
                    
                    resource = "anatomical_scan"
                    
                if ("rest" in scan_id) or ("rest" in filename) or \
                    ("func" in scan_id) or ("func" in filename):
                    
                    resource = "functional_scan"

                
                if (scan_type == "anat") and \
                    (resource == "anatomical_scan") and \
                        (subject_id is in inclusion_list):

                    if subject_id not in sub_dict.keys():
                        sub_dict[subject_id] = {}
                    
                    if session_id not in sub_dict[subject_id].keys():
                        sub_dict[subject_id][session_id] = {}
                    
                    if resource not in sub_dict[subject_id][session_id].keys():
                        sub_dict[subject_id][session_id][resource] = {}
                                                       
                    if scan_id not in sub_dict[subject_id][session_id][resource].keys():
                        sub_dict[subject_id][session_id][resource][scan_id] = fullpath
                    
                        
                if (scan_type == "func") and \
                     (resource == "functional_scan") and \
                         (subject_id is in inclusion_list):
                
                    if subject_id not in sub_dict.keys():
                        sub_dict[subject_id] = {}
                    
                    if session_id not in sub_dict[subject_id].keys():
                        sub_dict[subject_id][session_id] = {}
                    
                    if resource not in sub_dict[subject_id][session_id].keys():
                        sub_dict[subject_id][session_id][resource] = {}
                                    
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
                                 
    parser.add_argument("--include", type=str, \
                            help="text file containing subject IDs of the " \
                                 "subjects you want to include - leave " \
                                 "this out if you want to run all of them")

    args = parser.parse_args()

    # run it!
    gather_raw_data(args.site_folder, args.outfile_path, args.scan_type, \
                        args.include)



if __name__ == "__main__":
    main()
