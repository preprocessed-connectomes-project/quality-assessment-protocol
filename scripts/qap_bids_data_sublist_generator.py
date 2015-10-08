
def gather_bids_data(dataset_folder, yaml_outpath, subject_inclusion=None):

    import os
    import yaml
    from glob import glob
    
    sub_dict = {}
    inclusion_list = []
    
    subject_ids = [x[4:] for x in next(os.walk(dataset_folder))[1] if x.startswith("sub-")]
    
    if not subject_ids:
        raise Exception("This does not appear to be a BIDS dataset.")

    # create subject inclusion list
    if subject_inclusion != None:
        with open(subject_inclusion, "r") as f:
            inclusion_list = f.readlines()
        # remove any /n's
        inclusion_list = map(lambda s: s.strip(), inclusion_list)
        
        subject_ids = set(subject_ids).intersection(inclusion_list)
    
    
    for subject_id in subject_ids:
        #TODO: implement multisession support
        session_name = 'session_1'
        
        anatomical_scans = sorted(glob(os.path.join(dataset_folder, 
                                             "sub-%s"%subject_id, 
                                             "anat", 
                                             "sub-%s_*T1w.nii.gz"%subject_id,
                                             )))
                                             
        functional_scans = sorted(glob(os.path.join(dataset_folder, 
                                             "sub-%s"%subject_id, 
                                             "func", 
                                             "sub-%s_*bold.nii.gz"%subject_id,
                                             )))
        
        if anatomical_scans or functional_scans:
            sub_dict[subject_id] = {session_name:{}}
            if anatomical_scans:
                sub_dict[subject_id][session_name]["anatomical_scan"] = {}
                for i, anatomical_scan in enumerate(anatomical_scans):
                    sub_dict[subject_id][session_name]["anatomical_scan"]["anat_%d"%(i+1)] = anatomical_scan
            if functional_scans:
                sub_dict[subject_id][session_name]["functional_scan"] = {}
                for i, anatomical_scan in enumerate(functional_scans):
                    sub_dict[subject_id][session_name]["functional_scan"]["func_%d"%(i+1)] = anatomical_scan
    
    with open(yaml_outpath,"wt") as f:
        f.write(yaml.dump(sub_dict))

def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("dataset_folder", type=str, \
                            help="full path to the directory holding the " \
                                 "raw data, organized according to the BIDS specification" \
                                 "http://bids.neuroimaging.io")

    parser.add_argument("outfile_path", type=str, \
                            help="full path for the generated subject list")
                                 

    parser.add_argument("--include", type=str, \
                            help="text file containing subject IDs of the " \
                                 "subjects you want to include - leave " \
                                 "this out if you want to run all of them")

    args = parser.parse_args()

    # run it!
    gather_bids_data(args.dataset_folder, args.outfile_path, args.include)



if __name__ == "__main__":
    main()
    
    
