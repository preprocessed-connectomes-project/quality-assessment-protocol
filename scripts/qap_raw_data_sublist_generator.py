#!/usr/bin/env python

def main():

    import argparse

    from qap.script_utils import gather_raw_data, \
                                 write_inputs_dict_to_yaml_file

    parser = argparse.ArgumentParser()

    parser.add_argument("site_folder", type=str, \
                            help="full path to the directory holding the " \
                                 "raw data, organized by site/subject/" \
                                 "session/scan/file.nii.gz")

    parser.add_argument("outfile_path", type=str, \
                            help="full path for the generated subject list")
                                 
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
    sub_dict = gather_raw_data(args.site_folder, args.sites, args.include)

    write_inputs_dict_to_yaml_file(sub_dict, args.outfile_path)


if __name__ == "__main__":
    main()
    
    
