#!/usr/bin/env python

def main():

    import argparse

    from qap.script_utils import gather_filepath_list, \
                                 read_txt_file, \
                                 pull_s3_sublist, \
                                 parse_raw_data_list, \
                                 write_inputs_dict_to_yaml_file

    parser = argparse.ArgumentParser()

    parser.add_argument("data_folder", type=str,
                            help="full path to the directory holding the " \
                                 "raw data, organized by /site/subject/" \
                                 "session/scan/file.nii.gz or in BIDS " \
                                 "format")

    parser.add_argument("outfile_path", type=str,
                            help="filename for the generated subject list")

    parser.add_argument("--include", type=str,
                            help="text file containing participant IDs of " \
                                 "the subjects you want to include - leave " \
                                 "this out if you want to run all of them")

    parser.add_argument("--BIDS", action="store_true",
                            help="if the dataset is in BIDS format")

    parser.add_argument("--creds_path", type=str,
                            help="the path to the file containing your AWS " \
                                 "credentials")

    args = parser.parse_args()

    s3 = False

    # run it!
    if "s3://" in args.data_folder:
        s3 = True
        filepath_list, data_dir = \
            pull_s3_sublist(args.data_folder, args.creds_path)
    else:
        data_dir = args.data_folder
        filepath_list = gather_filepath_list(args.data_folder)

    if args.include:
        inclusion_list = read_txt_file(args.include)
    else:
        inclusion_list = None

    if not args.BIDS:
        sub_dict = parse_raw_data_list(filepath_list, data_dir,
                                       inclusion_list=inclusion_list,
                                       s3_bucket=s3)
    else:
        from qap.bids_utils import extract_bids_data
        sub_dict = extract_bids_data(filepath_list, inclusion_list)

    write_inputs_dict_to_yaml_file(sub_dict, args.outfile_path)


if __name__ == "__main__":
    main()

