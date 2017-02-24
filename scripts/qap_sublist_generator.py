#!/usr/bin/env python


def main():

    import os
    import argparse
    from qap.script_utils import gather_filepath_list, \
                                 read_txt_file, \
                                 pull_s3_sublist, \
                                 parse_raw_data_list, \
                                 write_inputs_dict_to_yaml_file

    parser = argparse.ArgumentParser()

    parser.add_argument("data_folder", type=str,
                            help="full path to the directory holding the "
                                 "raw data, organized by /site/subject/"
                                 "session/scan/file.nii.gz or in BIDS format")

    parser.add_argument("outfile_path", type=str,
                            help="filename for the generated subject list")

    parser.add_argument("--include", type=str,
                            help="text file containing participant IDs of "
                                 "the subjects you want to include - leave "
                                 "this out if you want to run all of them")

    parser.add_argument("--BIDS", action="store_true",
                            help="if the dataset is in BIDS format")

    parser.add_argument("--creds_path", type=str,
                            help="the path to the file containing your AWS "
                                 "credentials")

    args = parser.parse_args()

    # run it!
    if "s3://" in args.data_folder:
        data_dir = args.data_folder
        filepath_list = pull_s3_sublist(data_dir, args.creds_path)
    else:
        data_dir = os.path.abspath(args.data_folder)
        filepath_list = gather_filepath_list(data_dir)

    if args.include:
        inclusion_list = read_txt_file(args.include)
    else:
        inclusion_list = None

    if not args.BIDS:
        sub_dict = parse_raw_data_list(filepath_list, data_dir,
                                       inclusion_list=inclusion_list)
    else:
        from qap.bids_utils import bids_gen_qap_sublist
        sub_dict = bids_gen_qap_sublist(data_dir, filepath_list)

    write_inputs_dict_to_yaml_file(sub_dict, args.outfile_path)


if __name__ == "__main__":
    main()

