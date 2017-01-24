#!/usr/bin/env python

def main():

    import argparse
    from qap.script_utils import pull_s3_sublist, \
                                 read_txt_file, \
                                 parse_raw_data_list, \
                                 write_inputs_dict_to_yaml_file
    from qap.bids_utils import extract_bids_data

    parser = argparse.ArgumentParser()
 
    parser.add_argument("bucket_name", type=str, \
                            help="the name of your AWS S3 bucket")
 
    parser.add_argument("bucket_prefix", type=str, \
                            help="the filepath prefix to the top level of " \
                                 "your raw data directory on S3 storage")
 
    parser.add_argument("outfile_path", type=str, \
                            help="the full filepath for the S3 subject " \
                                 "YAML dictionary this script will create")

    parser.add_argument("--creds_path", type=str, \
                            help="the path to the file containing your AWS " \
                                 "credentials")

    parser.add_argument("--participant_list", type=str, \
                            help="filepath to a text file containing the " \
                                 "names of participants you want included, " \
                                 "one on each line")

    parser.add_argument("--BIDS", action="store_true", \
                            help="if the dataset is in BIDS format")
 
    args = parser.parse_args()

    # run it!
    s3_list = pull_s3_sublist(args.bucket_name, args.bucket_prefix,
                              args.creds_path)

    # create subject inclusion list
    if args.participant_list:
        inclusion_list = read_txt_file(args.participant_list)
    else:
        inclusion_list = None

    if not args.BIDS:
        s3_dict = parse_raw_data_list(s3_list, args.bucket_prefix,
                                      inclusion_list=inclusion_list,
                                      s3_bucket=True)
    else:
        s3_dict = extract_bids_data(s3_list, inclusion_list)

    write_inputs_dict_to_yaml_file(s3_dict, args.outfile_path)


if __name__ == "__main__":
    main()
