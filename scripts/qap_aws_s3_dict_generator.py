#!/usr/bin/env python

def main():

    import argparse
    from qap.script_utils import pull_s3_sublist, \
                                 create_subdict_from_s3_list, \
                                 write_inputs_dict_to_yaml_file

    parser = argparse.ArgumentParser()
 
    parser.add_argument("bucket_name", type=str, \
                            help="the name of your AWS S3 bucket")
 
    parser.add_argument("bucket_prefix", type=str, \
                            help="the filepath prefix to the top level of " \
                                 "your raw data directory on S3 storage")
 
    parser.add_argument("creds_path", type=str, \
                            help="the path to the file containing your AWS " \
                                 "credentials")

    parser.add_argument("outfile_path", type=str, \
                            help="the full filepath for the S3 subject " \
                                 "YAML dictionary this script will create")

    parser.add_argument('--include_sites', action='store_true', \
                            help="include this flag if you wish to include " \
                                 "site information in your subject list - " \
                                 "data must be organized as /site_name/" \
                                 "subject_id/session_id/scan_id/..")

    parser.add_argument("--session_list", type=str, \
                            help="filepath to a text file containing the " \
                                 "names of sessions you want included, one " \
                                 "on each line")

    parser.add_argument("--series_list", type=str, \
                            help="filepath to a text file containing the " \
                                 "names of series you want included, one " \
                                 "on each line")

    parser.add_argument("--BIDS", action="store_true", \
                            help="if the dataset is in BIDS format")
 
    args = parser.parse_args()

    # run it!
    s3_list = pull_s3_sublist(args.bucket_name, args.bucket_prefix, \
        args.creds_path)

    s3_dict = create_subdict_from_s3_list(s3_list, args.bucket_prefix, \
        args.session_list, args.series_list, args.BIDS)

    write_inputs_dict_to_yaml_file(s3_dict, args.outfile_path)


if __name__ == "__main__":
    main()
    
    
