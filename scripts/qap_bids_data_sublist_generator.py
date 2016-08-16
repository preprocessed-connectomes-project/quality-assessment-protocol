#!/usr/bin/env python
from qap.bids_utils import gather_nifti_file_paths, extract_bids_data


def main():

    import argparse
    import yaml

    parser = argparse.ArgumentParser()

    parser.add_argument("dataset_folder", type=str, \
                        help="full path to the directory holding the " \
                             "raw data, organized according to the BIDS specification" \
                             "http://bids.neuroimaging.io")

    parser.add_argument("outfile_path", type=str, \
                        help="full path for the generated subject list")

    parser.add_argument("--include", type=str,
                        help="text file containing subject IDs of the "
                             "subjects you want to include - leave "
                             "this out if you want to run all of them")

    parser.add_argument(
        '-t', '--type', type=str, choices=['functional', 'anatomical'],
        help='include only functional or anatomical')

    args = parser.parse_args()

    # create subject inclusion list
    inclusion_list=[]
    if args.include is not None:
        with open(args.include, "r") as f:
            inclusion_list = f.readlines()
        # remove any /n's
        inclusion_list = map(lambda s: s.strip(), inclusion_list)

    file_path_list = gather_nifti_file_paths(args.dataset_folder)

    # run it!
    sub_dict = extract_bids_data(file_path_list, inclusion_list)

    with open(args.outfile_path, "wt") as f:
        f.write(yaml.dump(sub_dict))

if __name__ == "__main__":
    main()
