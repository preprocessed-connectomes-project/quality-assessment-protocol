#!/usr/bin/env python

def gather_bids_data(dataset_folder, yaml_outpath, subject_inclusion=None,
                     scan_type=None):

    import os
    import os.path as op
    import yaml
    from glob import glob

    sub_dict = {}
    inclusion_list = []

    subject_ids = [x for x in next(os.walk(dataset_folder))[1] if x.startswith("sub-")]

    if scan_type is None:
        scan_type = 'functional anatomical'

    get_anat = 'anatomical' in scan_type
    get_func = 'functional' in scan_type

    if not subject_ids:
        raise Exception("This does not appear to be a BIDS dataset.")

    # create subject inclusion list
    if subject_inclusion is not None:
        with open(subject_inclusion, "r") as f:
            inclusion_list = f.readlines()
        # remove any /n's
        inclusion_list = map(lambda s: s.strip(), inclusion_list)

        subject_ids = set(subject_ids).intersection(inclusion_list)

    for subject_id in subject_ids:
        # TODO: implement multisession support
        session_name = 'session_1'

        anatomical_scans = sorted(glob(op.join(
            dataset_folder, subject_id, "anat",
            "%s_*T1w.nii.gz" % subject_id, )))

        functional_scans = sorted(glob(op.join(
            dataset_folder, subject_id, "func",
            "%s_*bold.nii.gz" % subject_id, )))

        if anatomical_scans or functional_scans:
            sub_dict[subject_id] = {session_name: {}}
            if anatomical_scans and get_anat:
                sub_dict[subject_id][session_name]["anatomical_scan"] = {}
                for i, anatomical_scan in enumerate(anatomical_scans):
                    sub_dict[subject_id][session_name]["anatomical_scan"]["anat_%d"%(i+1)] = op.abspath(anatomical_scan)
            if functional_scans and get_func:
                sub_dict[subject_id][session_name]["functional_scan"] = {}
                for i, anatomical_scan in enumerate(functional_scans):
                    sub_dict[subject_id][session_name]["functional_scan"]["func_%d"%(i+1)] = op.abspath(anatomical_scan)

    with open(yaml_outpath, "wt") as f:
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

    parser.add_argument("--include", type=str,
                        help="text file containing subject IDs of the "
                             "subjects you want to include - leave "
                             "this out if you want to run all of them")

    parser.add_argument(
        '-t', '--type', type=str, choices=['functional', 'anatomical'],
        help='include only functional or anatomical')

    args = parser.parse_args()

    # run it!
    gather_bids_data(args.dataset_folder, args.outfile_path, args.include,
                     args.type)


if __name__ == "__main__":
    main()
