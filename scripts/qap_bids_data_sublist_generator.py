#!/usr/bin/env python

def gather_bids_data(dataset_folder, subject_inclusion=None,
                     scan_type=None):

    import os
    import os.path as op
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
        sessions = [x for x in next(os.walk(op.join(dataset_folder,subject_id)))[1] if x.startswith("ses-")]

        if not sessions:
            sessions = ['']

        for session_dir in sessions:

            # the session name is just the session directory name
            session_name = session_dir

            # unless the session directory is missing, in which
            # case set the session_name to ses-1
            if not session_name:
                session_name = "ses-1"

            # for now restrict our analysis to T1 weighted anats and
            # BOLD weighted functionals
            anatomical_scans = sorted(glob(op.join(
                dataset_folder,  subject_id, session_dir, "anat",
                "%s_*T1w.nii.gz" % subject_id, )))

            functional_scans = sorted(glob(op.join(
                dataset_folder, subject_id, session_dir, "func",
                "%s_*bold.nii.gz" % subject_id, )))

            if anatomical_scans or functional_scans:
                sub_dict[subject_id] = {session_name: {}}
                if anatomical_scans and get_anat:
                    sub_dict[subject_id][session_name]["anatomical_scan"] = {}
                    for i, anatomical_scan in enumerate(anatomical_scans):
                        sd={s.split("-")[0]:s.split("-")[1] \
                            for s in anatomical_scan.split("_") \
                            if len(s.split("-")) > 1}
                        anat_key=(anatomical_scan.split("_")[-1]).split(".")[0]
                        if "acq" in sd.keys():
                            anat_key = "_".join([anat_key,sd["acq"]])
                        sub_dict[subject_id][session_name]["anatomical_scan"][anat_key] = op.abspath(anatomical_scan)

                if functional_scans and get_func:
                    sub_dict[subject_id][session_name]["functional_scan"] = {}
                    for i, functional_scan in enumerate(functional_scans):
                        sd={s.split("-")[0]:s.split("-")[1] \
                            for s in functional_scan.split("_") \
                            if len(s.split("-")) > 1}
                        func_key=(functional_scan.split("_")[-1]).split(".")[0]
                        if "task" in sd.keys():
                            func_key="_".join([func_key,sd["task"]])
                        else:
                            raise Exception("%s: Missing task key , this does"\
                                "not appear to be a BIDS dataset."\
                                %(functional_scan))
                        if "acq" in sd.keys():
                            func_key = "_".join([func_key,sd["acq"]])
                        sub_dict[subject_id][session_name]["functional_scan"][func_key] = op.abspath(anatomical_scan)

    return(sub_dict)



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

    # run it!
    sub_dict = gather_bids_data(args.dataset_folder, args.include,
                     args.type)

    with open(args.outfile_path, "wt") as f:
        f.write(yaml.dump(sub_dict))

if __name__ == "__main__":
    main()
