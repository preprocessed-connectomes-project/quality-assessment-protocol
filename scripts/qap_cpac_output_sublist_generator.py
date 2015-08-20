

def run(cpac_outdir, outfile_name, qap_type, session_format):

    import os
    import glob
    import yaml

 
    if qap_type == "anat":

        outputs = ["anatomical_reorient", "anatomical_csf_mask", \
                   "anatomical_gm_mask", "anatomical_wm_mask", \
                   "anatomical_to_mni_linear_xfm"]

    elif qap_type == "func":

        outputs = ["mean_functional", "functional_brain_mask", \
                   "motion_correct", "coordinate_transformation"]


    outputs_dict = {}


    for sub_dir in os.listdir(cpac_outdir):

        if not os.path.isdir(os.path.join(cpac_outdir, sub_dir)):
            continue

        sessions = []

        # if the folder structure is sub_id/session_id/scan_id/...
        if session_format == 1:
            for session in os.listdir(os.path.join(cpac_outdir, sub_dir)):
                if os.path.isdir(os.path.join(cpac_outdir, sub_dir, session)):
                    sessions.append(session)

        # if there is no session label in the folder structure
        if session_format == 2:
            # since there is no session, let's assign one
            sessions = ["session_1"]

        # if the session is embedded in the subject ID
        if session_format == 3:
            subject_session = sub_dir

            if "_" not in sub_dir:
                err = "\n\n[!] You said your CPAC output directory had the " \
                      "session IDs embedded in the subject IDs, but it " \
                      "doesn't seem that way for subject ID %s!\n\nIs it " \
                      " formatted like this?   ../pipeline_output/subject_" \
                      "session/output/..\n\nIf not, you're using the wrong " \
                      "option for session_format! Use the -h flag to see " \
                      "the documentation.\n\n%s not being included in the " \
                      "subject list.\n\n" % (sub_dir, sub_dir)
                print err
                continue

            session_id = sub_dir.split("_",1)[1]
            sub_dir = sub_dir.split("_",1)[0]
            sessions = [session_id]


        for session in sessions:

            for resource in outputs:

                resource_path = ""

                if session_format == 1:
                    resource_folder = os.path.join(cpac_outdir, sub_dir, \
                                                       session, resource)
                elif session_format == 2:
                    resource_folder = os.path.join(cpac_outdir, sub_dir, \
                                                       resource)

                elif session_format == 3:
                    resource_folder = os.path.join(cpac_outdir, \
                                                       subject_session, \
                                                       resource)

                # if this current resource/output does not exist for this
                # subject, go to next resource in list
                if not os.path.isdir(resource_folder):
                    continue


                if qap_type == "anat":

                    ''' until CPAC writes multiple anat scans in the '''
                    ''' output folder structure '''
                    scans = ["anat_1"]


                if qap_type == "func":
    
                    scans = []

                    for item in os.listdir(resource_folder):
                        if os.path.isdir(os.path.join(resource_folder, item)):
                            item = item.replace("_scan_","")
                            item = item.replace("_rest","")
                            scans.append(item)


                for scan in scans:

                    if qap_type == "anat":

                        if "mask" in resource:
                            resource_paths = glob.glob(os.path.join(resource_folder, "*", "*"))
                        else:
                            resource_paths = glob.glob(os.path.join(resource_folder, "*"))

                        if len(resource_paths) == 1:
                            resource_path = resource_paths[0]
                        else:
                            print "\nMultiple files for %s for subject %s!!" \
                                  % (resource, sub_dir)
                            print "Check the directory: %s" \
                                      % resource_folder
                            print "%s for %s has not been included in the " \
                                  "subject list.\n" % (resource, sub_dir)
                            continue

                    if qap_type == "func":

                        fullscan = "_scan_" + scan + "_rest"

                        resource_paths = glob.glob(os.path.join(resource_folder, fullscan, "*"))

                        if len(resource_paths) == 1:
                            resource_path = resource_paths[0]
                        else:
                            print "\nMultiple files for %s for subject %s!!" \
                                  % (resource, sub_dir)
                            print "Check the directory: %s" \
                                      % resource_folder
                            print "%s for %s has not been included in the " \
                                  "subject list.\n" % (resource, sub_dir)
                            continue


                    ''' put a catch here for multiple files '''


                    if sub_dir not in outputs_dict.keys():
                        outputs_dict[sub_dir] = {}

                    if session not in outputs_dict[sub_dir].keys():
                        outputs_dict[sub_dir][session] = {}

                    if resource not in outputs_dict[sub_dir][session].keys():
                        outputs_dict[sub_dir][session][resource] = {}

                    if scan not in outputs_dict[sub_dir][session][resource].keys():
                        outputs_dict[sub_dir][session][resource][scan] = resource_path




    # make up for QAP - CPAC resource naming discrepancy
    for subid in outputs_dict.keys():

        for session in outputs_dict[subid].keys():

            for resource in outputs_dict[subid][session].keys():

                if resource == "motion_correct":

                    filepath = outputs_dict[subid][session]["motion_correct"]

                    outputs_dict[subid][session]["func_motion_correct"] = \
                        filepath

                    del outputs_dict[subid][session]["motion_correct"]

                if resource == "anatomical_to_mni_linear_xfm":

                    filepath = outputs_dict[subid][session]["anatomical_to_mni_linear_xfm"]

                    outputs_dict[subid][session]["flirt_affine_xfm"] = \
                        filepath

                    del outputs_dict[subid][session]["anatomical_to_mni_linear_xfm"]



    outfile = os.path.join(os.getcwd(), outfile_name + ".yml")

    with open(outfile, 'w') as f:

        f.write(yaml.dump(outputs_dict, default_flow_style=True))



def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("cpac_output_dir", type=str, \
                            help="path to output directory of a CPAC " \
                                 "individual-analysis run - ex. the folder " \
                                 "named after the pipeline, containing the " \
                                 "subject folders")

    parser.add_argument("outfile_name", type=str, \
                            help="name for the generated subject list")

    parser.add_argument("qap_type", type=str, \
                            help="'anat' or 'func' - whether to extract " \
                                 "the anatomical or functional outputs, " \
                                 "depending on which type of QAP measures " \
                                 "you wish to run ('anat' for anatomical " \
                                 "spatial, 'func' for functional spatial " \
                                 "or temporal)")

    parser.add_argument("session_format", type=int, \
                            help="input as integer: '1' if your CPAC output "\
                                 "file structure is organized as /subject_id"\
                                 "/session_id/output/.., '2' if there are " \
                                 "no sessions and the file structure is " \
                                 "organized as /subject_id/output/.., '3' " \
                                 "if the session ID is embedded in the " \
                                 "subject ID like so: /subject_session/" \
                                 "output/..")

    args = parser.parse_args()

    # run it!
    run(args.cpac_output_dir, args.outfile_name, args.qap_type, \
            args.session_format)



if __name__ == "__main__":
    main()



