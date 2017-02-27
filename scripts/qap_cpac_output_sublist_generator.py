#!/usr/bin/env python

def main():

    from qap.script_utils import create_CPAC_outputs_dict, \
                                 write_inputs_dict_to_yaml_file

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
    outputs_dict = create_CPAC_outputs_dict(args.cpac_output_dir, 
        args.qap_type, args.session_format)
    write_inputs_dict_to_yaml_file(outputs_dict, args.outfile_name)


if __name__ == "__main__":
    main()



