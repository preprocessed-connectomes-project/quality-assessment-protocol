

def run(cpac_outdir, outfile_name, qap_type):

    import os
    import glob
    import yaml

    ''' SCRIPT UNDER CONSTRUCTION!!! '''

    if qap_type == "anat":

        outputs = ["anatomical_reorient", "anatomical_csf_mask", \
                   "anatomical_gm_mask", "anatomical_wm_mask", \
                   "ants_affine_xfm"]

    elif qap_type == "func":

        outputs = ["mean_functional", "functional_brain_mask", \
                   "motion_correct", "coordinate_transformation"]


    outputs_dict = {}


    for sub_dir in os.listdir(cpac_outdir):

        for resource in outputs:

            resource_path = ""

            filepaths = glob.glob(os.path.join(cpac_outdir, sub_dir, \
                                               resource))


            if len(filepaths) == 1:

                for root, folders, files in os.walk(filepaths[0]):
                
                    for filename in files:

                        if qap_type == "func":

                            '''scan_id was removed! now must support multiple scans, multiple strategies, etc.!'''
                            if scan_id in root:

                                filepath = os.path.join(root, filename)

                                resource_path = filepath

                        elif qap_type == "anat":

                            filepath = os.path.join(root, filename)

                            resource_path = filepath


                if sub_dir not in outputs_dict.keys():
                    outputs_dict[sub_dir] = {}

                if resource not in outputs_dict[sub_dir].keys():
                    outputs_dict[sub_dir][resource] = resource_path


    # make up for QAP - CPAC resource naming discrepancy
    for subid in outputs_dict.keys():

        for resource in outputs_dict[subid].keys():

            if resource == "motion_correct":

                filepath = outputs_dict[subid]["motion_correct"]

                outputs_dict[subid]["func_motion_correct"] = filepath



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
                                 "you wish to run ('anat' for Spatial, " \
                                 "'func' for Spatial EPI or Temporal)")

    parser.add_argument("-s", "--scan_id", type=str, \
                            help="(for qap_type 'func' only) - the name of " \
                                 "the scan you wish to run the measures for")

    args = parser.parse_args()

    # run it!
    run(args.cpac_output_dir, args.outfile_name, args.qap_type, args.scan_id)



if __name__ == "__main__":
    main()



