
# this script takes in all of the generated QAP numbers from the QAP pipeline
# output directory you provide it and aggregates them into one CSV file

# the QAP pipelines already create their own combined CSV file, this script is
# for when you need to do it manually

# script usage:
#     > qap_merge_outputs.py <qap pipeline output directory>

def csv_merge(csv_list, outfile):

    import csv

    fields = []

    sub_nums_list = []

    for sub in csv_list:

        with open(sub, "r") as in_f:

            csv_reader = csv.DictReader(in_f)

            if not fields:
                fields = csv_reader.fieldnames
            elif ("site" in csv_reader.fieldnames) and ("site" not in fields):
                fields = csv_reader.fieldnames

            try:
            
                sub_info_dict = csv_reader.next()
                
                if sub_info_dict["subject"] not in sub:
                    print "WARNING: data says sub id: %s" \
                          % sub_info_dict["subject"]
                    print "filepath says: %s" % sub
                    print "\n"
            
            
                sub_nums_list.append(sub_info_dict)

            except:

                print "didn't work: %s" % sub
                pass


    sub_nums_list = sorted(sub_nums_list, key=lambda k: k['subject']) 


    with open(outfile, "wt") as out_f:               

        csv_writer = csv.DictWriter(out_f, fields)

        csv_writer.writeheader()

        for row in sub_nums_list:
            csv_writer.writerow(row)



def merge_qap_outputs(output_directory):

    import glob

    spatial_outputs = glob.glob("%s/*/*/*/qap_anatomical_spatial*/*.csv" \
                                % output_directory)

    spatial_epi_outputs = glob.glob("%s/*/*/*/qap_functional_spatial*/*.csv" \
                                    % output_directory)

    temporal_outputs = glob.glob("%s/*/*/*/qap_functional_temporal*/*.csv" \
                                 % output_directory)


    if spatial_outputs:
        csv_merge(spatial_outputs, "qap_anatomical_spatial_%s.csv" % \
                                   output_directory.split("/")[-1])


    if spatial_epi_outputs:
        csv_merge(spatial_epi_outputs, "qap_functional_spatial_%s.csv" % \
                                       output_directory.split("/")[-1])


    if temporal_outputs:
        csv_merge(temporal_outputs, "qap_functional_temporal_%s.csv" % \
                                    output_directory.split("/")[-1])



def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("qap_output_directory", type=str, \
                            help="path to directory of subject output " \
                                 "folders containing QAP measures CSVs")

    args = parser.parse_args()

    # run it!
    merge_qap_outputs(args.qap_output_directory)



if __name__ == "__main__":
    main()


