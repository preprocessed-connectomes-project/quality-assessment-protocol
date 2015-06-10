
# this script takes in all of the generated QAP numbers from the QAP pipeline
# output directory you provide it and aggregates them into one CSV file

# the QAP pipelines already create their own combined CSV file, this script is
# for when you need to do it manually

# script usage:
#     > python merge_qap_outputs.py <qap pipeline output directory>

def csv_merge(csv_list, outfile):

    import csv

    fields = []

    sub_nums_list = []

    for sub in csv_list:

        with open(sub, "r") as in_f:

            csv_reader = csv.DictReader(in_f)

            if not fields:
                fields = csv_reader.fieldnames

            try:
                sub_nums_list.append(csv_reader.next())
            except:
                print "didn't work: %s" % sub
                pass


    with open(outfile, "wt") as out_f:               

        csv_writer = csv.DictWriter(out_f, fields)

        csv_writer.writeheader()

        for row in sub_nums_list:
            csv_writer.writerow(row)



def merge_qap_outputs(output_directory):

    import glob

    spatial_outputs = glob.glob("%s/*/*/*/qap_spatial/*" % output_directory)

    spatial_epi_outputs = glob.glob("%s/*/*/*/qap_spatial_epi/*" % \
                          output_directory)

    temporal_outputs = glob.glob("%s/*/*/*/qap_temporal/*" % output_directory)


    if spatial_outputs:
        csv_merge(spatial_outputs, "qap_spatial_%s.csv" % \
                                   output_directory.split("/")[-1])


    if spatial_epi_outputs:
        csv_merge(spatial_epi_outputs, "qap_spatial_epi_%s.csv" % \
                                       output_directory.split("/")[-1])


    if temporal_outputs:
        csv_merge(temporal_outputs, "qap_temporal_%s.csv" % \
                                    output_directory.split("/")[-1])



def run():

    import sys

    output_dir = sys.argv[1]

    merge_qap_outputs(output_dir)



run()


