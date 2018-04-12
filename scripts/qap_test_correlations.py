#!/usr/bin/env python

def main():

    import os
    import argparse

    from qap.script_utils import csv_to_pandas_df, read_txt_file, \
                                 qap_csv_correlations

    parser = argparse.ArgumentParser()

    parser.add_argument("old_csv", type=str, \
                            help="path to the QAP CSV output file from a " \
                                 "previous run or version")

    parser.add_argument("new_csv", type=str, \
                            help="path to the QAP CSV output file from the " \
                                 "current run or version")
                                
    parser.add_argument("--replacements", type=str, \
                            help="text file containing column name pairs " \
                                 "to be updated/renamed in the old CSV - " \
                                 "ex. 'subject,Participant' on one line to " \
                                 "rename the 'subject' column to " \
                                 "'Participant'")

    args = parser.parse_args()

    # run it!
    old_data = csv_to_pandas_df(args.old_csv)
    new_data = csv_to_pandas_df(args.new_csv)

    if args.replacements:
    	replacements = read_txt_file(args.replacements)
    else:
        replacements = None

    correlations_dict = qap_csv_correlations(old_data, new_data, replacements)

    print "\nQAP Results Correlations (Pearson's r)"
    print "Previous dataset: %s" % os.path.abspath(args.old_csv)
    print "Current dataset: %s\n" % os.path.abspath(args.new_csv)
    for key in sorted(correlations_dict.keys()):
    	print "%s: %f" % (key, correlations_dict[key][0])
    print "\n"


if __name__ == "__main__":
    main()
    
