#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("anatomical_scan", type=str,
                    help="path to the raw anatomical scan")
                    
parser.add_argument("head_template", type=str,
                    help="path to the anatomical template used for " \
                         "anatomical registration")
                    
parser.add_argument("participant_id", type=str,
                    help="participant ID of the anatomical scan (used for " \
                    	 "identification in the output CSV)")
                         
parser.add_argument("--session_id", type=str, default=None,
                    help="session ID of the anatomical scan (used for " \
                    	 "identification in the output CSV)")

parser.add_argument("--series_id", type=str, default=None,
                    help="series ID of the anatomical scan (used for " \
                    	 "identification in the output CSV)")

parser.add_argument("--site_name", type=str, default=None,
                    help="site name of the anatomical scan (used for " \
                    	 "identification in the output CSV)")

parser.add_argument("--output_directory", type=str, default=None,
                    help="where to write the output CSV (default: current " \
                    	 "directory)")
                    
args = parser.parse_args()



from qap import qap_workflows as qw

qw.run_everything_qap_anatomical_spatial(args.anatomical_scan,
    	                                 args.head_template,
    	                                 str(args.participant_id),
    	                                 args.session_id,
    	                                 args.series_id,
    	                                 args.site_name,
    	                                 args.output_directory)

