#!/usr/bin/env python


def main():

    import os
    import argparse
    from qap.script_utils import gather_json_info, json_to_csv
    from qap.qap_utils import raise_smart_exception

    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str,
                        help="the main output directory of the QAP run "
                        "which contains the participant directories")
    parser.add_argument("--with_group_reports", action='store_true',
                        default=False, help="Write a summary report in PDF "
                        "format.")
    parser.add_argument("--with_full_reports", action='store_true',
                        default=False, help="Write the summary report and "
                        "the individual participant reports as well.")

    args = parser.parse_args()

    json_dict = gather_json_info(args.output_dir)
    json_to_csv(json_dict)

    if args.with_group_reports or args.with_full_reports:
        from qap.viz.reports import workflow_report
        # TODO: read in combined results dictionary from log JSONs
        #logs_dict = gather_json_info("_".join([args.output_dir, "logs"]))

        qap_types = ["anatomical_spatial",
                     "functional_spatial",
                     "functional_temporal"]
        for qap_type in qap_types:
            qap_type = "_".join(["qap", qap_type])
            run_name = args.output_dir.split("/")[-1]
            in_csv = os.path.join(os.getcwd(), '%s.csv' % qap_type)
            if not os.path.isfile(in_csv):
                continue
            reports = workflow_report(in_csv, qap_type, run_name,
                                      out_dir=args.output_dir,
                                      full_reports=args.with_full_reports)


if __name__ == "__main__":
    main()