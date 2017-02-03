#!/usr/bin/env python


def main():

    import argparse
    from qap.script_utils import gather_json_info, json_to_csv

    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str,
                        help="the main output directory of the QAP run")
    parser.add_argument("--with_reports", action='store_true',
                        default=False, help="Write a summary report in PDF "
                        "format.")

    args = parser.parse_args()

    json_dict = gather_json_info(args.output_dir)
    json_to_csv(json_dict)

    if args.with_reports:
        from qap.viz.reports import workflow_report
        # TODO: read in combined results dictionary from log JSONs
        #logs_dict = gather_json_info("_".join([args.output_dir, "logs"]))
        # TODO: get a handle for the three CSVs, do checks

        qap_types = ["anatomical_spatial",
                     "functional_spatial",
                     "functional_temporal"]
        for qap_type in qap_types:
            qap_type = "_".join(["qap", qap_type])
            in_csv = op.join(config['output_directory'],
                             '%s.csv' % qap_type)
            reports = workflow_report(in_csv, qap_type, self._run_name,
                                      out_dir=config['output_directory'])


if __name__ == "__main__":
    main()