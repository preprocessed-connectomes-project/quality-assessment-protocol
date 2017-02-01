#!/usr/bin/env python


def main():

    import argparse
    from qap.script_utils import gather_json_info, json_to_csv

    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=str,
                        help="the main output directory of the QAP run")
    parser.add_argument("--with-reports", action='store_true',
                        default=False, help="Write a summary report in PDF "
                        "format.")

    args = parser.parse_args()

    json_dict = gather_json_info(args.output_dir)
    json_to_csv(json_dict)


if __name__ == "__main__":
    main()