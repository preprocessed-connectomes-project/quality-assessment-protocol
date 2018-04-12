#!/usr/bin/env python


def main():

    import os
    import argparse
    from qap.script_utils import check_csv_missing_subs, csv_to_pandas_df, \
        write_inputs_dict_to_yaml_file, read_yml_file
    from qap.qap_utils import raise_smart_exception

    parser = argparse.ArgumentParser()
    parser.add_argument("output_csv", type=str,
                        help="the main output directory of the QAP run "
                        "which contains the participant directories")
    parser.add_argument("data_config", type=str,
                        help="the main output directory of the QAP run "
                        "which contains the participant directories")
    parser.add_argument("data_type", type=str,
                        help="the main output directory of the QAP run "
                        "which contains the participant directories")

    args = parser.parse_args()

    csv_df = csv_to_pandas_df(args.output_csv)
    data_dict = read_yml_file(args.data_config)

    new_dict = check_csv_missing_subs(csv_df, data_dict, args.data_type)

    if new_dict:
        out_file = os.path.join(os.getcwd(),
                                "missing_%s_data.yml" % args.data_type)
        write_inputs_dict_to_yaml_file(new_dict, out_file)


if __name__ == "__main__":
    main()