#!/user/bin/env python

def main():

    from qap.script_utils import gather_filepath_list, \
                                 gather_custom_raw_data, \
                                 write_inputs_dict_to_yaml_file

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("base_directory", type=str, \
                            help="full path to the directory holding the " \
                                 "raw data")

    parser.add_argument("outfile_path", type=str, \
                            help="full path for the generated data list")
                            
    parser.add_argument("directory_format", type=str, \
                            help="a directory layout description - place " \
                                 "{site}, {participant}, {session}, and " \
                                 "{series} where they occur in your " \
                                 "directory format separated by slashes - " \
                                 "example: '/{site}/{participant}/{session}/"\
                                 "{series}' - NOTE: if you don't have some " \
                                 "of these organization levels, they can be "\
                                 "omitted in the description")
                                 
    parser.add_argument('--anatomical_keywords', type=str, \
                            help="a string of keywords that should be found "\
                                 "in the series/scan IDs or filenames of " \
                                 "the anatomical scans")

    parser.add_argument('--functional_keywords', type=str, \
                            help="a string of keywords that should be found "\
                                 "in the series/scan IDs or filenames of " \
                                 "the functional scans")

    args = parser.parse_args()

    # run the thing
    filepath_list = gather_filepath_list(args.base_directory)

    data_dict = gather_custom_raw_data(filepath_list, args.base_directory, 
        args.directory_format, args.anatomical_keywords, 
        args.functional_keywords)

    write_inputs_dict_to_yaml_file(data_dict, args.outfile_path)


if __name__ == "__main__":
    main()