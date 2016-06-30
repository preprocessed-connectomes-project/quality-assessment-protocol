#!/user/bin/env python

def populate_data_dict(data_dict, filepath, part_id, session_id, series_id,
	site_id, resource, default_series_label=None):

    new_dict = data_dict

    if not session_id:
        session_id = "session_1"
    if not series_id:
    	if default_series_label:
    		series_id = default_series_label
    	else:
    		series_id = "scan_1"

    if part_id not in new_dict.keys():
        new_dict[part_id] = {}
    if session_id not in new_dict[part_id].keys():
        new_dict[part_id][session_id] = {}
    if resource not in new_dict[part_id][session_id].keys():
        new_dict[part_id][session_id][resource] = {}
    if site_id not in new_dict[part_id][session_id].keys():
        new_dict[part_id][session_id]["site_name"] = site_id
    if series_id not in new_dict[part_id][session_id][resource].keys():
        new_dict[part_id][session_id][resource][series_id] = filepath    

    data_dict.update(new_dict)

    return data_dict



def gather_raw_data(base_folder, directory_format, anatomical_keywords=None,
	functional_keywords=None):

    import os

    if "{participant}" not in directory_format:
    	pass

    data_dict = {}

    format_list = [x for x in directory_format.split("/") if x != ""]

    indices = {}
    operators = ["{site}", "{participant}", "{session}", "{series}"]

    for op in operators:
    	if op in format_list:
    		indices[op] = format_list.index(op)

    if anatomical_keywords:
        anatomical_keywords = [x for x in anatomical_keywords.split(" ") if x != ""]

    if functional_keywords:
        functional_keywords = [x for x in functional_keywords.split(" ") if x != ""]

    for root, dirs, files in os.walk(base_folder):
        for filename in files:
            if ".nii" in filename:

                filepath = os.path.join(root, filename)

                second_half_list = [x for x in filepath.split(base_folder)[1].split("/") if x != ""]

                site_id = None
                part_id = None
                session_id = None
                series_id = None

                if "{site}" in indices.keys():
                	site_id = second_half_list[indices["{site}"]]
                if "{participant}" in indices.keys():
                	part_id = second_half_list[indices["{participant}"]]
                if "{session}" in indices.keys():
                	session_id = second_half_list[indices["{session}"]]
                if "{series}" in indices.keys():
                	series_id = second_half_list[indices["{series}"]]

                if anatomical_keywords:
                    for word in anatomical_keywords:
                        if (word in filename) or (word in session_id):

                            data_dict = populate_data_dict(data_dict,
                                                           filepath, part_id,
                                                           session_id,
                            	                           series_id,
                                                           site_id,
                                                           "anatomical_scan",
                                                           "anat_1")

                if functional_keywords:
                    for word in functional_keywords:
                        if (word in filename) or (word in session_id):

                            data_dict = populate_data_dict(data_dict,
                                                           filepath, part_id,
                                                           session_id,
                                                           series_id,
                                                           site_id,
                                                           "functional_scan",
                                                           "func_1")

    if len(data_dict) == 0:
    	err = "\n\n[!] No data files were found given the inputs you " \
    	      "provided! Double-check your setup.\n\n"
    	raise Exception(err)

    return data_dict



def write_yaml_file(data_dict, out_file):

    import yaml

    if "." in out_file:
        out_file = out_file.split(".")[0]

    out_file = out_file + ".yml"

    with open(out_file,"wt") as f:
        f.write(yaml.dump(data_dict))



def main():

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
    data_dict = gather_raw_data(args.base_directory, args.directory_format,
        args.anatomical_keywords, args.functional_keywords)

    write_yaml_file(data_dict, args.outfile_path)



if __name__ == "__main__":
    main()