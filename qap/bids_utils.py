
def gather_nifti_file_paths(dataset_folder, creds_path=None):
    """Gather the NIFTI filepaths present on either an Amazon S3 bucket, or
    on the local filesystem.

    Keyword Arguments:
      dataset_folder -- [string] the path to the directory containing the data
      creds_path -- [string] (default: None) the path to the file containing 
                    your AWS credentials

    Returns:
      file_path_list -- [Python list] a list of the gathered filepaths
    """

    import os

    s3_prefix="s3://s3.amazonaws.com"

    file_path_list=[]

    # paths that include s3:// are assumed to live in AWS Simple Storage Service
    if "s3://" in dataset_folder:
        try:
            from indi_aws import fetch_creds
        except Exception as e:
            print "Error ({0:s}): Could not import indi_aws package".format(e.message)
            raise(e)

        try:
            s3_path_vals=(dataset_folder.replace(s3_prefix,"")).split('/')
            bucket_name = s3_path_vals[1]
            data_path="/".join(s3_path_vals[2:])
        except Exception as e:
            print "Error ({0:s}): There is a problem with s3 path {1:s}".format(e.message,dataset_folder)
            raise(e)

        print "Extracting NIfTI paths from s3 bucket {0:s}::{1:s})".format(bucket_name,data_path)

        bucket = fetch_creds.return_bucket(creds_path, bucket_name)

        # Build S3-subjects to download
        for bk in bucket.objects.filter(Prefix=data_path):
            if str(bk.key).endswith(".nii") or str(bk.key).endswith(".nii.gz"):
                file_path_list.append(os.path.join(s3_prefix,bucket_name,str(bk.key)))

    else:

        print "Extracting NIfTI paths from local filesystem"
        for root, folders, files in os.walk(os.path.abspath(dataset_folder)):
            for filename in files:
                if filename.endswith('.nii') or filename.endswith('.nii.gz'):
                    file_path_list.append(os.path.join(root,filename))

    if not file_path_list:
        raise Exception( "Did not find any nifti files in %s"%(dataset_folder) )

    return(file_path_list)


def extract_bids_data( file_path_list, inclusion_list=None ):
    """Parse the BIDS-formatted data from a list of filepaths and create a
    dictionary mapping the data to participant, session, series, etc.

    Keyword Arguments:
      file_path_list -- [Python list] a list of existing NIFTI filepaths
      inclusion_list -- [Python list] (default: None) a list of participant 
                        IDs, to prune down the dictionary to only include 
                        these participants

    Returns:
      sub_dict -- [Python dictionary] a dictionary containing the NIFTI 
                  filepaths mapped to their participant information and type 
                  of scan

    Notes:
      - For more information on the BIDS data structure format, visit:
          http://bids.neuroimaging.io/
    """

    import os

    # iterate through the files and put them into a dictionary, all of the information that we need to do this
    # is in the filename
    sub_dict = {}

    for file_path in file_path_list:
        filename=os.path.basename(file_path)
        try:
            # discard the file extension and split filename into key-value chunks, the last chunk being the series type
            f_chunks = (filename.split(".")[0]).split("_")
            # make a dictionary from the key-value chunks
            f_dict = {chunk.split("-")[0]:"-".join(chunk.split("-")[1:]) for chunk in f_chunks[:-1]}
            f_dict["series"] = f_chunks[-1]
        except Exception as e:
            print "Error (%s): gather_bids_data, %s does not appear to be in BIDS format"%(e.message,filename)
            continue

        if "sub" not in f_dict.keys():
            print u"Error (missing 'sub-' key): {0:s} does not appear to be in BIDS format".format(filename)
            continue

        # add sub- onto the front of subject name to preserve bids-ness
        f_dict["sub"] = "-".join(["sub",f_dict["sub"]])

        # straighten out the session
        if "ses" not in f_dict.keys():
            f_dict["ses"] = "1"

        # add ses- onto the front of session name to preserve bids-ness
        f_dict["ses"] = "-".join(["ses",f_dict["ses"]])

        # determine whether the scan is anatomical or functional, we don't know how to handle anything but T1w
        # and BOLD for now
        if "t1w" in f_dict["series"].lower():
            scan_type = "anatomical_scan"
        elif "bold" in f_dict["series"].lower():
            scan_type = "functional_scan"
        else:
            print u"QAP currently does not support {0:s} scans".format(f_dict["series"])
            continue

        if not inclusion_list or f_dict["sub"] in inclusion_list:

            # make sure that our different levels of dictionaries exist
            if f_dict["sub"] not in sub_dict.keys():
                sub_dict[f_dict["sub"]] = {}
            if f_dict["ses"] not in sub_dict[f_dict["sub"]].keys():
                sub_dict[f_dict["sub"]][f_dict["ses"]]={}
            if scan_type not in sub_dict[f_dict["sub"]][f_dict["ses"]].keys():
                sub_dict[f_dict["sub"]][f_dict["ses"]][scan_type]={}

            # calculate a key for the scan file
            f_key=f_dict["series"]
            if "functional_scan" in scan_type:
                if not "task" in f_dict:
                    print "Error (missing 'task-' key), Functional scan {0:s}".format(filename) + \
                          " does not appear to be in BIDS format"
                    continue
                else:
                    f_key = "_".join([f_key,"-".join(["task",f_dict["task"]])])
            if "acq" in f_dict:
                f_key = "_".join([f_key,"-".join(["acq",f_dict["acq"]])])

            # insert the full path to the scan into dictionary
            sub_dict[f_dict["sub"]][f_dict["ses"]][scan_type][f_key] = file_path

    if len(sub_dict) == 0:
        err = "\nThe participant list came out empty! Double-check your " \
              "settings.\n"
        raise Exception(err)

    return sub_dict

    return(sub_dict)
