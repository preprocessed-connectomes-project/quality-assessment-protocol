
def pull_S3_sublist(yaml_outpath, img_type, bucket_name, bucket_prefix, creds_path):

    import os
    from CPAC.AWS import fetch_creds
    import yaml

    s3_list = []
    s3_dict = {}

    bucket = fetch_creds.return_bucket(creds_path, bucket_name)

    # Filter for anat/rest
    if img_type == 'anat':
        subkey_type = 'anatomical_scan'
    elif img_type == 'rest':
        subkey_type = 'functional_scan'


    # Build S3-subjects to download
    for bk in bucket.list(prefix=bucket_prefix):
        s3_list.append(str(bk.name))

    # Build dictionary of filepaths
    for sfile in s3_list:

        ssplit = sfile.split('/')

        sub_id = ssplit[-4]

        session_id = ssplit[-3]

        scan_id = ssplit[-2]

        if img_type in scan_id:
             
            # this ONLY handles raw data inputs, not CPAC-generated outputs!
            if not s3_dict.has_key((sub_id, session_id, scan_id)):

                resource_dict = {}
                resource_dict[subkey_type] = sfile

                s3_dict[(sub_id, session_id, scan_id)] = {}
                s3_dict[(sub_id, session_id, scan_id)].update(resource_dict)

            else:

                s3_dict[(sub_id, session_id, scan_id)].update(resource_dict)         

        else:

            continue
    
            
    if len(s3_dict) == 0:
        err = "\n[!] Filepaths have not been successfully gathered from " \
              "the S3 bucket!\n"
        raise Exception(err)

    dict_len = len(s3_dict)            
           
    # write yaml file
    with open(yaml_outpath,"wt") as f:
        f.write(yaml.dump(s3_dict))
        
    
    if os.path.isfile(yaml_outpath):
        print "\nS3 dictionary file successfully created: %s\n" % yaml_outpath
        print "Total number of subject-session-scans: %d\n" % dict_len
    else:
        err = "\n[!] Filepaths from the S3 bucket have not been " \
              "successfully saved to the YAML file!\nOutput filepath: %s\n" \
              % yaml_outpath
        raise Exception(err)



def main():

    import argparse

    parser = argparse.ArgumentParser()
                       
    parser.add_argument("scan_type", type=str, \
                            help="'anat' or 'func', depending on which QAP " \
                                 "measures you will be using the S3 subject "\
                                 "dictionary for")
 
    parser.add_argument("bucket_name", type=str, \
                            help="the name of your AWS S3 bucket")
 
    parser.add_argument("bucket_prefix", type=str, \
                            help="the filepath prefix to the top level of " \
                                 "your raw data directory on S3 storage")
 
    parser.add_argument("creds_path", type=str, \
                            help="the path to the file containing your AWS " \
                                 "credentials")

    parser.add_argument("outfile_path", type=str, \
                            help="the full filepath for the S3 subject " \
                                 "YAML dictionary this script will create")
 
    args = parser.parse_args()


    # run it!
    pull_S3_sublist(args.outfile_path, args.scan_type, args.bucket_name, \
                        args.bucket_prefix, args.creds_path)



if __name__ == "__main__":
    main()
    
    
