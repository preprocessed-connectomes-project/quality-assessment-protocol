#!/usr/bin/env python

def download_outputs(path_prefix, creds_path, bucket_name, qap_type, \
                          download_to):
                          
    import pickle
    from indi_aws import fetch_creds
    from indi_aws.aws_utils import s3_download

    src_list = []

    bucket = fetch_creds.return_bucket(creds_path, bucket_name)


    if qap_type == "anat_spatial":
        search_for = "anatomical_spatial"
    elif qap_type == "func_spatial":
        search_for = "functional_spatial"
    elif qap_type == "func_temporal":
        search_for = "functional_temporal"


    for k in bucket.list(prefix=path_prefix):

        k_name = str(k.name)
    
        if (search_for in k_name) and (".csv" in k_name):
    
            src_list.append(k_name)
        
 
    s3_download(bucket, src_list, download_to)
    
    

def main():

    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("path_prefix", type=str, \
                            help="path to the directory on your S3 bucket " \
                                 "containing the subject output directories")

    parser.add_argument("creds_path", type=str, \
                            help="path to file containing your AWS " \
                                 "credential keys")
                            
    parser.add_argument("bucket_name", type=str, \
                            help="name of your S3 bucket")
                                 
    parser.add_argument("qap_type", type=str, \
                            help="'anat_spatial', 'func_spatial', or 'func_" \
                                 "temporal', which type of QAP measures out" \
                                 "put you wish to download")
                                 
    parser.add_argument("download_to", type=str, \
                            help="full path to where you want to save the " \
                                 "downloaded output CSV files")

    args = parser.parse_args()

    # run it!
    download_outputs(args.path_prefix, args.creds_path, args.bucket_name, \
                        args.qap_type, args.download_to)



if __name__ == "__main__":
    main()
    
    
