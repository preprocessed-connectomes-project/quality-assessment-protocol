#!/usr/bin/env python
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
import os
import os.path as op
import sys


def build_anatomical_spatial_workflow(
        resource_pool, config, subject_info, run_name, site_name=None):

    # build pipeline for each subject, individually
    # ~ 29 minutes per subject with 1 core to ANTS

    import os
    import os.path as op
    import sys
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu
    import nipype.interfaces.fsl.maths as fsl

    import glob

    import time
    from time import strftime
    from nipype import config as nyconfig
    from nipype import logging
    logger = logging.getLogger('workflow')

    sub_id = str(subject_info[0])

    if subject_info[1]:
        session_id = subject_info[1]
    else:
        session_id = "session_0"

    if subject_info[2]:
        scan_id = subject_info[2]
    else:
        scan_id = "scan_0"

    # define and create the output directory
    output_dir = op.join(config["output_directory"], run_name,
                         sub_id, session_id, scan_id)

    try:
        os.makedirs(output_dir)
    except:
        if not op.isdir(output_dir):
            err = "[!] Output directory unable to be created.\n" \
                  "Path: %s\n\n" % output_dir
            raise Exception(err)
        else:
            pass

    log_dir = output_dir

    # set up logging
    nyconfig.update_config(
        {'logging': {'log_directory': log_dir, 'log_to_file': True}})
    logging.update_logging(nyconfig)

    # take date+time stamp for run identification purposes
    unique_pipeline_id = strftime("%Y%m%d%H%M%S")
    pipeline_start_stamp = strftime("%Y-%m-%d_%H:%M:%S")

    pipeline_start_time = time.time()

    logger.info(
        ("Pipeline start time: %s.\nContents of resource pool: %s.\n"
         "Configuration settings: %s.\n") % (
            pipeline_start_stamp, str(resource_pool), str(config)))

    # for QAP spreadsheet generation onlygit
    config.update({'subject_id': sub_id, 'session_id': session_id,
                   'scan_id': scan_id, 'run_name': run_name})

    if site_name:
        config["site_name"] = site_name

    # os.environ['ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS'] =
    #    str(config["num_ants_threads"])

    workflow = pe.Workflow(name=scan_id)
    workflow.base_dir = op.join(config["working_directory"], sub_id,
                                session_id)

    # set up crash directory
    workflow.config['execution'] = {
        'crashdump_dir': config["output_directory"]}

    # update that resource pool with what's already in the output directory
    for resource in os.listdir(output_dir):
        if op.isdir(op.join(output_dir, resource)) and \
                resource not in resource_pool.keys():

            resource_pool[resource] = glob.glob(op.join(output_dir,
                                                        resource, "*"))[0]

    # resource pool check
    invalid_paths = []

    for resource in resource_pool.keys():
        if not op.isfile(resource_pool[resource]):
            invalid_paths.append((resource, resource_pool[resource]))

    if len(invalid_paths) > 0:
        err = "\n\n[!] The paths provided in the subject list to the " \
              "following resources are not valid:\n"

        for path_tuple in invalid_paths:
            err = err + path_tuple[0] + ": " + path_tuple[1] + "\n"

        err = err + "\n\n"
        raise Exception(err)

    # start connecting the pipeline
    if "qap_anatomical_spatial" not in resource_pool.keys():
        from qap.qap_workflows import qap_anatomical_spatial_workflow
        workflow, resource_pool = qap_anatomical_spatial_workflow(
            workflow, resource_pool, config)

    # set up the datasinks
    new_outputs = 0

    if "write_all_outputs" not in config.keys():
        config["write_all_outputs"] = False

    if config["write_all_outputs"]:
        for output in resource_pool.keys():
            # we use a check for len()==2 here to select those items in the
            # resource pool which are tuples of (node, node_output), instead
            # of the items which are straight paths to files

            # resource pool items which are in the tuple format are the
            # outputs that have been created in this workflow because they
            # were not present in the subject list YML (the starting resource
            # pool) and had to be generated

            if len(resource_pool[output]) == 2:
                ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
                ds.inputs.base_directory = output_dir
                node, out_file = resource_pool[output]
                workflow.connect(node, out_file, ds, output)
                new_outputs += 1
    else:
        # write out only the output CSV (default)
        output = "qap_anatomical_spatial"
        if len(resource_pool[output]) == 2:
            ds = pe.Node(nio.DataSink(), name='datasink_%s' % output)
            ds.inputs.base_directory = output_dir
            node, out_file = resource_pool[output]
            workflow.connect(node, out_file, ds, output)
            new_outputs += 1

    # run the pipeline (if there is anything to do)
    if new_outputs > 0:
        workflow.write_graph(
            dotfilename=op.join(output_dir, run_name + ".dot"),
            simple_form=False)

        if config["num_cores_per_subject"] == 1:
            workflow.run(plugin='Linear')
        else:
            workflow.run(
                plugin='MultiProc',
                plugin_args={'n_procs': config["num_cores_per_subject"]})
    else:
        print "\nEverything is already done for subject %s." % sub_id

    # Remove working directory when done
    if not config["write_all_outputs"]:
        try:
            work_dir = op.join(workflow.base_dir, scan_id)
            if op.exists(work_dir):
                import shutil
                shutil.rmtree(work_dir)
        except:
            print "Couldn\'t remove the working directory!"
            pass

    pipeline_end_stamp = strftime("%Y-%m-%d_%H:%M:%S")
    logger.info("Pipeline end time: %s" % pipeline_end_stamp)
    pipeline_end_time = time.time()
    logger.info("Elapsed time (minutes) since last start: %s"
                % ((pipeline_end_time - pipeline_start_time) / 60))
    return workflow


def run(subject_list, config, cloudify=False):
    import os
    import os.path as op
    import yaml
    import time
    from multiprocessing import Process

    from nipype import logging
    logger = logging.getLogger('workflow')

    write_report = config.get('write_report', False)
    if write_report:
        logger.info('PDF Reports enabled')

    with open(subject_list, "r") as f:
        subdict = yaml.load(f)

    if not cloudify:
        flat_sub_dict = {}
        sites_dict = {}

        for subid in subdict.keys():
            # sessions
            for session in subdict[subid].keys():
                # resource files
                for resource in subdict[subid][session].keys():
                    if type(subdict[subid][session][resource]) is dict:
                        # then this has sub-scans defined
                        for scan in subdict[subid][session][resource].keys():
                            filepath = subdict[subid][session][resource][scan]
                            resource_dict = {}
                            resource_dict[resource] = filepath
                            sub_info_tuple = (subid, session, scan)

                            if sub_info_tuple not in flat_sub_dict.keys():
                                flat_sub_dict[sub_info_tuple] = {}

                            flat_sub_dict[sub_info_tuple].update(resource_dict)
                    elif resource == "site_name":
                        sites_dict[subid] = subdict[subid][session][resource]

                    else:
                        filepath = subdict[subid][session][resource]
                        resource_dict = {}
                        resource_dict[resource] = filepath
                        sub_info_tuple = (subid, session, None)

                        if sub_info_tuple not in flat_sub_dict.keys():
                            flat_sub_dict[sub_info_tuple] = {}

                        flat_sub_dict[sub_info_tuple].update(resource_dict)

    try:
        os.makedirs(config["output_directory"])
    except:
        if not op.isdir(config["output_directory"]):
            err = ("[!] Output directory unable to be created.\n"
                   "Path: %s\n\n" % config["output_directory"])
            raise Exception(err)
        else:
            pass

    try:
        os.makedirs(config["working_directory"])
    except:
        if not op.isdir(config["working_directory"]):
            err = "[!] Output directory unable to be created.\n" \
                  "Path: %s\n\n" % config["working_directory"]
            raise Exception(err)
        else:
            pass

    # get the pipeline config file name, use it as the run name
    run_name = config['pipeline_config_yaml'].split("/")[-1].split(".")[0]
    ns_at_once = config.get('num_subjects_at_once', 1)

    if not cloudify:
        # skip parallel machinery if we are running only one subject at once
        if ns_at_once == 1:
            for sub_info in flat_sub_dict.keys():
                if sites_dict:
                    site = sites_dict[sub_info[0]]
                else:
                    site = None
                if 'anatomical_scan' in flat_sub_dict[sub_info].keys():
                    build_anatomical_spatial_workflow(
                        flat_sub_dict[sub_info], config, sub_info,
                        run_name, site)
        else:
            if 'anatomical_scan' in flat_sub_dict[sub_info].keys():
                if len(sites_dict) > 0:
                    procss = [Process(
                        target=build_anatomical_spatial_workflow,
                        args=(flat_sub_dict[sub_info], config, sub_info,
                              run_name, sites_dict[sub_info[0]]))
                              for sub_info in flat_sub_dict.keys()]

                elif len(sites_dict) == 0:
                    procss = [Process(
                        target=build_anatomical_spatial_workflow,
                        args=(flat_sub_dict[sub_info], config, sub_info,
                              run_name, None))
                              for sub_info in flat_sub_dict.keys()]

                pid = open(op.join(
                    config["output_directory"], 'pid.txt'), 'w')

                # Init job queue
                job_queue = []

                # Stream the subject workflows for preprocessing.
                # At Any time in the pipeline c.numSubjectsAtOnce
                # will run, unless the number remaining is less than
                # the value of the parameter stated above

                idx = 0
                nprocs = len(procss)
                while idx < nprocs:
                    # Check every job in the queue's status
                    for job in job_queue:
                        # If the job is not alive
                        if not job.is_alive():
                            # Find job and delete it from queue
                            logger.info('found dead job: %s' % str(job))
                            loc = job_queue.index(job)
                            del job_queue[loc]

                    # Check free slots after prunning jobs
                    slots = ns_at_once - len(job_queue)

                    if slots > 0:
                        idc = idx
                        for p in procss[idc:idc + slots]:
                            # ..and start the next available process (subject)
                            p.start()
                            print >>pid, p.pid
                            # Append this to job queue and increment index
                            job_queue.append(p)
                            idx += 1

                    # Add sleep so while loop isn't consuming 100% of CPU
                    time.sleep(2)
                pid.close()

                # Join all processes if report must be written out
                if write_report:
                    for p in procss:
                        p.join()
    else:
        # run on cloud
        sub = subject_list.keys()[0]

        # get the site name!
        for resource_path in subject_list[sub]:
            if ".nii" in resource_path:
                filepath = resource_path
                break

        filesplit = filepath.split(config["bucket_prefix"])
        site_name = filesplit[1].split("/")[1]

        build_anatomical_spatial_workflow(
            subject_list[sub], config, sub, run_name, site_name)

    # PDF reporting
    if write_report:
        import os.path as op
        import qap.viz.reports as qvr

        in_csv = op.join(
            config['output_directory'], 'qap_anatomical_spatial.csv')
        out_file = op.join(
            config['output_directory'], 'qap_anatomical.pdf')

        qvr.report_anatomical(in_csv, out_file=out_file)


# Main routine
def main():
    import argparse
    import yaml

    parser = argparse.ArgumentParser()
    group = parser.add_argument_group("Regular Use Inputs (non-cloud runs)")
    cloudgroup = parser.add_argument_group("AWS Cloud Inputs (only required "
                                           "for AWS Cloud runs)")
    req = parser.add_argument_group("Required Inputs")
    cloudgroup.add_argument('--subj_idx', type=int,
                            help='Subject index to run')
    cloudgroup.add_argument('--s3_dict_yml', type=str,
                            help='Path to YAML file containing S3 input '
                            'filepaths dictionary')

    # Subject list (YAML file)
    group.add_argument("--sublist", type=str,
                       help="filepath to subject list YAML")
    req.add_argument("config", type=str,
                     help="filepath to pipeline configuration YAML")
    # Write PDF reports
    group.add_argument("--with-reports", action='store_true', default=False,
                       help="Write a summary report in PDF format.")

    args = parser.parse_args()

    # checks
    if args.subj_idx and not args.s3_dict_yml and not args.sublist:
        print("\n[!] You provided --subj_idx, but not --s3_dict_yml. When "
              "executing cloud-based runs, please provide both inputs.\n")

    elif args.s3_dict_yml and not args.subj_idx and not args.sublist:
        print("\n[!] You provided --s3_dict_yml, but not --subj_idx. When "
              "executing cloud-based runs, please provide both inputs.\n")

    elif not args.sublist and not args.subj_idx and not args.s3_dict_yml:
        print("\n[!] Either --sublist is required for regular runs, or both "
              "--subj_idx and --s3_dict_yml for cloud-based runs.\n")

    elif args.sublist and args.subj_idx and args.s3_dict_yml:
        print("\n[!] Either --sublist is required for regular runs, or both "
              "--subj_idx and --s3_dict_yml for cloud-based runs, but not "
              "all three. (I'm not sure which you are trying to do!)\n")

    elif args.sublist and (args.subj_idx or args.s3_dict_yml):
        print("\n[!] Either --sublist is required for regular runs, or both "
              "--subj_idx and --s3_dict_yml for cloud-based runs. (I'm not "
              "sure which you are trying to do!)\n")

    else:
        with open(args.config, "r") as f:
            config = yaml.load(f)

        config['pipeline_config_yaml'] = args.config

        if 'write_report' not in config:
            config['write_report'] = False

        if args.with_reports:
            config['write_report'] = True

        if args.subj_idx and args.s3_dict_yml:
            # ---- Cloud-ify! ----
            # Import packages
            from qap.cloud_utils import dl_subj_from_s3, upl_qap_output

            # Download and build a one-subject dictionary from S3
            sub_dict = dl_subj_from_s3(args.subj_idx, config,
                                       args.s3_dict_yml)

            if not sub_dict:
                err = ("\n[!] Subject dictionary was not successfully "
                       "downloaded from the S3 bucket!\n")
                raise Exception(err)

            # Run it
            run(sub_dict, config, cloudify=True)

            # Upload results
            upl_qap_output(args.config)

        elif args.sublist:
            # Run it
            run(args.sublist, config, cloudify=False)


# Make executable
if __name__ == "__main__":
    main()
