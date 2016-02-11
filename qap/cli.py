#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import os.path as op
import time
import argparse
import yaml

from nipype import logging
logger = logging.getLogger('workflow')


class QAProtocolCLI:
    """
    This class and the associated _run_workflow function implement what
    the former scripts (qap_anatomical_spatial.py, etc.) contained
    """

    def __init__(self):
        parser = argparse.ArgumentParser()

        group = parser.add_argument_group(
            "Regular Use Inputs (non-cloud runs)")
        cloudgroup = parser.add_argument_group(
            "AWS Cloud Inputs (only required for AWS Cloud runs)")
        req = parser.add_argument_group("Required Inputs")

        cloudgroup.add_argument('--subj_idx', type=int,
                                help='Subject index to run')
        cloudgroup.add_argument(
            '--s3_dict_yml', type=str,
            help='Path to YAML file containing S3 input filepaths dictionary')

        # Subject list (YAML file)
        group.add_argument(
            "--sublist", type=str, help="filepath to subject list YAML")
        req.add_argument(
            "config", type=str, help="filepath to pipeline configuration YAML")

        # Write PDF reports
        group.add_argument(
            "--with-reports", action='store_true', default=False,
            help="Write a summary report in PDF format.")

        args = parser.parse_args()

        # checks
        if args.subj_idx and not args.s3_dict_yml and not args.sublist:
            raise RuntimeError(
                "\n[!] You provided --subj_idx, but not --s3_dict_yml. "
                "When executing cloud-based runs, please provide both "
                "inputs.\n")

        elif args.s3_dict_yml and not args.subj_idx and not args.sublist:
            raise RuntimeError(
                "\n[!] You provided --s3_dict_yml, but not --subj_idx. "
                "When executing cloud-based runs, please provide both "
                "inputs.\n")

        elif not args.sublist and not args.subj_idx and not args.s3_dict_yml:
            raise RuntimeError(
                "\n[!] Either --sublist is required for regular runs, or "
                "both --subj_idx and --s3_dict_yml for cloud-based runs.\n")

        elif args.sublist and args.subj_idx and args.s3_dict_yml:
            raise RuntimeError(
                "\n[!] Either --sublist is required for regular runs, or "
                "both --subj_idx and --s3_dict_yml for cloud-based runs, "
                "but not all three. (I'm not sure which you are trying to "
                "do!)\n")

        elif args.sublist and (args.subj_idx or args.s3_dict_yml):
            raise RuntimeError(
                "\n[!] Either --sublist is required for regular runs, or "
                "both --subj_idx and --s3_dict_yml for cloud-based runs. "
                "(I'm not sure which you are trying to do!)\n")

        self._cloudify = False

        if args.subj_idx and args.s3_dict_yml:
            self._cloudify = True

            # ---- Cloud-ify! ----
            # Import packages
            from cloud_utils import dl_subj_from_s3, upl_qap_output
            # Download and build a one-subject dictionary from S3
            self._sub_dict = dl_subj_from_s3(
                args.subj_idx, args.config, args.s3_dict_yml)

            if not self._sub_dict:
                err = "\n[!] Subject dictionary was not successfully " \
                      "downloaded from the S3 bucket!\n"
                raise RuntimeError(err)

        elif args.sublist:
            self._sub_dict = args.sublist

        else:
            raise RuntimeError(
                "\n[!] Arguments were parsed, but no appropriate run found")

        # Load config
        with open(args.config, "r") as f:
            self._config = yaml.load(f)

        self._config['pipeline_config_yaml'] = args.config
        self._config['qap_type'] = parser.prog[4:-3]

        if args.with_reports:
            self._config['write_report'] = True

    def _run_here(self, run_name):
        ns_at_once = self._config.get('num_subjects_at_once', 1)
        with open(self._sub_dict, "r") as f:
            subdict = yaml.load(f)

        flat_sub_dict = {}
        sites_dict = {}

        # Preamble: generate flat_subdict
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

        # in case some subjects have site names and others don't
        if len(sites_dict.keys()) > 0:
            for subid in subdict.keys():
                if subid not in sites_dict.keys():
                    sites_dict[subid] = None

        # integrate site information into flat_sub_dict
        #     it was separate in the first place to circumvent the fact that
        #     even though site_name doesn't get keyed with scan names, that
        #     doesn't necessarily mean scan names haven't been specified for
        #     that participant
        for sub_info_tuple in flat_sub_dict.keys():
            site_info = {}
            site_info["site_name"] = sites_dict[sub_info_tuple[0]]
            flat_sub_dict[sub_info_tuple].update(site_info)

        # Start the magic
        logger.info('There are %d subjects in the pool' %
                    len(flat_sub_dict.keys()))

        results = _run_workflow(flat_sub_dict, flat_sub_dict.keys(), \
                                    self._config, run_name)


        return results



    def _run_cloud(self, run_name, subject_list):

        from cloud_utils import upl_qap_output
        # get the site name!

        sub = subject_list.keys()[0]

        for resource_path in subject_list[sub].values():
            if ".nii" in resource_path:
                filepath = resource_path
                break

        try:
            filesplit = filepath.split(self._config["bucket_prefix"])
            site_name = filesplit[1].split("/")[1]
            subject_list[sub]["site_name"] = site_name
        except:
            pass

        rt = _run_workflow(subject_list, sub, self._config, run_name)
	
        # make uploading results to S3 bucket the default if not specified
        if "upload_to_s3" not in self._config.keys():
            self._config["upload_to_s3"] = False

        # upload results
        if self._config["upload_to_s3"]:
            upl_qap_output(self._config)

        return rt

    def run(self):
        # Get configurations and settings
        config = self._config
        subject_list = self._sub_dict
        cloudify = self._cloudify
        ns_at_once = config.get('num_subjects_at_once', 1)
        write_report = config.get('write_report', False)

        # Create output directory
        try:
            os.makedirs(config["output_directory"])
        except:
            if not op.isdir(config["output_directory"]):
                err = "[!] Output directory unable to be created.\n" \
                      "Path: %s\n\n" % config["output_directory"]
                raise Exception(err)
            else:
                pass

        # Create working directory
        try:
            os.makedirs(config["working_directory"])
        except:
            if not op.isdir(config["working_directory"]):
                err = "[!] Output directory unable to be created.\n" \
                      "Path: %s\n\n" % config["working_directory"]
                raise Exception(err)
            else:
                pass

        run_name = config['pipeline_config_yaml'].split("/")[-1].split(".")[0]

        results = None

        #try:
        if not cloudify:
            results = self._run_here(run_name)
        else:
            results = self._run_cloud(run_name, subject_list)
        #except Exception as e:
        #    raise Exception("\n\n%s\n\n%s\n\n" % (e, locals()))

        # PDF reporting
        if write_report:
            from qap.viz.reports import workflow_report
            logger.info('Writing PDF reports')
            qap_type = 'qap_' + config['qap_type']
            in_csv = op.join(config['output_directory'], '%s.csv' % qap_type)

            reports = workflow_report(in_csv, qap_type, run_name, results,
                                      out_dir=config['output_directory'])

            for k, v in reports.iteritems():
                if v['success']:
                    logger.info('Written report (%s) in %s' % (k, v['path']))



def starter_node_func(starter):
    return starter



def _run_workflow(resource_pool_list, sub_info_list, config, run_name):

    # build pipeline for each subject, individually
    # ~ 5 min 20 sec per subject
    # (roughly 320 seconds)

    import os
    import os.path as op
    import sys

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu

    import nipype.interfaces.utility as util
    import nipype.interfaces.fsl.maths as fsl

    import qap

    import glob

    import time
    from time import strftime
    from nipype import config as nyconfig


    qap_type = config['qap_type']

    # Read and apply general settings in config
    keep_outputs = config.get('write_all_outputs', False)

    log_dir = op.join(config["output_directory"], run_name)

    try:
        os.makedirs(log_dir)
    except:
        if not op.isdir(log_dir):
            err = "[!] Output directory unable to be created.\n" \
                    "Path: %s\n\n" % log_dir
            raise Exception(err)
        else:
            pass


    # set up logging
    nyconfig.update_config(
        {'logging': {'log_directory': log_dir, 'log_to_file': True}})
    logging.update_logging(nyconfig)

    # take date+time stamp for run identification purposes
    unique_pipeline_id = strftime("%Y%m%d%H%M%S")
    pipeline_start_stamp = strftime("%Y-%m-%d_%H:%M:%S")

    pipeline_start_time = time.time()

    logger.info("QAP version %s" % qap.__version__)
    logger.info("Pipeline start time: %s" % pipeline_start_stamp)


    workflow = pe.Workflow(name=run_name)

    # set up crash directory
    workflow.config['execution'] = \
        {'crashdump_dir': config["output_directory"]}

    # create the one node all participants will start from
    starter_node = pe.Node(niu.Function(input_names=['starter'], 
                                        output_names=['starter'], 
                                        function=starter_node_func),
                           name='starter_node')

    # set a dummy variable
    starter_node.inputs.starter = ""


    new_outputs = 0

    for sub_info in sub_info_list:

        resource_pool = resource_pool_list[sub_info]

        # resource pool check

        invalid_paths = []

        for resource in resource_pool.keys():
            if not op.isfile(resource_pool[resource]) and \
                resource != "site_name":
                invalid_paths.append((resource, resource_pool[resource]))

        if len(invalid_paths) > 0:
            err = "\n\n[!] The paths provided in the subject list to the " \
                  "following resources are not valid:\n"

            for path_tuple in invalid_paths:
                err = err + path_tuple[0] + ": " + path_tuple[1] + "\n"

            err = err + "\n\n"
            raise Exception(err)


        resource_pool["starter"] = (starter_node, 'starter')

        # process subject info

        sub_id = str(sub_info[0])
        # for nipype
        if "-" in sub_id:
            sub_id = sub_id.replace("-","_")
        if "." in sub_id:
            sub_id = sub_id.replace(".","_")

        if sub_info[1]:
            session_id = sub_info[1]
            # for nipype
            if "-" in session_id:
                session_id = session_id.replace("-","_")
            if "." in session_id:
                session_id = session_id.replace(".","_")
        else:
            session_id = "session_0"

        if sub_info[2]:
            scan_id = sub_info[2]
            # for nipype
            if "-" in scan_id:
                scan_id = scan_id.replace("-","_")
            if "." in scan_id:
                scan_id = scan_id.replace(".","_")
        else:
            scan_id = "scan_0"


        # set output directory

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


        # individual workflow and logger setup

        logger.info("Contents of resource pool:\n" + str(resource_pool))
        logger.info("Configuration settings:\n" + str(config))

        # for QAP spreadsheet generation only
        config.update({"subject_id": sub_id, "session_id": session_id,
                       "scan_id": scan_id, "run_name": run_name})

        workflow.base_dir = op.join(config["working_directory"], sub_id,
                                session_id)


        # update that resource pool with what's already in the output directory

        for resource in os.listdir(output_dir):
            if (op.isdir(op.join(output_dir, resource)) and
                    resource not in resource_pool.keys()):
                resource_pool[resource] = glob.glob(op.join(output_dir,
                                                            resource, "*"))[0]


        # start connecting the pipeline
        if 'qap_' + qap_type not in resource_pool.keys():
            from qap import qap_workflows as qw
            wf_builder = getattr(qw, 'qap_' + qap_type + '_workflow')
            workflow, resource_pool = wf_builder(workflow, resource_pool, config)


        # set up the datasinks

        out_list = ['qap_' + qap_type]
        if keep_outputs:
            out_list = resource_pool.keys()

        # Save reports to out_dir if necessary
        if config.get('write_report', False):
            out_list += ['qap_mosaic']

            # The functional temporal also has an FD plot
            if 'functional_temporal' in qap_type:
                out_list += ['qap_fd']

        for output in out_list:
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

        rt = {'id': sub_id, 'session': session_id, 'scan': scan_id,
              'status': 'started'}


    # run the pipeline (if there is anything to do)
    if new_outputs > 0:
        if config.get('write_graph', False):
            workflow.write_graph(
                dotfilename=op.join(config["output_directory"], run_name + ".dot"),
                simple_form=False)

        runargs = {'plugin': 'Linear', 'plugin_args': {}}
        if nc_per_subject > 1:
            runargs['plugin'] = 'MultiProc'
            runargs['plugin_args'] = {'n_procs': num_subjects_at_once}

        try:
            workflow.run(**runargs)
            rt['status'] = 'finished'
        except Exception as e:  # TODO We should be more specific here ...
            rt.update({'status': 'failed', 'msg': e})
            # ... however this is run inside a pool.map: do not raise Exception

    else:
        rt['status'] = 'cached'
        logger.info("\nEverything is already done for subject %s." % sub_id)

    # Remove working directory when done
    if not keep_outputs:
        try:
            work_dir = op.join(workflow.base_dir, scan_id)

            if op.exists(work_dir):
                import shutil
                shutil.rmtree(work_dir)
        except:
            logger.warn("Couldn\'t remove the working directory!")
            pass

    pipeline_end_stamp = strftime("%Y-%m-%d_%H:%M:%S")
    pipeline_end_time = time.time()
    logger.info("Elapsed time (minutes) since last start: %s"
                % ((pipeline_end_time - pipeline_start_time) / 60))
    logger.info("Pipeline end time: %s" % pipeline_end_stamp)
    return rt
