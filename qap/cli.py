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
        cloudgroup.add_argument('--bundle_idx', type=int,
                                help='Bundle index to run')
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

        if args.bundle_idx and not args.s3_dict_yml and not args.sublist:
            raise RuntimeError(
                "\n[!] You provided --bundle_idx, but not --s3_dict_yml. "
                "When executing cloud-based runs, please provide both "
                "inputs.\n")

        elif not args.sublist and not args.subj_idx and not args.s3_dict_yml:
            raise RuntimeError(
                "\n[!] Either --sublist is required for regular runs, or "
                "both --subj_idx and --s3_dict_yml for cloud-based runs.\n")

        elif args.sublist and args.s3_dict_yml:
            raise RuntimeError(
                "\n[!] Either --sublist is required for regular runs, or "
                "both --s3_dict_yml and either --subj_idx or --bundle_idx " \
                "for cloud-based runs, but not all three. (I'm not sure " \
                "which you are trying to do!)\n")

        '''
        elif args.s3_dict_yml and not args.sublist:
            if not args.subj_idx and not args.bundle_idx:
                raise RuntimeError(
                    "\n[!] You provided --s3_dict_yml, but no --subj_idx or "\
                    "--bundle_idx. When executing cloud-based runs, please " \
                    "provide both required inputs.\n")
            if args.subj_idx and args.bundle_idx:
                raise RuntimeError(
                    "\n[!] You provided both --subj_idx and --bundle_idx. " \
                    "You only need one for the run.\n")
        '''

        '''
        elif args.sublist and (args.subj_idx or args.bundle_idx or \
            args.s3_dict_yml):
            raise RuntimeError(
                "\n[!] Either --sublist is required for regular runs, or "
                "both --s3_dict_yml and either --subj_idx or --bundle_idx " \
                "for cloud-based runs. (I'm not sure which you are trying " \
                "to do!)\n")
        '''

        # Load config
        with open(args.config, "r") as f:
            self._config = yaml.load(f)

        self._config['pipeline_config_yaml'] = args.config
        self._config['qap_type'] = parser.prog[4:-3]

        if args.with_reports:
            self._config['write_report'] = True

        if "num_subjects_per_bundle" not in self._config.keys():
            self._config["num_subjects_per_bundle"] = 1

        if "num_bundles_at_once" not in self._config.keys():
            self._config["num_bundles_at_once"] = 1


        self._cloudify = False

        if args.s3_dict_yml and (args.subj_idx or args.bundle_idx):

            # ---- Cloud-ify! ----
            self._cloudify = True

        elif args.s3_dict_yml:

            # do this alone in case we are sending in an S3 dict YAML but no
            # YAML entry indexes (like bundle_idx) - this happens when we
            # first kick off a cluster run
            self._s3_dict_yml = args.s3_dict_yml
            self._config["subject_list"] = None

        elif args.sublist:

            self._config["subject_list"] = args.sublist
            self._s3_dict_yml = None

        else:
            raise RuntimeError(
                "\n[!] Arguments were parsed, but no appropriate run found")

        
        if args.bundle_idx:

            self._bundle_idx = args.bundle_idx
            self._subj_idx = None

        elif args.subj_idx:

            self._subj_idx = args.subj_idx
            self._bundle_idx = None



    def _prepare_cluster_batch_file(self, run_name, num_bundles):

        from INDI_utils.indi_schedulers import cluster_templates

        # Cluster batch file setup code borrowed from dclark87's CPAC cluster
        # setup code
        #     https://github.com/FCP-INDI/C-PAC/blob/0.4.0_development/CPAC/..
        #         ..pipeline/cpac_runner.py
        #     https://github.com/dclark87

        # Create cluster log dir
        cluster_files_dir = \
            os.path.join(self._config["output_directory"], "cluster_files")
        if not os.path.exists(cluster_files_dir):
            os.makedirs(cluster_files_dir)

        # Batch file variables
        timestamp = str(strftime("%Y_%m_%d_%H_%M_%S"))
        shell = commands.getoutput('echo $SHELL')
        user_account = getpass.getuser()

        # Set up config dictionary
        config_dict = {'timestamp' : timestamp,
                       'shell' : shell,
                       'pipeline_name' : run_name,
                       'num_tasks' : num_bundles,
                       'queue' : "all.q",
                       'par_env' : "mpi_smp",
                       'cores_per_task' : self._num_subjects_per_bundle,
                       'user' : user_account,
                       'work_dir' : cluster_files_dir,
                       'plugin_args' : plugin_args}


        # Get string template for job scheduler
        if self._platform == "PBS":
            env_arr_idx = '$PBS_ARRAYID'
            batch_file_contents = cluster_templates.pbs_template
            confirm_str = '(?<=Your job-array )\d+'
            exec_cmd = 'qsub'
        elif self._platform == "SGE":
            env_arr_idx = '$SGE_TASK_ID'
            batch_file_contents = cluster_templates.sge_template
            confirm_str = '(?<=Your job-array )\d+'
            exec_cmd = 'qsub'
        elif self._platform == "SLURM":
            hrs_limit = 8*len(sublist)
            time_limit = '%d:00:00' % hrs_limit
            config_dict["time_limit"] = time_limit
            env_arr_idx = '$SLURM_ARRAY_TASK_ID'
            batch_file_contents = cluster_templates.slurm_template
            confirm_str = '(?<=Submitted batch job )\d+'
            exec_cmd = 'sbatch'


        config_dict['run_cmd'] = 'echo "Running task: %s"' % env_arr_idx

        # Populate string from config dict values
        batch_file_contents = batch_file_contents % config_dict


        if self._s3_dict_yml:
            subdict_arg = "--s3_dict_yml"
            subdict = self._s3_dict_yml
        elif self._config["subject_list"]:
            subdict_arg = "--sublist"
            subdict = self._config["subject_list"]

        if self._bundle_idx:
            idx_arg = "--bundle_idx"
            idx = self._bundle_idx
        elif self._subj_idx:
            idx_arg = "--subj_idx"
            idx = self._subj_idx

        run_str = "%s.py %s %s %s %d %s" % \
                      (self._config["qap_type"], \
                       subdict_arg, subdict, \
                       idx_arg, idx, \
                       self._config["pipeline_config_yaml"])

        batch_file_contents = "\n".join([batch_file_contents], run_str)

        batch_filepath = os.path.join(cluster_files_dir, 'cpac_submit_%s.%s' \
                                      % (timestamp, self._platform))


        return batch_file_contents, batch_filepath, exec_cmd, \
                   confirm_str, cluster_files_dir



    def _run_on_cluster(self, batch_file_contents, batch_filepath, exec_cmd, \
                            confirm_str, cluster_files_dir):

        with open(batch_filepath, 'w') as f:
            f.write(batch_file_contents)

        # Get output response from job submission
        out = commands.getoutput('%s %s' % (exec_cmd, batch_filepath))

        # Check for successful qsub submission
        if re.search(confirm_str, out) == None:
            err_msg = 'Error submitting QAP pipeline run to %s queue' \
                      % job_scheduler
            raise Exception(err_msg)

        # Get pid and send to pid file
        pid = re.search(confirm_str, out).group(0)
        pid_file = os.path.join(cluster_files_dir, 'pid.txt')
        with open(pid_file, 'w') as f:
            f.write(pid)



    def create_flat_sub_dict_dict(self, subdict):

        # "collapses" the subject resource pools so that each sub-session-scan
        # combination has its own entry

        # input
        #   subdict: resource pool dictionary loaded from a sublist YAML
        #   format:  
        #     {'sub_01': {'session_01': 
        #                    {'anatomical_scan': {'scan_01': <filepath>,
        #                                         'scan_02': <filepath>},
        #                     'site_name': 'Site_1'} },
        #      'sub_02': {..} }

        # output
        #   flat_sub_dict_dict is a dictionary of dictionaries. format:
        #     { (sub01,session01,scan01): {"anatomical_scan": <filepath>,
        #                                  "anatomical_brain": <filepath>} }

        flat_sub_dict_dict = {}
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
                            if sub_info_tuple not in flat_sub_dict_dict.keys():
                                flat_sub_dict_dict[sub_info_tuple] = {}

                            flat_sub_dict_dict[sub_info_tuple].update(resource_dict)

                    elif resource == "site_name":
                        sites_dict[subid] = subdict[subid][session][resource]

                    else:
                        filepath = subdict[subid][session][resource]
                        resource_dict = {}
                        resource_dict[resource] = filepath
                        sub_info_tuple = (subid, session, None)

                        if sub_info_tuple not in flat_sub_dict_dict.keys():
                            flat_sub_dict_dict[sub_info_tuple] = {}

                        flat_sub_dict_dict[sub_info_tuple].update(resource_dict)

        if len(flat_sub_dict_dict) == 0:
            # this error message meant more for devs than user
            msg = "The participant dictionary is empty."
            raise_smart_exception(locals(),msg)

        # in case some subjects have site names and others don't
        if len(sites_dict.keys()) > 0:
            for subid in subdict.keys():
                if subid not in sites_dict.keys():
                    sites_dict[subid] = None

            # integrate site information into flat_sub_dict_dict
            #     it was separate in the first place to circumvent the fact
            #     that even though site_name doesn't get keyed with scan names
            #     names, that doesn't necessarily mean scan names haven't been
            #     specified for that participant
            for sub_info_tuple in flat_sub_dict_dict.keys():
                site_info = {}
                site_info["site_name"] = sites_dict[sub_info_tuple[0]]
                flat_sub_dict_dict[sub_info_tuple].update(site_info)


        return flat_sub_dict_dict



    def _load_and_flatten_sublist(self):

        import yaml

        with open(self._config["subject_list"], "r") as f:
            subdict = yaml.load(f)

        # subdict is in the format:
        #   {'sub_01': {'session_01': 
        #                  {'anatomical_scan': {'scan_01': <filepath>,
        #                                       'scan_02': <filepath>},
        #                   'site_name': 'Site_1'} },
        #    'sub_02': {..} }

        if len(subdict) == 0:
            msg = "The participant list provided is either empty or could " \
                  "not be read properly!"
            raise_smart_exception(locals(),msg)


        # Generate flat_sub_dict_dict

        # flat_sub_dict_dict is a dictionary of dictionaries. format:
        #   { (sub01,session01,scan01): {"anatomical_scan": <filepath>,
        #                                "anatomical_brain": <filepath>} }
        flat_sub_dict_dict = self.create_flat_sub_dict_dict(subdict)


        return flat_sub_dict_dict



    def _create_bundles(self, flat_sub_dict_dict):

        # input
        #   flat_sub_dict_dict is a dictionary of dictionaries. format:
        #     { (sub01,session01,scan01): {"anatomical_scan": <filepath>,
        #                                  "anatomical_brain": <filepath>} }

        # output
        #   bundles: a list of bundles - each bundle being a dictionary that
        #            is a starting resource pool for N sub-session-scan combos
        #            with N being the number of subjects per bundle
        #            (set by the user)

        i = 0
        bundles = []

        if len(flat_sub_dict_dict) < self._num_subjects_per_bundle:
            bundles.append(flat_sub_dict_dict)
        else:
            for sub_info_tuple in flat_sub_dict_dict.keys():
                if i == 0:
                    new_bundle = {}
                new_bundle[sub_info_tuple] = flat_sub_dict_dict[sub_info_tuple]
                i += 1
                if i == self._num_subjects_per_bundle:
                    bundles.append(new_bundle)
                    i = 0

        if len(bundles) == 0:
            msg = "No bundles created."
            raise_smart_exception(locals(),msg)


        return bundles



    def _run_here(self, run_name):

        # input
        #   run_name: the filename of the pipeline config YAML file

        # get flattened sublist
        flat_sub_dict_dict = self._load_and_flatten_sublist()

        logger.info('There are %d subjects in the pool' %
                    len(flat_sub_dict_dict.keys()))

        # Create bundles

        # bundles is a list of "bundles" - each bundle being a dictionary that
        # is a starting resource pool for N sub-session-scan combos with N
        # being the number of subjects per bundle (set by the user)
        bundles = self._create_bundles(flat_sub_dict_dict)

        # Stack workflow args, make a list of tuples containing run args
        #     one tuple for each bundle
        #     len(wfargs) = number of bundles
        wfargs = [(data_bundle, data_bundle.keys(), self._config, run_name,
                   self.runargs) for data_bundle in bundles]

        results = []

        # skip parallel machinery if we are running only one bundle at once
        if self._num_bundles_at_once == 1:
            for a in wfargs:
                results.append(_run_workflow(a))
        # or use Pool if running multiple bundles simultaneously
        else:
            from multiprocessing import Pool

            try:
                pool = Pool(processes=self._num_bundles_at_once, \
                            masktasksperchild=50)
            except TypeError:  # Make python <2.7 compatible
                pool = Pool(processes=self._num_bundles_at_once)

            results = pool.map(_run_workflow, wfargs)
            pool.close()
            pool.terminate()


        return results



    def _run_one_bundle(self, run_name):

        # kick off a single-bundle run
        #   this will be called multiple times throughout the execution of
        #   a batch script for SGE on the cloud

        from cloud_utils import dl_subj_from_s3, upl_qap_output

        self._sub_dict = {}

        # if the user is using S3 storage, download the bundle or subject
        if self._s3_dict_yml:

            # s3_dict_yml is a dictionary of dictionaries, keyed with
            # sub-session-scan tuples, and is generated by the 
            # qap_aws_s3_dict_generator.py script, containing filepath info of
            # data on S3 storage
            #   format:
            #     { (sub01,session01,scan01):
            #           {"anatomical_scan": <S3 filepath>},
            #       (sub02,session01,scan01):
            #           {"anatomical_scan": <S3 filepath>}, ..}

            # bundle_idx is an integer (minimum of 1) which will be used to
            # calculate where in the s3_dict_yml to begin downloading subjects
            # for that bundle, based on the number of subjects per bundle
            # (specified by the user)

            # subj_idx is an integer (minimum of 1) which will select the
            # entry from s3_dict_yml to retrieve using the dl_subj_from_s3
            # function

            if self._bundle_idx:

                # download data from S3 by the bundle
                #   ex. if we are doing 4 subjects per bundle and we are on
                #       our 3rd bundle (bundle_idx = 3), then this should be
                #       subjects 9-12 from the s3_dict_yml

                first_subj_idx = \
                    (bundle_idx-1)*self._config["num_subjects_per_bundle"] + 1

                last_subj_idx = \
                    first_subj_idx + self._config["num_subjects_per_bundle"]

                for idx in range(first_subj_idx, last_subj_idx):
                    single_sub_dict = dl_subj_from_s3(idx, self._config, \
                                                          self._s3_dict_yml)
                    self._sub_dict.update(single_sub_dict)

            elif self._subj_idx:

                self._sub_dict = dl_subj_from_s3(self._subj_idx, \
                    self._config, self._s3_dict_yml)


            if len(self._sub_dict) == 0:
                err = "\n[!] Subject dictionary was not successfully " \
                      "downloaded from the S3 bucket!\n"
                raise RuntimeError(err)


            try:
                # integrate site information into the subject list
                #   it was separate in the first place to circumvent the fact
                #   that even though site_name doesn't get keyed with scan
                #   names, that doesn't necessarily mean scan names haven't
                #   been specified for that participant

                for sub in self._sub_dict.keys():
                    for resource_path in subject_list[sub].values():
                        if ".nii" in resource_path:
                            filepath = resource_path
                            break

                    filesplit = filepath.split(self._config["bucket_prefix"])
                    site_name = filesplit[1].split("/")[1]
                    self._sub_dict[sub]["site_name"] = site_name

            except:
                pass


        elif self._config["subject_list"]:

            # run on a cluster without pulling data from S3

            # get flat sublist
            flat_sub_dict_dict = self._load_and_flatten_sublist()

            if self._bundle_idx:

                first_subj_idx = \
                    (bundle_idx-1)*self._config["num_subjects_per_bundle"] + 1

                last_subj_idx = \
                    first_subj_idx + self._config["num_subjects_per_bundle"]

                for idx in range(first_subj_idx, last_subj_idx):

                    # Get list of subject keys for indexing
                    sd_keys = flat_sub_dict_dict.keys()
                    sd_keys.sort()

                    # Grab subject dictionary of interest
                    subj_key = sd_keys[idx-1]
                    single_sub_dict = flat_sub_dict_dict[subj_key]

                    self._sub_dict.update(single_sub_dict)

            elif self._subj_idx:

                # Get list of subject keys for indexing
                sd_keys = flat_sub_dict_dict.keys()
                sd_keys.sort()

                # Grab subject dictionary of interest
                subj_key = sd_keys[self._subj_idx-1]
                self._sub_dict = flat_sub_dict_dict[subj_key]
        

        # let's go!
        rt = _run_workflow((self._sub_dict, self._sub_dict.keys(), \
                               self._config, run_name, self.runargs))
    
        # make not uploading results to S3 bucket the default if not specified
        if "upload_to_s3" not in self._config.keys():
            self._config["upload_to_s3"] = False

        # upload results
        if self._config["upload_to_s3"]:
            upl_qap_output(self._config)


        return rt



    def run(self):

        from qap.workflow_utils import raise_smart_exception

        # Get configurations and settings
        config = self._config
        self._num_subjects_per_bundle = config.get('num_subjects_per_bundle', 1)
        self._num_bundles_at_once = int(self._config["num_bundles_at_once"])
        write_report = config.get('write_report', False)

        if "resource_manager" in config.keys():
            platforms = ["SGE","PBS","SLURM"]
            self._platform = str(config["resource_manager"]).upper()
            if self._platform not in platforms:
                msg = "The resource manager %s provided in the pipeline " \
                      "configuration file is not one of the valid choices. " \
                      "It must be one of the following:\n%s" \
                      % (self._platform,str(platforms))
                raise_smart_exception(locals(),msg)
        else:
            self._platform = None

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


        # settle run arguments (plugins)
        self.runargs = {'plugin': 'Linear', 'plugin_args': {}}
        if self._num_subjects_per_bundle > 1:
            self.runargs['plugin'] = 'MultiProc'
            self.runargs['plugin_args'] = \
                {'n_procs': self._num_subjects_per_bundle}


        # Start the magic
        if self._cloudify:

            results = self._run_one_bundle()

        elif not self._platform:

            results = self._run_here(run_name)

        elif self._platform:

            import yaml
            import math

            # get num_bundles
            bundle_size = self._config["num_subjects_per_bundle"]

            if self._s3_dict_yml:

                with open(self._s3_dict_yml,"r") as f:
                    s3_dict = yaml.load(f)

                num_bundles = len(s3_dict) / bundle_size

            elif self._config["subject_list"]:

                # get flattened sublist
                flat_sub_dict_dict = self._load_and_flatten_sublist()

                num_bundles = \
                    len(flat_sub_dict_dict) / bundle_size


            # round up if it is a float
            num_bundles = int(math.ciel(num_bundles))

            batch_file_contents, batch_filepath, exec_cmd, \
                   confirm_str, cluster_files_dir = \
                   self._prepare_cluster_batch_file(run_name, num_bundles)

            self._run_on_cluster(batch_file_contents, batch_filepath, \
                                     exec_cmd, confirm_str, cluster_files_dir)


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



def _run_workflow(args):

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

    # unpack args
    resource_pool_list, sub_info_list, config, run_name, runargs = args

    qap_type = config['qap_type']


    # Read and apply general settings in config
    keep_outputs = config.get('write_all_outputs', False)

    num_subjects_at_once = config.get('num_subjects_at_once',1)

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

    workflow.base_dir = op.join(config["working_directory"])

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


        name = "_" + sub_id + "_" + session_id + "_" + scan_id


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
            workflow, resource_pool = wf_builder(workflow, resource_pool, \
                                                 config, name)


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
            if (len(resource_pool[output]) == 2) and (output != "starter"):
                ds = pe.Node(nio.DataSink(), name='datasink_%s%s' \
                    % (output,name))
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
                dotfilename=op.join(config["output_directory"], \
                                    run_name + ".dot"),
                simple_form=False)
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
