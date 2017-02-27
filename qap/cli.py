#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import os.path as op
import argparse

from nipype import config
log_dir=os.path.join("tmp","nipype","logs")
config.update_config({'logging': {'log_directory': log_dir,
                                  'log_to_file': True}})

from nipype import logging
logger = logging.getLogger('workflow')


class QAProtocolCLI:
    """
    This class and the associated run_workflow function implement what
    the former scripts (qap_anatomical_spatial.py, etc.) contained
    """

    def __init__(self, parse_args=True):

        if parse_args:
            self._parse_args()
        else:
            self._cloudify = False
            self._bundle_idx = None

    def _parse_args(self):

        parser = argparse.ArgumentParser()

        group = parser.add_argument_group(
            "Regular Use Inputs (non-cloud runs)")
        cloudgroup = parser.add_argument_group(
            "AWS Cloud Inputs (only required for AWS Cloud runs)")
        req = parser.add_argument_group("Required Inputs")

        cloudgroup.add_argument('--bundle_idx', type=int,
                                help='Bundle index to run')
        cloudgroup.add_argument('--log_dir', type=str,
                                help='Directory for workflow logging')

        # Subject list (YAML file)
        group.add_argument(
            "data_config_file", type=str,
            help="filepath to participant list YAML")
        req.add_argument(
            "pipeline_config_file", type=str,
            help="filepath to pipeline configuration YAML")

        # Write PDF reports
        group.add_argument(
            "--with-reports", action='store_true', default=False,
            help="Write a summary report in PDF format.")

        args = parser.parse_args()

        # Load config
        from qap.script_utils import read_yml_file
        self._config = read_yml_file(args.pipeline_config_file)
        self.validate_config_dict()

        self._config['pipeline_config_yaml'] = \
            os.path.realpath(args.pipeline_config_file)
        self._run_name = self._config['pipeline_name']

        if args.with_reports:
            self._config['write_report'] = True

        if "num_sessions_at_once" not in self._config.keys():
            self._config["num_sessions_at_once"] = 1

        if "num_bundles_at_once" not in self._config.keys():
            self._config["num_bundles_at_once"] = 1

        self._config["subject_list"] = os.path.realpath(args.data_config_file)

        if args.bundle_idx:
            self._bundle_idx = args.bundle_idx
        else:
            self._bundle_idx = None

        if args.log_dir:
            self._run_log_dir = args.log_dir
        else:
            self._run_log_dir = self._config["output_directory"]

    def submit_cluster_batch_file(self, num_bundles):
        """Write the cluster batch file for the appropriate scheduler.

        - Batch file setup code borrowed from dclark87's CPAC cluster setup
          code:
              - https://github.com/FCP-INDI/C-PAC/blob/0.4.0_development/CPAC/pipeline/cpac_runner.py
              - https://github.com/dclark87
        - This function will write the batch file appropriate for the
          scheduler being used, and then this CLI will be run again on each
          node/slot through the run_one_bundle function.

        :type num_bundles: int
        :param num_bundles: The number of bundles total being run.
        """

        import os
        import re
        import getpass
        import commands
        from time import strftime
        from indi_schedulers import cluster_templates

        print "Submitting cluster job to %s.." % self._platform

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
        config_dict = {'timestamp': timestamp,
                       'shell': shell,
                       'job_name': self._run_name,
                       'num_tasks': num_bundles,
                       'queue': "all.q",
                       'par_env': "mpi_smp",
                       'cores_per_task': self._config["num_processors"],
                       'user': user_account,
                       'work_dir': cluster_files_dir}

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
            hrs_limit = 8 * num_bundles
            time_limit = '%d:00:00' % hrs_limit
            config_dict["time_limit"] = time_limit
            env_arr_idx = '$SLURM_ARRAY_TASK_ID'
            batch_file_contents = cluster_templates.slurm_template
            confirm_str = '(?<=Submitted batch job )\d+'
            exec_cmd = 'sbatch'

        config_dict['env_arr_idx'] = env_arr_idx
        config_dict['run_cmd'] = 'echo "Running task: %s"' % env_arr_idx

        # Populate string from config dict values
        batch_file_contents = batch_file_contents % config_dict

        run_str = "qap_measures_pipeline.py --bundle_idx %s --log_dir %s %s "\
                  "%s" % (env_arr_idx, self._run_log_dir,
                          self._config["subject_list"],
                          self._config["pipeline_config_yaml"])

        batch_file_contents = "\n".join([batch_file_contents, run_str])

        batch_filepath = os.path.join(cluster_files_dir, 'cpac_submit_%s.%s'
                                      % (timestamp, self._platform))

        with open(batch_filepath, 'w') as f:
            f.write(batch_file_contents)

        print "Batch file written to %s.." % batch_filepath

        # Get output response from job submission
        out = commands.getoutput('%s %s' % (exec_cmd, batch_filepath))

        # Check for successful qsub submission
        if re.search(confirm_str, out) == None:
            err_msg = 'Error submitting QAP pipeline run to %s queue' \
                      % self._platform
            raise Exception(err_msg)

        print "Batch job submitted to %s queue." % self._platform

        # Get pid and send to pid file
        pid = re.search(confirm_str, out).group(0)
        pid_file = os.path.join(cluster_files_dir, 'pid.txt')
        with open(pid_file, 'w') as f:
            f.write(pid)

    def validate_config_dict(self):
        """Validate the pipeline configuration dictionary to ensure the
        parameters are properly set.
        """
        config_options = ["pipeline_name",
                          "num_processors",
                          "num_sessions_at_once",
                          "available_memory",
                          "cluster_system",
                          "output_directory",
                          "working_directory",
                          "template_head_for_anat",
                          "exclude_zeros",
                          "start_idx",
                          "stop_idx",
                          "write_report",
                          "write_graph",
                          "write_all_outputs",
                          "upload_to_s3",
                          "bucket_prefix",
                          "bucket_out_prefix",
                          "local_prefix",
                          "bucket_name",
                          "creds_path"]
        invalid = []
        for param in self._config.keys():
            if param not in config_options:
                invalid.append(param)
        if len(invalid) > 0:
            err = "\n[!] The following parameters in your configuration " \
                  "file are not recognized. Double-check the pipeline " \
                  "configuration template.\n"
            err += "\n".join([x for x in invalid])
            raise Exception(err)
        else:
            return 0

    def create_session_dict(self, subdict):
        """Collapse the participant resource pools so that each participant-
        session combination has its own entry.

        - input subdict format:
              {'sub_01': {'session_01':
                             {'anatomical_scan': {'scan_01': <filepath>,
                                                  'scan_02': <filepath>},
                              'site_name': 'Site_1'} },
               'sub_02': {..} }

        - output dict format:
              { (sub01,session01): {"scan_01": {
                                            "anatomical_scan": <filepath>},
                                   {"scan_02": {
                                            "anatomical_scan": <filepath>} } }

        :type subdict: dict
        :param subdict: A dictionary containing the filepaths of input files
                        for each participant, sorted by session and scan.
        :rtype: dict
        :return: A dictionary of dictionaries where each participant-session
                 combination has its own entry, and input file filepaths are
                 defined.
        """

        from qap.qap_utils import raise_smart_exception

        flat_sub_dict_dict = {}
        sites_dict = {}

        for subid in subdict.keys():
            subid = str(subid)
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
                            sub_info_tuple = (subid, session)
                            if sub_info_tuple not in flat_sub_dict_dict.keys():
                                flat_sub_dict_dict[sub_info_tuple] = {}
                            if scan not in flat_sub_dict_dict[sub_info_tuple].keys():
                                flat_sub_dict_dict[sub_info_tuple][scan] = {}

                            flat_sub_dict_dict[sub_info_tuple][scan].update(
                                resource_dict)

                    elif resource == "site_name":
                        sites_dict[subid] = subdict[subid][session][resource]

                    else:
                        filepath = subdict[subid][session][resource]
                        resource_dict = {}
                        resource_dict[resource] = filepath
                        sub_info_tuple = (subid, session)

                        if sub_info_tuple not in flat_sub_dict_dict.keys():
                            flat_sub_dict_dict[sub_info_tuple] = {}

                        flat_sub_dict_dict[sub_info_tuple].update(
                            resource_dict)

        if len(flat_sub_dict_dict) == 0:
            # this error message meant more for devs than user
            msg = "The participant dictionary is empty."
            raise_smart_exception(locals(), msg)

        # in case some subjects have site names and others don't
        if len(sites_dict.keys()) > 0:
            for subid in subdict.keys():
                subid = str(subid)
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

    def load_sublist(self):
        """Load the participant list YAML file into a dictionary and check.

        - subdict format:
              {'sub_01': {'session_01':
                            {'anatomical_scan': {'scan_01': <filepath>,
                                                 'scan_02': <filepath>},
                             'site_name': 'Site_1'} },
              'sub_02': {..} }

        :rtype: dictionary
        :return: The participant list in a dictionary.
        """

        import yaml
        from qap.qap_utils import raise_smart_exception

        if "subject_list" in self._config.keys():
            with open(self._config["subject_list"], "r") as f:
                subdict = yaml.load(f)
        else:
            msg = "\n\n[!] There is no participant list YML to read.\n\n"
            raise_smart_exception(locals(),msg)

        if len(subdict) == 0:
            msg = "The participant list provided is either empty or could " \
                  "not be read properly!"
            raise_smart_exception(locals(),msg)

        return subdict

    def create_bundles(self):
        """Create a list of participant "bundles".

        :rtype: list
        :return: A list of bundles - each bundle being a dictionary that is a
                 starting resource pool for N sub-session-scan combos with N
                 being the number of participants per bundle (set by the user)
        """

        from qap.qap_utils import raise_smart_exception

        i = 0
        bundles = []

        for session_tuple in self._sub_dict.keys():
            if i == 0:
                new_bundle = {}
            sub = session_tuple[0]
            ses = session_tuple[1]
            site_name = None
            if "site_name" in self._sub_dict[session_tuple].keys():
                site_name = self._sub_dict[session_tuple]["site_name"]
            for scan in self._sub_dict[session_tuple].keys():
                if type(self._sub_dict[session_tuple][scan]) is dict:
                    # to avoid fields in sub_dict[session_tuple] that are
                    # strings (such as site_name or creds_path)
                    sub_info_tuple = (sub, ses, scan)
                    new_bundle[sub_info_tuple] = \
                        self._sub_dict[session_tuple][scan]
                    if site_name:
                        new_bundle[sub_info_tuple].update({"site_name": site_name})
            i += 1
            if i == self._config["num_sessions_at_once"]:
                bundles.append(new_bundle)
                i = 0

        if i > 0:
            bundles.append(new_bundle)

        if len(bundles) == 0:
            msg = "No bundles created."
            raise_smart_exception(locals(),msg)

        return bundles

    def run_one_bundle(self, bundle_idx, run=True):
        """Execute one bundle's workflow on one node/slot of a cluster/grid.

        - Compatible with Amazon AWS cluster runs, and S3 buckets.

        :type bundle_idx: int
        :param bundle_idx: The bundle ID number - used to calculate which
                           entries in the participant list to pull into the
                           current bundle, based on the number of participants
                           per bundle (participants at once).
        :type run: bool
        :param run: Run flag, set to False for testing.
        :rtype: dictionary
        :return: A dictionary with information about the workflow run, its
                  status, and results.
        """

        import os
        from qap_utils import write_json
        from cloud_utils import download_single_s3_path

        self._config["workflow_log_dir"] = self._run_log_dir

        bundle_dict = self._bundles_list[bundle_idx-1]
        num_bundles = len(self._bundles_list)

        # check for s3 paths
        for sub in bundle_dict.keys():
            # in case we're dealing with string entries in the data dict
            try:
                bundle_dict[sub].keys()
            except AttributeError:
                continue
            for resource in bundle_dict[sub].keys():
                value = bundle_dict[sub][resource]
                if "s3://" in value:
                    bundle_dict[sub][resource] = \
                        download_single_s3_path(value, self._config)

        wfargs = (bundle_dict, bundle_dict.keys(),
                  self._config, self._run_name, self.runargs,
                  bundle_idx, num_bundles)

        if run:
            # let's go!
            rt = run_workflow(wfargs)

            # write bundle results to JSON file
            write_json(rt, os.path.join(rt["bundle_log_dir"],
                                        "workflow_results.json"))

            # make not uploading results to S3 bucket the default if not
            # specified
            if "upload_to_s3" not in self._config.keys():
                self._config["upload_to_s3"] = False

            # upload results
            if self._config["upload_to_s3"]:
                from cloud_utils import upl_qap_output
                upl_qap_output(self._config)

            return rt
        else:
            return wfargs

    def run(self, config_file=None, partic_list=None):
        """Establish where and how we're running the pipeline and set up the
        run. (Entry point)

        - This is the entry point for pipeline building and connecting.
          Depending on the inputs, the appropriate workflow runner will
          be selected and executed.

        :type config_file: str
        :param config_file: Filepath to the pipeline configuration file in
                            YAML format.
        :type partic_list: str
        :param partic_list: Filepath to the participant list file in YAML
                            format.
        """

        from time import strftime
        from qap_utils import raise_smart_exception, \
                              check_config_settings

        # in case we are overloading
        if config_file:
            from qap.script_utils import read_yml_file
            self._config = read_yml_file(config_file)
            self.validate_config_dict()
            self._config["pipeline_config_yaml"] = config_file
      
        if not self._config:
             raise Exception("config not found!")

        if partic_list:
            self._config["subject_list"] = partic_list

        # Get configurations and settings
        check_config_settings(self._config, "num_processors")
        check_config_settings(self._config, "num_sessions_at_once")
        check_config_settings(self._config, "available_memory")
        check_config_settings(self._config, "output_directory")
        check_config_settings(self._config, "working_directory")

        self._num_bundles_at_once = 1
        write_report = self._config.get('write_report', False)

        if "cluster_system" in self._config.keys() and not self._bundle_idx:
            res_mngr = self._config["cluster_system"]
            if (res_mngr == None) or ("None" in res_mngr) or \
                ("none" in res_mngr):
                self._platform = None
            else:
                platforms = ["SGE", "PBS", "SLURM"]
                self._platform = str(res_mngr).upper()
                if self._platform not in platforms:
                    msg = "The resource manager %s provided in the pipeline "\
                          "configuration file is not one of the valid " \
                          "choices. It must be one of the following:\n%s" \
                          % (self._platform, str(platforms))
                    raise_smart_exception(locals(), msg)
        else:
            self._platform = None

        # Create output directory
        try:
            os.makedirs(self._config["output_directory"])
        except:
            if not op.isdir(self._config["output_directory"]):
                err = "[!] Output directory unable to be created.\n" \
                      "Path: %s\n\n" % self._config["output_directory"]
                raise Exception(err)
            else:
                pass

        # Create working directory
        try:
            os.makedirs(self._config["working_directory"])
        except:
            if not op.isdir(self._config["working_directory"]):
                err = "[!] Output directory unable to be created.\n" \
                      "Path: %s\n\n" % self._config["working_directory"]
                raise Exception(err)
            else:
                pass

        results = []

        # set up callback logging
        import logging
        from nipype.pipeline.plugins.callback_log import log_nodes_cb

        cb_log_filename = os.path.join(self._config["output_directory"],
                                       "callback.log")
        # Add handler to callback log file
        cb_logger = logging.getLogger('callback')
        cb_logger.setLevel(logging.DEBUG)
        handler = logging.FileHandler(cb_log_filename)
        cb_logger.addHandler(handler)

        # settle run arguments (plugins)
        self.runargs = {}
        self.runargs['plugin'] = 'MultiProc'
        self.runargs['plugin_args'] = \
            {'memory_gb': int(self._config["available_memory"]),
             'status_callback': log_nodes_cb}
        n_procs = {'n_procs': self._config["num_processors"]}
        self.runargs['plugin_args'].update(n_procs)

        # load the participant list file into dictionary
        subdict = self.load_sublist()

        # flatten the participant dictionary
        self._sub_dict = self.create_session_dict(subdict)

        # create the list of bundles
        self._bundles_list = self.create_bundles()
        num_bundles = len(self._bundles_list)

        if not self._bundle_idx:
            # want to initialize the run-level log directory (not the bundle-
            # level) only the first time we run the script, due to the
            # timestamp. if sub-nodes are being kicked off by a batch file on
            # a cluster, we don't want a new timestamp for every new node run
            self._run_log_dir = op.join(self._config['output_directory'],
                                        '_'.join([self._run_name, "logs"]),
                                        '_'.join([strftime("%Y%m%d_%H_%M_%S"),
                                                 "%dbundles" % num_bundles]))

        if self._run_log_dir:
            if not os.path.isdir(self._run_log_dir):
                try:
                    os.makedirs(self._run_log_dir)
                except:
                    if not op.isdir(self._run_log_dir):
                        err = "[!] Log directory unable to be created.\n" \
                              "Path: %s\n\n" % self._run_log_dir
                        raise Exception(err)
                    else:
                        pass

        if num_bundles == 1:
            self._config["num_sessions_at_once"] = \
                len(self._bundles_list[0])

        # Start the magic
        if not self._platform and not self._bundle_idx:
            # not a cluster/grid run
            for idx in range(1, num_bundles+1):
                results.append(self.run_one_bundle(idx))

        elif not self._bundle_idx:
            # there is a self._bundle_idx only if the pipeline runner is run
            # with bundle_idx as a parameter - only happening either manually,
            # or when running on a cluster
            self.submit_cluster_batch_file(num_bundles)

        else:
            # if there is a bundle_idx supplied to the runner
            results = self.run_one_bundle(self._bundle_idx)


def starter_node_func(starter):
    """Pass a dummy string through to provide a basic function for the first
    Nipype workflow node.

    - This is used for a Nipype utility function node to serve as a starting
      node to connect to multiple unrelated Nipype workflows. Each of these
      workflows runs QAP for one participant in the current bundle being run.
    - Connecting the multiple non-interdependent participant workflows as
      one workflow allows the Nipype resource scheduler to maximize
      performance.

    :type starter: str
    :param starter: A dummy string.
    :rtype: str
    :return: The same string.
    """
    return starter


def run_workflow(args, run=True):
    """Connect and execute the QAP Nipype workflow for one bundle of data.

    - This function will update the resource pool with what is found in the
      output directory (if it already exists). If the final expected output
      of the pipeline is already found, the pipeline will not run and it
      will move onto the next bundle. If the final expected output is not
      present, the pipeline begins to build itself backwards.

    :type args: tuple
    :param args: A 7-element tuple of information comprising of the bundle's
                 resource pool, a list of participant info, the configuration
                 options, the pipeline ID run name and miscellaneous run args.
    :rtype: dictionary
    :return: A dictionary with information about the workflow run, its status,
             and results.
    """

    import os
    import os.path as op

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as niu

    import qap
    from qap_utils import read_json

    import glob

    import time
    from time import strftime
    from nipype import config as nyconfig

    # unpack args
    resource_pool_dict, sub_info_list, config, run_name, runargs, \
        bundle_idx, num_bundles = args

    # Read and apply general settings in config
    keep_outputs = config.get('write_all_outputs', False)

    # take date+time stamp for run identification purposes
    pipeline_start_stamp = strftime("%Y-%m-%d_%H:%M:%S")
    pipeline_start_time = time.time()

    if "workflow_log_dir" not in config.keys():
        config["workflow_log_dir"] = config["output_directory"]

    bundle_log_dir = op.join(config["workflow_log_dir"],
                             '_'.join(["bundle", str(bundle_idx)]))

    try:
        os.makedirs(bundle_log_dir)
    except:
        if not op.isdir(bundle_log_dir):
            err = "[!] Bundle log directory unable to be created.\n" \
                    "Path: %s\n\n" % bundle_log_dir
            raise Exception(err)
        else:
            pass

    # set up logging
    nyconfig.update_config(
        {'logging': {'log_directory': bundle_log_dir, 'log_to_file': True}})
    logging.update_logging(nyconfig)

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

    # iterate over each subject in the bundle
    logger.info("Starting bundle %s out of %s.." % (str(bundle_idx),
                                                    str(num_bundles)))
    # results dict
    rt = {'status': 'Started', 'bundle_log_dir': bundle_log_dir}

    for sub_info in sub_info_list:

        resource_pool = resource_pool_dict[sub_info]

        # in case we're dealing with string entries in the data dict
        try:
            resource_pool.keys()
        except AttributeError:
            continue

        # resource pool check
        invalid_paths = []

        for resource in resource_pool.keys():
            try:
                if not op.isfile(resource_pool[resource]) and resource != "site_name":
                    invalid_paths.append((resource, resource_pool[resource]))
            except:
                err = "\n\n[!]"
                raise Exception(err)

        if len(invalid_paths) > 0:
            err = "\n\n[!] The paths provided in the subject list to the " \
                  "following resources are not valid:\n"

            for path_tuple in invalid_paths:
                err = "%s%s: %s\n" % (err, path_tuple[0], path_tuple[1])

            err = "%s\n\n" % err
            raise Exception(err)

        # process subject info
        sub_id = str(sub_info[0])
        # for nipype
        if "-" in sub_id:
            sub_id = sub_id.replace("-","_")
        if "." in sub_id:
            sub_id = sub_id.replace(".","_")

        if sub_info[1]:
            session_id = str(sub_info[1])
            # for nipype
            if "-" in session_id:
                session_id = session_id.replace("-","_")
            if "." in session_id:
                session_id = session_id.replace(".","_")
        else:
            session_id = "session_0"

        if sub_info[2]:
            scan_id = str(sub_info[2])
            # for nipype
            if "-" in scan_id:
                scan_id = scan_id.replace("-","_")
            if "." in scan_id:
                scan_id = scan_id.replace(".","_")
        else:
            scan_id = "scan_0"

        name = "_".join(["", sub_id, session_id, scan_id])

        rt[name] = {'id': sub_id, 'session': session_id, 'scan': scan_id,
                    'resource_pool': str(resource_pool)}

        logger.info("Participant info: %s" % name)

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

        # for QAP spreadsheet generation only
        config.update({"subject_id": sub_id, "session_id": session_id,
                       "scan_id": scan_id, "run_name": run_name})

        if "site_name" in resource_pool:
            config.update({"site_name": resource_pool["site_name"]})

        logger.info("Configuration settings:\n%s" % str(config))

        qap_types = ["anatomical_spatial", 
                     "functional_spatial", 
                     "functional_temporal"]

        # update that resource pool with what's already in the output
        # directory
        for resource in os.listdir(output_dir):
            if (op.exists(op.join(output_dir, resource)) and
                    resource not in resource_pool.keys()):
                try:
                    resource_pool[resource] = \
                        glob.glob(op.join(output_dir, resource, "*"))[0]
                except IndexError:
                    if ".json" in resource:
                        # load relevant json info into resource pool
                        json_file = op.join(output_dir, resource)
                        json_dict = read_json(json_file)
                        sub_json_dict = json_dict["%s %s %s" % (sub_id,
                                                                session_id,
                                                                scan_id)]

                        if "anatomical_header_info" in sub_json_dict.keys():
                            resource_pool["anatomical_header_info"] = \
                                sub_json_dict["anatomical_header_info"]

                        if "functional_header_info" in sub_json_dict.keys():
                            resource_pool["functional_header_info"] = \
                                sub_json_dict["functional_header_info"]

                        for qap_type in qap_types:
                            if qap_type in sub_json_dict.keys():
                                resource_pool["_".join(["qap",qap_type])] = \
                                    sub_json_dict[qap_type]
                except:
                    # a stray file in the sub-sess-scan output directory
                    pass

        # create starter node which links all of the parallel workflows within
        # the bundle together as a Nipype pipeline
        resource_pool["starter"] = (starter_node, 'starter')

        # individual workflow and logger setup
        logger.info("Contents of resource pool for this participant:\n%s"
                    % str(resource_pool))

        # start connecting the pipeline
        qw = None
        for qap_type in qap_types:
            if "_".join(["qap", qap_type]) not in resource_pool.keys():
                if qw is None:
                    from qap import qap_workflows as qw
                wf_builder = \
                    getattr(qw, "_".join(["qap", qap_type, "workflow"]))
                workflow, resource_pool = wf_builder(workflow, resource_pool,
                                                     config, name)

        if ("anatomical_scan" in resource_pool.keys()) and \
            ("anatomical_header_info" not in resource_pool.keys()):
            if qw is None:
                from qap import qap_workflows as qw
            workflow, resource_pool = \
                qw.qap_gather_header_info(workflow, resource_pool, config,
                    name, "anatomical")

        if ("functional_scan" in resource_pool.keys()) and \
            ("functional_header_info" not in resource_pool.keys()):
            if qw is None:
                from qap import qap_workflows as qw
            workflow, resource_pool = \
                qw.qap_gather_header_info(workflow, resource_pool, config,
                    name, "functional")

        # set up the datasinks
        out_list = []
        for output in resource_pool.keys():
            for qap_type in qap_types:
                if qap_type in output:
                    out_list.append("_".join(["qap", qap_type]))

        # write_all_outputs (writes everything to the output directory, not
        # just the final JSON files)
        if keep_outputs:
            out_list = resource_pool.keys()
        logger.info("Outputs we're keeping: %s" % str(out_list))
        logger.info('Resource pool keys after workflow connection: '
                    '{}'.format(str(resource_pool.keys())))

        # Save reports to out_dir if necessary
        if config.get('write_report', False):

            if ("qap_mosaic" in resource_pool.keys()) and  \
                    ("qap_mosaic" not in out_list):
                out_list += ['qap_mosaic']

            # The functional temporal also has an FD plot
            if 'qap_functional_temporal' in resource_pool.keys():
                if ("qap_fd" in resource_pool.keys()) and \
                        ("qap_fd" not in out_list):
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
                ds = pe.Node(nio.DataSink(), name='datasink_%s%s'
                                                  % (output,name))
                ds.inputs.base_directory = output_dir
                node, out_file = resource_pool[output]
                workflow.connect(node, out_file, ds, output)
                new_outputs += 1
            elif ".json" in resource_pool[output]:
                new_outputs += 1

    logger.info("New outputs: %s" % str(new_outputs))

    # run the pipeline (if there is anything to do)
    if new_outputs > 0:
        if config.get('write_graph', False):
            workflow.write_graph(
                dotfilename=op.join(config["output_directory"],
                                    "".join([run_name, ".dot"])),
                simple_form=False)
            workflow.write_graph(
                graph2use="orig",
                dotfilename=op.join(config["output_directory"],
                                    "".join([run_name, ".dot"])),
                simple_form=False)
            workflow.write_graph(
                graph2use="hierarchical",
                dotfilename=op.join(config["output_directory"],
                                    "".join([run_name, ".dot"])),
                simple_form=False)
        if run:
            try:
                logger.info("Running with plugin %s" % runargs["plugin"])
                logger.info("Using plugin args %s" % runargs["plugin_args"])
                workflow.run(plugin=runargs["plugin"],
                             plugin_args=runargs["plugin_args"])
                rt['status'] = 'finished'
                logger.info("Workflow run finished for bundle %s."
                            % str(bundle_idx))
            except Exception as e:  # TODO We should be more specific here ...
                errmsg = e
                rt.update({'status': 'failed'})
                logger.info("Workflow run failed for bundle %s."
                            % str(bundle_idx))
                # ... however this is run inside a pool.map: do not raise
                # Exception
        else:
            return workflow

    else:
        rt['status'] = 'cached'
        logger.info("\nEverything is already done for bundle %s."
                    % str(bundle_idx))

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

    if rt["status"] == "failed":
        logger.error(errmsg)
    else:
        pipeline_end_stamp = strftime("%Y-%m-%d_%H:%M:%S")
        pipeline_end_time = time.time()
        logger.info("Elapsed time (minutes) since last start: %s"
                    % ((pipeline_end_time - pipeline_start_time) / 60))
        logger.info("Pipeline end time: %s" % pipeline_end_stamp)

    return rt
