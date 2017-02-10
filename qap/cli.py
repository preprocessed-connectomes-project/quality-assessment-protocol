#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import os
import os.path as op
import argparse
from qap_pipeline import build_and_run_qap_pipeline

# from nipype import config
# log_dir=os.path.join("tmp","nipype","logs")
# config.update_config({'logging': {'log_directory': log_dir, 'log_to_file': True}})
#
# from nipype import logging
# logger = logging.getLogger('workflow')



def submit_cluster_batch_file( num_bundles):
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

    print "Submitting cluster job to %s.." % _platform

    # Create cluster log dir
    cluster_files_dir = \
        os.path.join(_config["output_directory"], "cluster_files")
    if not os.path.exists(cluster_files_dir):
        os.makedirs(cluster_files_dir)

    # Batch file variables
    timestamp = str(strftime("%Y_%m_%d_%H_%M_%S"))
    shell = commands.getoutput('echo $SHELL')
    user_account = getpass.getuser()

    # Set up config dictionary
    config_dict = {'timestamp': timestamp,
                   'shell': shell,
                   'job_name': _run_name,
                   'num_tasks': num_bundles,
                   'queue': "all.q",
                   'par_env': "mpi_smp",
                   'cores_per_task': _num_processors,
                   'user': user_account,
                   'work_dir': cluster_files_dir}

    # Get string template for job scheduler
    if _platform == "PBS":
        env_arr_idx = '$PBS_ARRAYID'
        batch_file_contents = cluster_templates.pbs_template
        confirm_str = '(?<=Your job-array )\d+'
        exec_cmd = 'qsub'
    elif _platform == "SGE":
        env_arr_idx = '$SGE_TASK_ID'
        batch_file_contents = cluster_templates.sge_template
        confirm_str = '(?<=Your job-array )\d+'
        exec_cmd = 'qsub'
    elif _platform == "SLURM":
        hrs_limit = 8*len(subdict)
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
              "%s" % (env_arr_idx, _run_log_dir,
                      _config["subject_list"],
                      _config["pipeline_config_yaml"])

    batch_file_contents = "\n".join([batch_file_contents, run_str])

    batch_filepath = os.path.join(cluster_files_dir, 'cpac_submit_%s.%s'
                                  % (timestamp, _platform))

    with open(batch_filepath, 'w') as f:
        f.write(batch_file_contents)

    print "Batch file written to %s.." % batch_filepath

    # Get output response from job submission
    out = commands.getoutput('%s %s' % (exec_cmd, batch_filepath))

    # Check for successful qsub submission
    if re.search(confirm_str, out) == None:
        err_msg = 'Error submitting QAP pipeline run to %s queue' \
                  % _platform
        raise Exception(err_msg)

    print "Batch job submitted to %s queue." % _platform

    # Get pid and send to pid file
    pid = re.search(confirm_str, out).group(0)
    pid_file = os.path.join(cluster_files_dir, 'pid.txt')
    with open(pid_file, 'w') as f:
        f.write(pid)

def validate_config_dict():
    """Validate the pipeline configuration dictionary to ensure the
    parameters are properly set.
    """
    config_options = ["pipeline_name",
                      "num_processors",
                      "num_participants_at_once",
                      "memory_allocated_per_participant",
                      "resource_manager",
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
    for param in _config.keys():
        if param not in config_options:
            invalid.append(param)
    if len(invalid) > 0:
        err = "\n[!] The following parameters in your configuration " \
              "file are not recognized. Double-check the pipeline " \
              "configuration template.\n"
        err += "\n".join([x for x in invalid])
        raise Exception(err)

def create_flat_sub_dict_dict( subdict):
    """Collapse the participant resource pools so that each participant-
    session-scan combination has its own entry.

    - input subdict format:
          {'sub_01': {'session_01':
                         {'anatomical_scan': {'scan_01': <filepath>,
                                              'scan_02': <filepath>},
                          'site_name': 'Site_1'} },
           'sub_02': {..} }

    - output dict format:
          { (sub01,session01,scan01): {"anatomical_scan": <filepath>,
                                       "anatomical_brain": <filepath>} }

    :type subdict: dictionary
    :param subdict: A dictionary containing the filepaths of input files
                    for each participant, sorted by session and scan.
    :rtype: dictionary
    :return: A dictionary of dictionaries where each participant-session-
             scan combination has its own entry, and input file filepaths
             are defined.
    """

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

def load_sublist():
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

    if "subject_list" in _config.keys():
        with open(_config["subject_list"], "r") as f:
            subdict = yaml.load(f)
    else:
        msg = "\n\n[!] There is no participant list YML to read.\n\n"
        raise_smart_exception(locals(),msg)

    if len(subdict) == 0:
        msg = "The participant list provided is either empty or could " \
              "not be read properly!"
        raise_smart_exception(locals(),msg)

    return subdict

def create_bundles():
    """Create a list of participant "bundles".

    :type flat_sub_dict_dict: dictionary
    :param flat_sub_dict_dict: A dictionary of dictionaries where each
                               participant-session-scan combination has
                               its own entry, and input file filepaths
                               are defined.
    :rtype: list
    :return: A list of bundles - each bundle being a dictionary that is a
             starting resource pool for N sub-session-scan combos with N
             being the number of participants per bundle (set by the user)
    """

    i = 0
    bundles = []

    if len(_sub_dict) < _config["num_participants_at_once"]:
        bundles.append(flat_sub_dict_dict)
    else:
        for sub_info_tuple in _sub_dict.keys():
            if i == 0:
                new_bundle = {}
            new_bundle[sub_info_tuple] = _sub_dict[sub_info_tuple]
            i += 1
            if i == _config["num_participants_at_once"]:
                bundles.append(new_bundle)
                i = 0

        if i > 0:
            bundles.append(new_bundle)

    if len(bundles) == 0:
        msg = "No bundles created."
        raise_smart_exception(locals(),msg)

    return bundles

def run_one_bundle( bundle_idx):
    """Execute one bundle's workflow on one node/slot of a cluster/grid.

    - Compatible with Amazon AWS cluster runs, and S3 buckets.

    :type bundle_idx: int
    :param bundle_idx: The bundle ID number - used to calculate which
                       entries in the participant list to pull into the
                       current bundle, based on the number of participants
                       per bundle (participants at once).
    :rtype: dictionary
    :return: A dictionary with information about the workflow run, its
              status, and results.
    """

    import os
    from qap.qap_workflows_utils import write_json
    from cloud_utils import download_single_s3_path

    _config["workflow_log_dir"] = _run_log_dir

    bundle_dict = _bundles_list[bundle_idx]
    num_bundles = len(_bundles_list)

    # check for s3 paths
    for sub in bundle_dict.keys():
        for resource in bundle_dict[sub].keys():
            value = bundle_dict[sub][resource]
            if "s3://" in value:
                bundle_dict[sub][resource] = download_single_s3_path(value,
                                                                     _config)

    wfargs = (bundle_dict, bundle_dict.keys(),
              _config, _run_name, runargs,
              bundle_idx, num_bundles)

    # let's go!
    rt = build_and_run_qap_pipeline(wfargs)

    # write bundle results to JSON file
    write_json(rt, os.path.join(rt["bundle_log_dir"],
                                "workflow_results.json"))

    # make not uploading results to S3 bucket the default if not specified
    if "upload_to_s3" not in _config.keys():
        _config["upload_to_s3"] = False

    # upload results
    if _config["upload_to_s3"]:
        from cloud_utils import upl_qap_output
        upl_qap_output(_config)

    return rt

def main(config_file=None, partic_list=None):
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


    import sys
    import yaml
    import os
    from time import strftime
    from qap_cfg import default_cfg, write_config, config_output_str
    #from qap.workflow_utils import raise_smart_exception, \
    #                               check_config_settings


    parser = argparse.ArgumentParser()

    in_data_config_group = parser.add_argument_group('Input Data Configuration')

    # Subject list (YAML file)
    in_data_config_group.add_argument('--data_config_file', type=str,
                        help='Path to a data configuration YAML file. This file '
                             'can be automatically generated using the '
                             '--bids_dir, --anat_str, or --func_str options. '
                             'This option will take precedence over any of those'
                             ' options and enables multiple scripts to operate '
                             'on the same dataset without having to regenerate '
                             'the data configuration file for each one. This '
                             'reduces overhead for cluster runs.')

    in_data_config_group.add_argument('--bids_dir', type=str,
                        help='The directory with the input dataset formatted '
                             'according to the BIDS standard. Use the format '
                             's3://bucket/path/to/bidsdir to read data directly '
                             'from an S3 bucket. This may require AWS S3 '
                             'credentials specificied via the --aws_input_creds option.')

    in_data_config_group.add_argument('--anat_str', type=str,
                        help='Template string for generating paths to '
                             'anatomical images. This builds a data '
                             'configuration file for data that is not '
                             'in the BIDS format. Incorporates the '
                             'keywords {site}, {participant}, {session} '
                             'and {run}. Of these only {participant} is '
                             'required. For example: '
                             '/data/{site}/{participant}/{session}/anat/anat_{run}.nii.gz. '
                             'Similar to --bids_dir, cloud data can be specified '
                             'by preprending s3:\\\\<bucketname>\ to the path')

    in_data_config_group.add_argument('--func_str', type=str,
                        help='Template string for generating paths to '
                             'functional images. This builds a data '
                             'configuration file for data that is not '
                             'in the BIDS format. Incorporates the '
                             'keywords {site}, {participant}, {session} '
                             'and {run}. Of these only {participant} is '
                             'required. For example: '
                             '/data/{site}/{participant}/{session}/rest/rest_{run}.nii.gz.'
                             'Similar to --bids_dir, cloud data can be specified '
                             'by preprending s3:\\\\<bucketname>\ to the path'
                        )

    in_data_config_group.add_argument('--aws_input_creds', type=str,
                        help='Credentials for reading from S3. If not provided and s3 paths'
                             'are specified in the data config  we will try to access the '
                             'bucket anonymously',
                        default=None)

    pipeline_config_group = parser.add_argument_group('Pipeline Configuration')
    pipeline_config_group.add_argument('--pipeline_config_file', type=str,
                        help='YAML file specifying QAP configuration (i.e. the arguments '
                             ' to this command). Command line arguments will overload file values.')

    pipeline_config_group.add_argument('--working_dir', type=str,
                        help='The directory to be used for intermediary files. There is no default '
                             'working directory, one must be specified by a pipline configuration '
                             'file or using this arguement.')

    pipeline_config_group.add_argument('--save_working_dir', action='store_true',
                        help='Save the contents of the working directory.', default=False)

    pipeline_config_group.add_argument('--mni_template', type=str,
                        help='MNI template that will be registered to anatomical images '
                             'to exclude neck and lower jaw from image metric '
                             'calculations. default will try to find MNI using $FSLDIR '
                             'environment variable.')

    pipeline_config_group.add_argument('--func_disc_acqs', type=str,
                        help='Number of images that should be discarded from the beginning'
                             'of a fMRI scan. To account for T1 equilibrium effects. default = 0')

    pipeline_config_group.add_argument('--exclude_zeros', action='store_true',
                        help='exclude zero-value voxels from the background of the anatomical scan '
                        'this is meant for images that have been manually altered (ex. faces removed '
                        'for privacy considerations), where the artificial inclusion of zeros into '
                        'the image would skew the QAP metric results, disabled by default')

    output_config_group = parser.add_argument_group('Output Configuration')

    output_config_group.add_argument('--output_dir', type=str,
                        help='The directory where the output files should be stored. There is no default '
                             'working directory, one must be specified by a pipeline configuration '
                             'file or using this arguement.')

    output_config_group.add_argument('--s3_output_creds_path', type=str,
                         help='Credentials for writing to S3.'
                              ' If not provided and s3 paths are specified in the output directory'
                              ' we will try to access the bucket anonymously',
                        default=None)

    output_config_group.add_argument('--report', action='store_true', help='Generates pdf '
                         'for graphically assessing data quality.', default=False)

    output_config_group.add_argument('analysis_level',
                        help='Level of the analysis that will be performed. Multiple'
                             'participant level analyses can be run independently '
                             '(in parallel) using the same output_dir.  Group level '
                             'analysis compiles multiple participant level quality metrics into '
                             'group-level csv files.',
                        choices=['participant', 'group'])

    exec_config_group = parser.add_argument_group('Execution Configuration')
    exec_config_group.add_argument('--n_cpus', type=str,
                        help='Number of execution resources available for the '
                             'pipeline, default=1')

    exec_config_group.add_argument('--mem', help='Amount of RAM available '
                                      ' to the pipeline in GB, default = 4')

    exec_config_group.add_argument('--bundle_size', type=str,
                        help='QAP seperates the work to be performed into '
                             'bundles each of which are executed in parallel '
                             'on a workstation or cluster node. This allows the '
                             'user to balance parallelization on a single node '
                             'with parallizing across nodes to maximize resource '
                             'utilization. Most users will set this to 1.'
                             'default = 1')

    exec_config_group.add_argument('--cluster', choices=['sge','pbs','slurm'],
                        help='Configures QAP to generate and submit a cluster job. '
                             'Valid options are SGE, PBS, or SLURM.')

    exec_config_group.add_argument('--parallel_env', type=str,
                        help='Parallel environment for SGE (see QAP '
                             'web page for configuration recommendation). '
                             'default = mpi_smp.')

    exec_config_group.add_argument('--bundle_idx', type=int,
                            help='Index of bundle to run. This option is '
                                 'primarily for cluster support, but '
                                 'restricts QAP calculation to a subset '
                                 'of the data in the data_config file.')

    args = parser.parse_args()

    #### First lets deal with the pipeline configuration details, if a
    # configuration file is provide read it in, other wise use the
    # default
    if args.pipeline_config_file:
        pipe_config = yaml.load(open(args.pipeline_config_file,'r'))
    else:
        pipe_config = default_cfg

    # get the output directory, make sure its accessible, and if on AWS look for the credentials
    if not args.output_dir and not pipe_config["output_directory"]:
        raise ValueError("Output directory must be specified either in pipeline config file or on the command line.")
    elif args.output_dir:
        pipe_config["output_directory"] = args.output_dir

    if args.s3_output_creds_path:
        pipe_config["s3_output_creds_path"] = args.s3_output_creds_path

    if "s3://" in pipe_config["output_directory"].lower():
        if pipe_config["s3_output_creds_path"]:
            if not os.path.isfile(pipe_config["s3_output_creds_path"]):
                raise ValueError("S3 output credentials ({0}) could not be found.".format(
                    pipe_config["s3_output_creds_path"]))
        else:
            print("Output directory is on AWS, but no credentials provided. Will try write to the bucket anonymously.")

    elif not os.path.exists(pipe_config["output_directory"]):
            raise ValueError("Output directory ({0}) could not be found".format(pipe_config["output_directory"]))

    # get the working directory, make sure its accessible, and that it is not on S3
    if not args.working_dir and not pipe_config["working_directory"]:
        raise ValueError("Working directory must be specified either in pipeline config file or on the command line.")
    elif args.working_dir:
        pipe_config["working_directory"] = args.working_dir

    if "s3://" in pipe_config["working_directory"].lower():
        raise ValueError("QAP does not support writing the working directory to S3 ({0})".format(
                         pipe_config["working_directory"]))

    elif not os.path.exists(pipe_config["working_directory"]):
            raise ValueError("Working directory ({0}) could not be found".format(pipe_config["working_directory"]))

    ### now lets add in the other parameters
    if args.bundle_size:
        pipe_confg["bundle_size"] = int(args.bundle_size)

    if args.n_cpus:
        pipe_config["num_processors"] = int(args.n_cpus)

    if args.mem:
        pipe_config["memory"] = int(args.mem)

    if args.func_disc_acqs:
        pipe_config["start_idx"] = int(args.func_disc_acqs)

    if args.exclude_zeros:
        pipe_config["exclude_zeros"] = True

    if args.mni_template:
        pipe_config["template_head_for_anat"] = args.mni_template

    if not os.path.isfile(pipe_config["template_head_for_anat"]):
        raise ValueError("Could not find MNI template ({0})".format(pipe_config["template_head_for_anat"]))

    if args.report:
        pipe_config["write_report"] = True

    if args.save_working_dir:
        pipe_config["save_working_dir"] = True

    if args.cluster:
        pipe_config['run_on_cluster'] = True
        pipe_config['cluster_system'] = (args.cluster).upper()

    if args.parallel_env:
        pipe_config['parallel_env'] = args.parallel_env

    if 'SGE' in pipe_config['cluster_system'] and not pipe_config['parallel_env']:
        raise ValueError("Running on SGE requires a parallel environment, but none specified")

    print(config_output_str.format(**pipe_config))

    ### write out the pipeline configuration file, nice for debugging and otherstuff
    from datetime import datetime

    pipeline_config_fname = "pipe-config-{}.yml".format(datetime.now().strftime("%Y%m%d%H%M"))

    write_config(pipeline_config_fname, pipe_config)

    sys.exit()

    # the easiest case is when a data configuration file is passed in
    if args.data_config_file:
        data_config = read_yml_file(args.data_config_file)

    elif bids_dir:
        import bids_utils
        data_config = bids_utils.bids_gen_qap_sublist_qap()

    (img_files, config) = collect_bids_files_configs(bids_dir, creds_path)
    if dbg:
        print("Found %d config files for %d image files" % (len(config),
                                                            len(img_files)))

    subdict = bids_gen_cpac_sublist_qap(bids_dir, img_files, creds_path=creds_path, dbg=dbg)


    sys.exit()

    # Load config
    from qap.script_utils import read_yml_file
    _config = read_yml_file(args.pipeline_config_file)
    validate_config_dict()

    _config['pipeline_config_yaml'] = os.path.realpath(args.pipeline_config_file)
    _run_name = _config['pipeline_name']

    if args.with_reports:
        _config['write_report'] = True

    if "num_participants_at_once" not in _config.keys():
        _config["num_participants_at_once"] = 1

    if "num_bundles_at_once" not in _config.keys():
        _config["num_bundles_at_once"] = 1

    _config["subject_list"] = os.path.realpath(args.data_config_file)

    if args.bundle_idx:
        _bundle_idx = args.bundle_idx
    else:
        _bundle_idx = None

    if args.log_dir:
        _run_log_dir = args.log_dir
    else:
        _run_log_dir = None

    # in case we are overloading
    if config_file:
        from qap.script_utils import read_yml_file
        _config = read_yml_file(config_file)
        validate_config_dict()
        _config["pipeline_config_yaml"] = config_file

    if not _config:
         raise Exception("config not found!")

    if partic_list:
        _config["subject_list"] = partic_list

    # Get configurations and settings
    config = _config
    check_config_settings(config, "num_processors")
    check_config_settings(config, "num_participants_at_once")
    check_config_settings(config, "memory_allocated_per_participant")
    check_config_settings(config, "output_directory")
    check_config_settings(config, "working_directory")

    _num_processors = config["num_processors"]
    _num_participants_at_once = config.get('num_participants_at_once', 1)
    _num_bundles_at_once = 1
    write_report = config.get('write_report', False)

    if "resource_manager" in config.keys() and not _bundle_idx:
        res_mngr = config["resource_manager"]
        if (res_mngr == None) or ("None" in res_mngr) or \
            ("none" in res_mngr):
            _platform = None
        else:
            platforms = ["SGE","PBS","SLURM"]
            _platform = str(res_mngr).upper()
            if _platform not in platforms:
                msg = "The resource manager %s provided in the pipeline "\
                      "configuration file is not one of the valid " \
                      "choices. It must be one of the following:\n%s" \
                      % (_platform,str(platforms))
                raise_smart_exception(locals(),msg)
    else:
        _platform = None

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

    results = []

    # set up callback logging
    import logging
    from nipype.pipeline.plugins.callback_log import log_nodes_cb

    cb_log_filename = os.path.join(config["output_directory"],
                                   "callback.log")
    # Add handler to callback log file
    cb_logger = logging.getLogger('callback')
    cb_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(cb_log_filename)
    cb_logger.addHandler(handler)

    # settle run arguments (plugins)
    runargs = {}
    runargs['plugin'] = 'MultiProc'
    memory_for_entire_run = \
        int(config["memory_allocated_per_participant"] * config["num_participants_at_once"])
    runargs['plugin_args'] = {'memory_gb': memory_for_entire_run,
                                   'status_callback': log_nodes_cb}
    n_procs = {'n_procs': _num_processors}
    runargs['plugin_args'].update(n_procs)

    # load the participant list file into dictionary
    subdict = load_sublist()

    try:
        # integrate site information into the subject list
        #   it was separate in the first place to circumvent the fact
        #   that even though site_name doesn't get keyed with scan
        #   names, that doesn't necessarily mean scan names haven't
        #   been specified for that participant

        for sub in subdict.keys():
            sub = str(sub)
            for resource_path in subdict[sub].values():
                if ".nii" in resource_path:
                    filepath = resource_path
                    break

            site_name = filepath.split("/")[-5]
            sub_dict[sub]["site_name"] = site_name

    except:
        pass

    # flatten the participant dictionary
    _sub_dict = create_flat_sub_dict_dict(subdict)

    # create the list of bundles
    _bundles_list = create_bundles()
    num_bundles = len(_bundles_list)

    if not _bundle_idx:
        # want to initialize the run-level log directory (not the bundle-
        # level) only the first time we run the script, due to the
        # timestamp. if sub-nodes are being kicked off by a batch file on
        # a cluster, we don't want a new timestamp for every new node run
        _run_log_dir = op.join(config['output_directory'],
                                    '_'.join([_run_name, "logs"]),
                                    '_'.join([strftime("%Y%m%d_%H_%M_%S"),
                                             "%dbundles" % num_bundles]))

    if _run_log_dir:
        if not os.path.isdir(_run_log_dir):
            try:
                os.makedirs(_run_log_dir)
            except:
                if not op.isdir(_run_log_dir):
                    err = "[!] Log directory unable to be created.\n" \
                          "Path: %s\n\n" % _run_log_dir
                    raise Exception(err)
                else:
                    pass

    if num_bundles == 1:
        _config["num_participants_at_once"] = len(_bundles_list[0])

    # Start the magic
    if not _platform:
        # not a cluster/grid run

        if _num_bundles_at_once == 1:
            # this is always the case
            for idx in range(0, num_bundles):
                results.append(run_one_bundle(idx))

        else:
            # or use Pool if running multiple bundles simultaneously
            # for testing/experimentation only!!!
            from multiprocessing import Pool
            try:
                pool = Pool(processes=_num_bundles_at_once,
                            masktasksperchild=50)
            except TypeError:  # Make python <2.7 compatible
                pool = Pool(processes=_num_bundles_at_once)

            results = pool.map(run_one_bundle, range(0, num_bundles))
            pool.close()
            pool.terminate()

    elif not _bundle_idx:
        # there is a _bundle_idx only if the pipeline runner is run
        # with bundle_idx as a parameter - only happening either manually,
        # or when running on a cluster
        submit_cluster_batch_file(num_bundles)

    else:
        # if there is a bundle_idx supplied to the runner
        results = run_one_bundle(_bundle_idx)

if __name__ == '__main__':
    main()