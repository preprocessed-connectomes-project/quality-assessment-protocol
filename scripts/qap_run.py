#!/usr/bin/env python
# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:

import argparse
from qap import bids_utils
from qap import qap_cfg

from qap.qap_pipeline import build_and_run_qap_pipeline


def create_bundles(data_configuration_dict, bundle_size):
    """

    Create a list of participant "bundles".

    :param data_configuration_dict: A nested dictionary session->participant->session->scan->{resource pool}
    :param bundle_size: the number of sub-session-session combos (tranches) in a bundle

    :return: A list of bundles - each bundle being a dictionary that is a starting resource pool for N
    sub-session-session combos with bundle_size being set by the user
    """

    tranche_count = 0
    bundles = []
    new_bundle = {}

    for site in data_configuration_dict.keys():
        for participant in data_configuration_dict[site].keys():
            for session in data_configuration_dict[site][participant].keys():

                # each tranche corresponds to a scanning session, a plank of
                # the tranche corresponds to a scan collected during the
                # session
                if tranche_count % bundle_size == 0:
                    if new_bundle:
                        bundles.append(new_bundle)
                    new_bundle = {}

                for scan in data_configuration_dict[site][participant][session].keys():
                    tranche_plank_identifier = (site, participant, session, scan)

                    if tranche_plank_identifier in new_bundle:
                        print("WARNING: tranche {0}, {1}, {2}, {3} already "
                              "found in bundle, overwriting".format(tranche_plank_identifier[0],
                                                                    tranche_plank_identifier[1],
                                                                    tranche_plank_identifier[2],
                                                                    tranche_plank_identifier[3]))

                    new_bundle[tranche_plank_identifier] = data_configuration_dict[site][participant][session][scan]

                tranche_count += 1

    # if there is a bundle still hanging around, add it to the list
    if new_bundle:
        bundles.append(new_bundle)

    return bundles


def write_bundles(bundles, bundle_filename):

    import yaml
    if not bundle_filename.endswith('yml'):
        bundle_filename += '.yml'

    with open(bundle_filename, 'w') as ofd:
        yaml.dump(bundles, ofd, encoding='utf-8')


def run_qap(pipeline_configuration, data_bundles, bundle_index=None):
    """

    Build and execute the qap pipeline based on the configuration parameters.

    :param pipeline_configuration: dictionary containing execution parameters
    :param data_bundles: list of data dictionaries that specify data to be processed
    :param bundle_index: index of a particular bundle to run

    :return: A list of dictionaries with information about the workflow runs, their
              status, and results.
    """

    import os
    from qap.qap_utils import write_json
    from qap.cloud_utils import download_single_s3_path, copy_directory_to_s3

    # first deal with logging, if the plan is to write the log files to S3, we
    # first write them locally to the working directory and then move them to
    # S3
    s3_log_path = ''
    if "s3://" in pipeline_configuration["log_directory"]:
        s3_log_path = pipeline_configuration["log_directory"]
        pipeline_configuration["log_directory"] = os.path.join(pipeline_configuration["working_directory"], "logs")

    if not os.path.isdir(pipeline_configuration["log_directory"]):
        os.makedirs(pipeline_configuration["log_directory"])

    # set up callback logging
    import logging
    # from nipype.pipeline.plugins.callback_log import log_nodes_cb

    callback_log_filename = os.path.join(pipeline_configuration["log_directory"], "callback.log")

    # Add handler to callback log file
    callback_logger = logging.getLogger('callback')
    callback_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(callback_log_filename)
    callback_logger.addHandler(handler)

    num_bundles = len(data_bundles)

    if bundle_index:
        if bundle_index < 0 or bundle_index >= num_bundles:
            raise ValueError(
                "Invalid bundle_index {0}, should be in the range 0 < x < {1}".format(bundle_index, num_bundles))

        data_bundles = [data_bundles[int(bundle_index)]]

    results = []
    for index, bundle_dict in enumerate(data_bundles):

        #  the QAP pipeline uses the bundle index in the names of various files
        #  and directories, we want to preserve the actual bundle index to make
        #  sure that values in these names correctly map back out to the actual
        #  data, so, if we were passed a bundle index, use it, otherwise use
        #  our internally calculated index
        if not bundle_index:
            bundle_index = index

        # if input values are in s3, go get 'em
        for tranche_key, tranche_dict in bundle_dict.iteritems():
            for resource_key, resource_path in tranche_dict.iteritems():
                if "s3://" in resource_path.lower():
                    bundle_dict[tranche_key][resource_key] = \
                        download_single_s3_path(resource_path,
                                                pipeline_configuration)

        qap_pipeline_arguments = (bundle_dict, bundle_dict.keys(), pipeline_configuration, "qap_run", bundle_index,
                                  num_bundles)

        pipeline_execution_output = build_and_run_qap_pipeline(qap_pipeline_arguments)

        # save the output to a json file in the log directory
        write_json(pipeline_execution_output, os.path.join(pipeline_execution_output["bundle_log_dir"],
                                                           "workflow_results.json"))

        if s3_log_path:

            s3_path = pipeline_execution_output["bundle_log_dir"].replace(pipeline_configuration["log_directory"],
                                                                          s3_log_path)

            print("Copying log files from {0} to {1}".format(pipeline_execution_output["bundle_log_dir"], s3_path))
            copy_directory_to_s3(pipeline_execution_output["bundle_log_dir"], s3_path, pipeline_configuration)

        results.append(pipeline_execution_output)

    return results


def main():
    """

    Function to parse the command line arguments and locally execute QAP using the 'multiproc' plugin. This command
    line interface can be executed on cluster systems (SGE, SLURM, etc), but this CLI does not explicitly handle
    creating the job submission files or submitting them. This file is made consistent with the BIDS-App specification
    so that it is the entry point of the QAP BIDS-App container (woo-hoo!).

    """

    import sys
    import yaml
    import os
    import pkg_resources as p
    from datetime import datetime

    parser = argparse.ArgumentParser()

    in_data_config_group = parser.add_argument_group('Input Data Configuration')

    # Subject list (YAML file)
    in_data_config_group.add_argument('--data_config_file', type=str,
                                      help='Path to a data configuration YAML file. This file can be automatically'
                                           ' generated using the --bids_dir, --anat_str, or --func_str options. This'
                                           ' option will take precedence over any of those options and enables multiple'
                                           ' scripts to operate on the same dataset without having to regenerate the'
                                           ' data configuration file for each one. This reduces overhead for cluster'
                                           ' runs.',
                                      default=None)

    in_data_config_group.add_argument('bids_dir', type=str,
                                      help='The directory with the input dataset formatted according to the BIDS'
                                           ' standard. Use the format s3://bucket_name/path/to/bidsdir to read data'
                                           ' directly from an S3 bucket. This may require AWS S3 credentials'
                                           ' to be specified using the --aws_input_credentials option.',
                                      default=None)

    in_data_config_group.add_argument('--anatomical_str', type=str,
                                      help='Template string for generating paths to anatomical images. This builds '
                                           'a data configuration file for data that is not in the BIDS format. '
                                           'Incorporates the keywords {site}, {participant}, {session} and {run}. Of '
                                           'these only {participant} is required. For example: '
                                           '/data/{site}/{participant}/{session}/anat/anat_{run}.nii.gz. Similar to '
                                           '--bids_dir, to specify cloud data prepend s3://<bucket_name>/ to '
                                           'the path',
                                      default=None)

    in_data_config_group.add_argument('--functional_str', type=str,
                                      help='Template string for generating paths to functional images. This builds '
                                           'a data configuration file for data that is not \in the BIDS format. '
                                           'Incorporates the keywords {site}, {participant}, {session} and {run}. Of '
                                           'these only {participant} is required. For example: '
                                           '/data/{site}/{participant}/{session}/rest/rest_{run}.nii.gz.'
                                           'Similar to --bids_dir, to specify cloud data prepend s3://<bucket_name>/ '
                                           'to the path',
                                      default=None)

    in_data_config_group.add_argument('--s3_read_credentials', type=str,
                                      help='Path to file containing credentials for reading from S3. If not provided '
                                           'and s3 paths are specified in the data config we will try to access the '
                                           'bucket anonymously',
                                      default=None)

    pipeline_config_group = parser.add_argument_group('Pipeline Configuration')
    pipeline_config_group.add_argument('--pipeline_config_file', type=str,
                                       help='Path to YAML file specifying QAP configuration (i.e. the arguments '
                                            'to this command). Command line arguments will overload file values. To '
                                            'specify a file in s3, prepend s3://<bucket_name>/ to the path.',
                                       default=p.resource_filename("qap", os.path.join("configs", "qap_pipe_config_template.yml")))

    pipeline_config_group.add_argument('--working_dir', type=str,
                                       help='The directory to be used for intermediary files. Can be specified by a '
                                            'pipeline configuration file or using this argument. Default = /tmp',
                                       default=None)

    pipeline_config_group.add_argument('--recompute_all_derivatives', action='store_true',
                                       help='Recompute all outputs, even if they already exist in output directory.',
                                       default=False)

    pipeline_config_group.add_argument('--save_working_dir', action='store_true',
                                       help='Save the contents of the working directory.',
                                       default=False)

    pipeline_config_group.add_argument('--mni_template', type=str,
                                       help='MNI template that will be registered to anatomical images to exclude neck '
                                            'and lower jaw from image metric calculations. Default will search the '
                                            '$PATH for the MNI_avg152T1+tlrc template distributed with AFNI.',
                                       default=None)

    pipeline_config_group.add_argument('--functional_discard_volumes', type=str,
                                       help='Number of volumes that should be discarded from the beginning'
                                            'of a fMRI scan to account for T1 equilibrium effects. default = 0',
                                       default=None)

    pipeline_config_group.add_argument('--exclude_zeros', action='store_true',
                                       help='Exclude zero-value voxels from the background of the anatomical scan. '
                                            'This is meant for images that have been manually altered (ex. faces '
                                            'removed for privacy considerations), where the artificial inclusion of '
                                            'zeros into the image would skew the QAP metric results. '
                                            'Disabled by default.',
                                       default=False)

    output_config_group = parser.add_argument_group('Output Configuration')

    output_config_group.add_argument('output_dir', type=str,
                                     help='The directory where the output files should be stored. There is no default '
                                          'output directory, one must be specified by a pipeline configuration '
                                          'file or using this argument.',
                                     default=None)

    output_config_group.add_argument('--s3_write_credentials', type=str,
                                     help='Path to file containing credentials for writing to S3. If not provided and '
                                          's3 paths are specified in the output directory we will try to access the '
                                          'bucket anonymously',
                                     default=None)

    output_config_group.add_argument('--log_dir', type=str,
                                     help='The directory where log files should be stored. If not specified, log files '
                                          'will be written to [output_dir]/logs.',
                                     default=None)

    output_config_group.add_argument('--report', action='store_true',
                                     help='Generates pdf for graphically assessing data quality.',
                                     default=False)

    output_config_group.add_argument('analysis_level',
                                     help='Level of the analysis that will be performed. Multiple participant level '
                                          'analyses can be run independently (in parallel) using the same output_dir. '
                                          'Group level analysis compiles multiple participant level quality metrics '
                                          'into group-level csv files. test_config will check all of the specified '
                                          'parameters for consistency and write out participant and qap configuration '
                                          'files but will not execute the pipeline.',
                                     choices=['participant', 'group', 'test_config'])

    exec_config_group = parser.add_argument_group('Execution Configuration')
    exec_config_group.add_argument('--n_cpus', type=str,
                                   help='Number of execution resources available for the pipeline. If not specified in '
                                        'the pipeline configuration file, the default will be 1.',
                                   default=None)

    exec_config_group.add_argument('--mem_mb', type=str,
                                   help='Amount of RAM available to the pipeline in MB. If not specified in the '
                                        'pipeline configuration file the default will be 4096. This is '
                                        'included to be consistent with the BIDS-APP specification, but mem_gb '
                                        'is preferred. This flag is ignored if mem_gb is also specified.',
                                   default=None)

    exec_config_group.add_argument('--mem_gb', type=str,
                                   help='Amount of RAM available to the pipeline in GB. If not specified in the '
                                        'pipeline configuration file the default will be 4. This argument '
                                        'takes precedence over --mem_mb.',
                                   default=None)

    exec_config_group.add_argument('--bundle_size', type=str,
                                   help='QAP separates the work to be performed into bundles each of which are '
                                        'executed in parallel on a workstation or cluster node. This allows the '
                                        'user to balance parallelization on a single node with parallelizing across '
                                        'nodes to maximize resource utilization. Most users will set this to 1.'
                                        'If this is not specified in the pipeline configuration file the default will '
                                        'be 1',
                                   default=None)

    exec_config_group.add_argument('--bundle_index', type=int,
                                   help='Index of bundle to run. This option is primarily for cluster support and '
                                        'restricts QAP calculation to a subset of the data in the data_config file.',
                                   default=None)

    args = parser.parse_args()

    # First lets deal with the pipeline configuration details, if a
    # configuration file is provide read it in, other wise use the
    # default
    pipeline_configuration = qap_cfg.default_pipeline_configuration

    if args.pipeline_config_file:

        if "s3://" in args.pipeline_config_file.lower():
            from indi_aws import fetch_creds

            s3_bucket_name = args.pipeline_config_file.split('/')[2]
            s3_prefix = '/'.join(args.pipeline_config_file.split('/')[:3])
            pipeline_config_key = args.pipeline_config_file.replace(s3_prefix, '').lstrip('/')

            bucket = fetch_creds.return_bucket(args.s3_read_credentials, s3_bucket_name)
            pipeline_configuration.update(yaml.load(bucket.Object(pipeline_config_key).get()["Body"].read()))
        else:
            pipeline_configuration.update(yaml.load(open(args.pipeline_config_file, 'r')))

    # get the output directory, make sure its accessible, and if on AWS look
    # for the credentials
    if not args.output_dir and not pipeline_configuration["output_directory"]:
        raise ValueError("Output directory must be specified either in "
                         "pipeline config file or on the command line.")
    elif args.output_dir:
        output_dir = os.path.abspath(args.output_dir)
        pipeline_configuration["output_directory"] = output_dir

    if "s3://" in pipeline_configuration["output_directory"].lower():

        if pipeline_configuration["s3_write_credentials"]:
            if not os.path.isfile(pipeline_configuration["s3_write_credentials"]):
                raise ValueError("S3 output credentials ({0}) could not be "
                                 "found.".format(pipeline_configuration["s3_write_credentials"]))
        else:
            print('Output directory is on S3, but no write credentials were '
                  'provided. Will try write to the bucket anonymously.')

    elif not os.path.exists(pipeline_configuration["output_directory"]):
            raise ValueError('Output directory ({0}) could not be '
                             'found'.format(pipeline_configuration["output_directory"]))

    # get the logging directory, this can go to S3, but will have to be
    # handled differently than outputs which are written to S3 by the datasink
    if not args.log_dir and not pipeline_configuration["log_directory"]:
        pipeline_configuration["log_directory"] = pipeline_configuration["output_directory"]+"/logs"
    elif args.log_dir:
        pipeline_configuration["log_directory"] = args.log_dir

    # get the working directory, make sure its accessible, and that it is not
    # on S3
    if not args.working_dir and not pipeline_configuration["working_directory"]:
        raise ValueError("Working directory must be specified either in "
                         "pipeline config file or on the command line.")
    elif args.working_dir:
        pipeline_configuration["working_directory"] = args.working_dir

    if "s3://" in pipeline_configuration["working_directory"].lower():
        raise ValueError("QAP does not support writing the working directory to S3 ({0})".format(
            pipeline_configuration["working_directory"]))

    elif not os.path.exists(pipeline_configuration["working_directory"]):
           os.makedirs(pipeline_configuration["working_directory"])

    # now lets add in the other parameters
    if args.bundle_size:
        pipeline_configuration["bundle_size"] = int(args.bundle_size)

    if args.n_cpus:
        pipeline_configuration["num_processors"] = int(args.n_cpus)

    if args.mem_gb:
        pipeline_configuration["available_memory"] = float(args.mem_gb)
    elif args.mem_mb:
        pipeline_configuration["available_memory"] = float(args.mem_mb)/1024.0

    if args.functional_discard_volumes:
        pipeline_configuration["functional_start_index"] = int(args.functional_discard_volumes)
        if pipeline_configuration["functional_start_index"] < 0:
            raise ValueError("functional_start_index cannot be a negative number {}".format(
                pipeline_configuration["functional_start_index"]))

    if args.exclude_zeros:
        pipeline_configuration["exclude_zeros"] = True

    if args.mni_template:
        pipeline_configuration["anatomical_template"] = args.mni_template
    else:
        # if the user does not provide a anatomical template, search the path
        # environment variable to see if we can find the default AFNI template
        # MNI_avg152T1
        pipeline_configuration['anatomical_template'] = ''

        for path in os.environ['PATH'].split(':'):
            if os.path.isfile(os.path.join(path, 'MNI_avg152T1+tlrc.HEAD')):
                pipeline_configuration['anatomical_template'] = os.path.join(path, 'MNI_avg152T1+tlrc.BRIK.gz')
                break
                
        if not pipeline_configuration['anatomical_template']:
            raise ValueError('No MNI template specified and could not '
                             'find MNI_avg152T1+tlrc on the system PATH.')

    if not os.path.isfile(pipeline_configuration["anatomical_template"]) and \
            not os.path.isfile(pipeline_configuration["anatomical_template"]+".HEAD"):
        raise ValueError("Could not find MNI template ({0})".format(pipeline_configuration["anatomical_template"]))

    if args.report:
        pipeline_configuration["write_report"] = True

    if args.save_working_dir:
        pipeline_configuration["save_working_dir"] = True

    if args.recompute_all_derivatives:
        pipeline_configuration["recompute_all_derivatives"] = True

    qap_cfg.validate_pipeline_configuration(pipeline_configuration)

    print(qap_cfg.configuration_output_string.format(**pipeline_configuration))

    # write out the pipeline configuration file, nice for debugging and other
    # stuff

    timestamp_string = datetime.now().strftime("%Y%m%d%H%M")

    pipeline_configuration_filename = os.path.join(args.output_dir, "qap-pipe-config-{}.yml".format(timestamp_string))
    qap_cfg.write_pipeline_configuration(pipeline_configuration_filename, pipeline_configuration)

    # the easiest case is when a data configuration file is passed in
    data_configuration_dict = {}
    if args.data_config_file:
        data_configuration_dict = yaml.load(open(args.data_config_file, 'r'))

    elif args.bids_dir:

        bids_dir = os.path.abspath(args.bids_dir)

        (bids_image_files, bids_derivative_files, bids_parameters_dictionary) = \
            bids_utils.collect_bids_files_configs(bids_dir,
                                                  args.s3_read_credentials)

        if not pipeline_configuration['recompute_all_derivatives']:

            # since bids derivatives also live in the bids_dir, it is entirely
            # possible that the user set the out_dir to the same value as
            # bids_dir
            out_derivative_files = []
            if bids_dir != pipeline_configuration["output_directory"]:
                (out_image_files, out_derivative_files, out_parameters_dictionary) = \
                    bids_utils.collect_bids_files_configs(pipeline_configuration["output_directory"],
                                                          args.s3_write_credentials,
                                                          raise_error=False)

            print("Found {0} parameter dictionaries, {1} image files, and {2} derivatives".format(
                len(bids_parameters_dictionary), len(bids_image_files),
                len(bids_derivative_files) + len(out_derivative_files)))

            # add the derivatives to the files so that they will all be added
            # to the data configuration
            bids_image_files = bids_image_files + bids_derivative_files + out_derivative_files
        else:
            print("Found {0} parameter dictionaries, and {1} image files".format(
                len(bids_parameters_dictionary), len(bids_image_files)))

        data_configuration_dict = \
            bids_utils.bids_generate_qap_data_configuration(bids_dir, bids_image_files,
                                                            credentials_path=args.s3_read_credentials, debug=True)

    data_configuration_filename = os.path.join(args.output_dir, "qap-data-config-{}.yml".format(timestamp_string))
    bids_utils.write_data_configuration(data_configuration_filename, data_configuration_dict)

    if "test_config" in args.analysis_level:
        print("Configuration has been tested and appears to be OK. Pipeline "
              "configuration and data configuration can be found in {0} and "
              "{1} respectively.".format(pipeline_configuration_filename,
                                         data_configuration_filename))
        sys.exit(0)

    elif "participant" in args.analysis_level:
        print("Configuration has been tested and appears to be OK. Executing "
              "pipeline with pipeline configuration {0} and data "
              "configuration {1}".format(pipeline_configuration_filename, data_configuration_filename))

        # prepare the bundles
        data_bundles = create_bundles(data_configuration_dict,
                                      pipeline_configuration["bundle_size"])
        # data_bundle_filename = "qap-data-bundles-{}.yml".format(timestamp_string)
        # write_bundles(data_bundles, data_bundle_filename)

        num_bundles = len(data_bundles)
        print("Created {} bundle(s)".format(num_bundles))

        # --- Now Execute
        bundle_index = ''
        if args.bundle_index:
            bundle_index = int(args.bundle_index)

        run_qap(pipeline_configuration, data_bundles, bundle_index)


if __name__ == '__main__':
    main()
