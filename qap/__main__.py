#!/usr/bin/env python


def main():
    """

    Function to parse the command line arguments and locally execute QAP using the 'multiproc' plugin. This command
    line interface can be executed on cluster systems (SGE, SLURM, etc), but this CLI does not explicitly handle
    creating the job submission files or submitting them. This file is made consistent with the BIDS-App specification
    so that it is the entry point of the QAP BIDS-App container (woo-hoo!).

    """

    import os
    import sys
    import yaml
    import argparse
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


if __name__ == "__main__":
    main()