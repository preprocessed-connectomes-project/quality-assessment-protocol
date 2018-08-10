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


def build_and_run_qap_pipeline(args, run=True):
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
    :type run: bool
    :param run: whether or not the workflow should be run after it is built. If False a pointer to the workflow is
        returned, if True the output from executing the workflow is returned
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
    from qap.qap_utils import read_json

    import glob

    import time
    from time import strftime
    from nipype import config as nyconfig

    # set up Nipype Logger
    from nipype import logging as np_logging
    logger = np_logging.getLogger("workflow")

    # unpack args
    resource_pool_dict, sub_info_list, config, run_name, bundle_idx, num_bundles = args

    # Read and apply general settings in config
    keep_outputs = config["save_working_dir"]

    # take date+time stamp for run identification purposes
    pipeline_start_stamp = strftime("%Y-%m-%d_%H:%M:%S")
    pipeline_start_time = time.time()

    # set up callback logging
    import logging
    from nipype.pipeline.plugins.callback_log import log_nodes_cb

    callback_log_filename = os.path.join(config["log_directory"],
                                         "callback.log")

    # Add handler to callback log file
    callback_logger = logging.getLogger('callback')
    callback_logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(callback_log_filename)
    callback_logger.addHandler(handler)

    # runargs = {'plugin': 'Linear',
    runargs = {'plugin': 'Linear',
               'plugin_args': {'memory_gb': config['available_memory'],
                               'n_procs': config['num_processors'],
                               'status_callback': log_nodes_cb}}

    bundle_log_dir = op.join(config["log_directory"],
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

    # update Nipype logging (not callback)
    nyconfig.update_config({'logging': {'log_directory': bundle_log_dir, 'log_to_file': True}})
    np_logging.update_logging(nyconfig)
    logger.info("QAP version %s" % qap.__version__)
    logger.info("Pipeline start time: %s" % pipeline_start_stamp)

    # create the workflow
    workflow = pe.Workflow(name=run_name)
    workflow.base_dir = config["working_directory"]

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

    logger.info("Starting bundle {0} out of {1} ..".format(str(bundle_idx), str(num_bundles)))

    # initialize results dict
    rt = {'status': 'Started', 'bundle_log_dir': bundle_log_dir}

    # iterate over each subject in the bundle
    for sub_info in sub_info_list:
        resource_pool = resource_pool_dict[sub_info]

        # resource pool check, make sure that we can access each of the files
        invalid_paths = []
        for resource in resource_pool.keys():
            if resource != "site_name" and "s3://" not in resource_pool[resource]:
                if not op.isfile(resource_pool[resource]):
                    invalid_paths.append((resource, resource_pool[resource]))

        if len(invalid_paths) > 0:
            err = "\n\n[!] The following paths in the subject list could not be found:\n"

            for path_tuple in invalid_paths:
                err = "{0}{1}: {2}\n".format(err, path_tuple[0], path_tuple[1])

            err = "{0}\n\n".format(err)
            raise ValueError(err)

        # extract identifiers for this plank of the bundle and conform their names
        def clean_piece(s):
            return str(s).replace("_", "").replace("-", "").replace(".", "")

        site_id = clean_piece(sub_info[0])
        sub_id = clean_piece(sub_info[1])

        if sub_info[2]:
            session_id = clean_piece(sub_info[2])
        else:
            session_id = "session0"

        if sub_info[3]:
            scan_id = clean_piece(sub_info[3])
        else:
            scan_id = "scan0"

        name = "_".join(["", site_id, sub_id, session_id, scan_id])

        rt[name] = {'id': sub_id, 'session': session_id, 'scan': scan_id, 'resource_pool': str(resource_pool),
                    'site': site_id}

        logger.info("Participant info: %s" % name)

        # set output directory
        output_dir = op.join(config["output_directory"], run_name, site_id, sub_id, session_id)

        try:
            os.makedirs(output_dir)
        except:
            # TODO: catch path exists exception and leave all the rest
            if not op.isdir(output_dir):
                err = "[!] Output directory unable to be created.\nPath: %s\n\n".format(output_dir)
                raise Exception(err)
            else:
                pass

        # for QAP spreadsheet generation only
        config.update({"subject_id": sub_id, "session_id": session_id, "scan_id": scan_id, "run_name": run_name})

        if "site_name" in resource_pool:
            config.update({"site_name": resource_pool["site_name"]})
        else:
            config.update({"site_name": site_id})

        logger.info("Configuration settings:\n%s" % str(config))

        qap_types = ["anatomical_spatial", "functional"]

        # update that resource pool with what's already in the output
        # directory
        for resource in os.listdir(output_dir):
            if op.exists(op.join(output_dir, resource)) and resource not in resource_pool.keys():
                try:
                    resource_pool[resource] = glob.glob(op.join(output_dir, resource, "*"))[0]
                except IndexError:
                    if (".json" in resource) and (scan_id in resource):
                        # load relevant json info into resource pool
                        json_file = op.join(output_dir, resource)
                        json_dict = read_json(json_file)
                        try:
                            sub_json_dict = json_dict["%s %s %s".format(sub_id, session_id, scan_id)]
                        except KeyError:
                            err = "Reading wrong JSON file?\n{0}\nResource: {1} \nScan: {2}".format(json_file, resource,
                                                                                                    scan_id)
                            raise Exception(err)

                        if "anatomical_header_info" in sub_json_dict.keys():
                            resource_pool["anatomical_header_info"] = sub_json_dict["anatomical_header_info"]

                        if "functional_header_info" in sub_json_dict.keys():
                            resource_pool["functional_header_info"] = sub_json_dict["functional_header_info"]

                        for qap_type in qap_types:
                            if qap_type in sub_json_dict.keys():
                                resource_pool["_".join(["qap", qap_type])] = sub_json_dict[qap_type]
                except:
                    # a stray file in the sub-sess-scan output directory
                    pass

        # add the starter node to the resource pool, which connects workflows for each plank of the bundle
        # into a larger pipeline
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
                wf_builder = getattr(qw, "_".join(["qap", qap_type, "workflow"]))
                workflow, resource_pool = wf_builder(workflow, resource_pool, config, name)

        if ("anat_scan" in resource_pool.keys()) and ("anat_header_info" not in resource_pool.keys()):
            if qw is None:
                from qap import qap_workflows as qw
            workflow, resource_pool = qw.qap_gather_header_info(workflow, resource_pool, config, name, "anat")

        if ("func_scan" in resource_pool.keys()) and ("func_header_info" not in resource_pool.keys()):
            if qw is None:
                from qap import qap_workflows as qw
            workflow, resource_pool = qw.qap_gather_header_info(workflow, resource_pool, config, name, "func")

        # set up the datasinks
        out_list = []
        for output in resource_pool.keys():
            for qap_type in qap_types:
                if qap_type in output:
                    out_list.append("_".join(["qap", qap_type]))

        filepath_keep = ["anat_reorient", "anat_fav_artifacts_background", "func_mean", "func_estimated_nuisance",
                         "func_temporal_std_map", "func_SFS"]

        if keep_outputs:
            # write_all_outputs (writes everything to the output
            # directory, not just the final JSON files)
            out_list = resource_pool.keys()
        else:
            for resource in filepath_keep:
                if resource in resource_pool.keys():
                    out_list.append(resource)

        logger.info("Outputs we're keeping: %s" % str(out_list))
        logger.info('Resource pool keys after workflow connection: '
                    '{}'.format(str(resource_pool.keys())))

        # Save reports to out_dir if necessary
        if config.get('write_report', False):

            if ("qap_mosaic" in resource_pool.keys()) and \
                    ("qap_mosaic" not in out_list):
                out_list += ['qap_mosaic']

            # The functional temporal also has an FD plot
            if 'qap_functional' in resource_pool.keys():
                if ("qap_fd" in resource_pool.keys()) and \
                        ("qap_fd" not in out_list):
                    out_list += ['qap_fd']

        for output in out_list:
            # we use a check for len()==2 here to select those items in
            # the resource pool which are tuples of (node, node_output),
            # instead of the items which are straight paths to files

            # resource pool items which are in the tuple format are the
            # outputs that have been created in this workflow because they
            # were not present in the subject list YML (the starting
            # resource pool) and had to be generated
            if (len(resource_pool[output]) == 2) and (output != "starter"):

                node, out_file = resource_pool[output]

                # create the datasink
                ds = pe.Node(nio.DataSink(), name='datasink_{0}{1}'.format(output, name))
                ds.inputs.base_directory = output_dir

                # rename file to BIDS format
                rename = pe.Node(niu.Rename(), name='rename_{0}{1}'.format(output, name))
                rename.inputs.keep_ext = True
                # replace the underscores in 'output' (which are the keys
                # of the resource pool) with dashes - need to be
                # underscores for the workflow names, but need to be
                # dashes for BIDS output file naming format
                rename.inputs.format_string = "{0}_{1}_{2}_{3}".format(sub_id, session_id, scan_id,
                                                                       output.replace("_", "-"))
                workflow.connect(node, out_file, rename, 'in_file')
                workflow.connect(rename, 'out_file', ds, '{0}.@{1}'.format(output.split("_")[0], output))
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
        if run is True:
            try:
                logger.info("Running with plugin %s" % runargs["plugin"])
                logger.info("Using plugin args %s" % runargs["plugin_args"])
                workflow.run(plugin=runargs["plugin"],
                             plugin_args=runargs["plugin_args"])
                rt['status'] = 'finished'
                logger.info("Workflow run finished for bundle %s." % str(bundle_idx))
            except Exception as e:  # TODO We should be more specific here ...
                rt.update({'status': 'failed', 'msg': e})
                logger.info("Workflow run failed for bundle %s." % str(bundle_idx))
                # ... however this is run inside a pool.map: do not raise
                # Exception
        else:
            return workflow

    else:
        rt['status'] = 'cached'
        logger.info("\nEverything is already done for bundle %s." % str(bundle_idx))

    # Remove working directory when done
    if not config['save_working_dir']:
        try:
            if op.exists(config["working_directory"]):
                import shutil
                shutil.rmtree(config["working_directory"])
        except:
            logger.warn("Couldn\'t remove the working directory!")
            pass

    if rt["status"] == "failed":
        logger.error(rt["msg"])
    else:
        pipeline_end_stamp = strftime("%Y-%m-%d_%H:%M:%S")
        pipeline_end_time = time.time()
        logger.info("Elapsed time (minutes) since last start: %s"
                    % ((pipeline_end_time - pipeline_start_time) / 60))
        logger.info("Pipeline end time: %s" % pipeline_end_stamp)

    return rt
