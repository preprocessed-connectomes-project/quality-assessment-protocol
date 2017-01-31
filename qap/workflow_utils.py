

def raise_smart_exception(local_vars, msg=None):
    """Raise an exception with more information about the traceback, and 
    enforce inclusion of the locals().

    :type local_vars: dict
    :param local_vars: Input for locals().
    :type msg: str
    :param msg: (default: None) The custom error message for the exception in
                question.
    """

    import traceback
    e = "\n\nLocal variables:\n%s\n\n%s\n\n" \
        % (str(local_vars), str(traceback.format_exc()))
    if msg:
        e = "%s\n\n%s\n\n" % (e, str(msg))
    raise Exception(e)


def check_input_resources(resource_pool, resource_name):
    """Check to make sure a specific resource/file is present in the 
    resource pool.

    :type resource_pool: dict
    :param resource_pool: The resource pool of resources (which includes
                          output files of sub-workflows, and connection
                          pointers for Nipype nodes/workflows).
    :type resource_name: str
    :param resource_name: The name of the output/intermediary file to check
                          for within the resource pool.
    """

    import os

    if resource_name not in resource_pool.keys():
        err = "Resource pool: %s\n\n[!] The resource '%s' is missing " \
              "from the resource pool, and it is needed in one of the steps "\
              "of the pipeline. Please make sure it is specified " \
              "properly." % (resource_pool, resource_name)
              
        raise_smart_exception(locals(),err)

    else:
        if len(resource_pool[resource_name]) > 2:
            if not os.path.isfile(resource_pool[resource_name]):
                err = "[!] The path provided for the resource '%s' " \
                      "does not exist!\nPath provided: %s" % \
                      (resource_name, resource_pool[resource_name])
              
                raise_smart_exception(locals(),err)
            
            
def check_config_settings(config, parameter):
    """Check to make sure a configuration setting/parameter is present in the 
    pipeline configuration dictionary.

    :type config: dict
    :param config: A dictionary keying configuration options to their chosen
                   selections.
    :type parameter: str
    :param parameter: The key of the configuration parameter to be checked.
    """

    if parameter not in config.keys():
        err = "[!] The parameter '%s' is missing from your pipeline " \
              "configuration .YML file. Please make sure this is specified " \
              "properly." % parameter
        raise_smart_exception(locals(),err)


def generate_nipype_workflow_graphs(workflow, out_dir=None):
    """Generate the Nipype workflow dependency graphs given the workflow 
    object.

    :type workflow: Nipype workflow object
    :param workflow: The connected workflow object.
    :type out_dir: str
    :param out_dir: (default: None) The directory where to write the
                    dependency graph .dot and .png files to.
    """

    if not out_dir:
        pass

    """
    workflow.write_graph(
        dotfilename=op.join(config["output_directory"], \
                            "".join([run_name, ".dot"])),
        simple_form=False)
    workflow.write_graph(
        graph2use="orig",
        dotfilename=op.join(config["output_directory"], \
                            "".join([run_name, ".dot"])),
        simple_form=False)
    workflow.write_graph(
        graph2use="hierarchical",
        dotfilename=op.join(config["output_directory"], \
                            "".join([run_name, ".dot"])),
        simple_form=False)
    """