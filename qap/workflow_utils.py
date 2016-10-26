

def raise_smart_exception(local_vars, msg=None):
    """Raise an exception with more information about the traceback, and 
    enforce inclusion of the locals().

    Keyword Arguments:
      local_vars -- [Python dictionary] input for locals()
      msg -- [string] (default: None) the custom error message for the 
             exception in question

    Returns:
      N/A
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

    Keyword Arguments:
      resource_pool -- [Python dictionary] the resource pool of resources 
                       (which includes output files of sub-workflows, and 
                       connection pointers for Nipype nodes/workflows)
      resource_name -- [string] the name of the output/intermediary file to 
                       check for within the resource pool

    Returns:
      N/A
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

    Keyword Arguments:
      config -- [Python dictionary] a dictionary keying configuration options 
                to their chosen selections
      parameter -- [string] the key of the configuration parameter to be 
                   checked

    Returns:
      N/A
    """

    if parameter not in config.keys():
        err = "[!] The parameter '%s' is missing from your pipeline " \
              "configuration .YML file. Please make sure this is specified " \
              "properly." % parameter
        raise_smart_exception(locals(),err)
