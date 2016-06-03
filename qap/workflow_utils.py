

def raise_smart_exception(local_vars, msg=None):

    import traceback
    e = "\n\nLocal variables:\n%s\n\n%s\n\n" \
        % (str(local_vars), str(traceback.format_exc()))
    if msg:
        e = e + "\n\n%s\n\n" % str(msg)
    raise Exception(e)



def check_input_resources(resource_pool, resource_name):

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

    if parameter not in config.keys():
    
        err = "[!] The parameter '%s' is missing from your pipeline " \
              "configuration .YML file. Please make sure this is specified " \
              "properly." % parameter
              
        raise_smart_exception(locals(),err)
