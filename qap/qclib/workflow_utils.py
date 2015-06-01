
def check_input_resources(resource_pool, resource_name):

    import os

    if resource_name not in resource_pool.keys():
    
        err = "\n\n[!] The resource '%s' is missing from the resource " \
              "pool, and it is needed in one of the steps of the pipeline. " \
              "Please make sure it is specified properly.\n\n" % \
              resource_name
              
        raise Exception(err)

    else:

        if len(resource_pool[resource_name]) > 2:

            if not os.path.isfile(resource_pool[resource_name]):
        
                err = "\n\n[!] The path provided for the resource '%s' " \
                      "does not exist!\nPath provided: %s\n\n" % \
                      (resource_name, resource_pool[resource_name])
              
                raise Exception(err)
            
            
            
def check_config_settings(config, parameter):

    if parameter not in config.keys():
    
        err = "\n\n[!] The parameter '%s' is missing from your pipeline " \
              "configuration .YML file. Please make sure this is specified " \
              "properly.\n\n" % parameter
              
        raise Exception(err)
