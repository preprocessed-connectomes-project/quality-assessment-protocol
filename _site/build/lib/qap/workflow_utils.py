
def check_input_resources(resource_pool, resource_name):

    import os

    if resource_name not in resource_pool.keys():
    
        err = "\n\nResource pool: %s\n\n[!] The resource '%s' is missing " \
              "from the resource pool, and it is needed in one of the steps "\
              "of the pipeline. Please make sure it is specified " \
              "properly.\n\n" % (resource_pool, resource_name)
              
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



def build_test_case(workflow, ref_inputs, ref_graph, wf_inputs_string):

    import os

    # get the reference workflow input string
    with open(ref_inputs,"r") as f:
        ref_inputs_string = f.read()


    ref_inputs_string = ref_inputs_string.replace("\n","")


    # write the dependency graph of the workflow we are testing
    out_graph = os.path.join(os.getcwd(), "graph.dot")
    workflow.write_graph(dotfilename=out_graph, simple_form=False)
    
    
    # load the both the reference and the to-test dependency graphs
    with open(ref_graph,"r") as f:
        ref_graph_lines = sorted(f.readlines())

    with open(out_graph,"r") as f:
        out_graph_lines = sorted(f.readlines())
                   
    
    # clear generated graph.dot
    try:
        os.system("rm %s" % os.path.join(os.getcwd(), "graph.dot"))
    except:
        pass


    # test the case
    flag = 0
    
    err = ""

    if ref_graph_lines == out_graph_lines:
        flag += 1
    else:
        err = err + "\n" + str(ref_graph_lines)
        err = err + "\n\n" + str(out_graph_lines)
        err = err + "\n\nDependency graph TEST FAIL\n\n"
        
    if ref_inputs_string == wf_inputs_string:
        flag += 1
    else:
        err = err + "\n" + ref_inputs_string
        err = err + "\n\n" + wf_inputs_string
        err = err + "\n\nWorkflow inputs TEST FAIL\n\n"


    return flag, err



def generate_test_files(workflow, test_data_folder):

    import os

    workflow_name = test_data_folder.split("/")[-1]

    graph_path = os.path.join(test_data_folder, "graph_%s.dot" % \
                                  workflow_name)

    workflow.write_graph(dotfilename=graph_path, simple_form=False)

    inputs_path = os.path.join(test_data_folder, "wf_inputs.txt")


    inputs_string_list = str(workflow.inputs).split("\n")

    new_wf_inputs = ""

    for input_string in inputs_string_list:

        if ("/" in input_string) and ("function_str" not in input_string):

            new_string_list = input_string.split(" = ")

            new_string = new_string_list[0] + " = "

            new_string = new_string + new_string_list[0] + "_here"

            new_wf_inputs = new_wf_inputs + "\n" + new_string

        else:

            new_wf_inputs = new_wf_inputs + "\n" + input_string

    new_wf_inputs = new_wf_inputs.replace("\n","",1)

    with open(inputs_path,"wt") as f:
        print >>f, new_wf_inputs



def run_all_tests():

    from test_anatomical_preproc import run_all_tests_anatomical_preproc
    from test_dvars import run_all_tests_dvars
    from test_functional_preproc import run_all_tests_functional_preproc
    from test_qap_workflows import run_all_tests_qap_workflows
    from test_qap_workflows_utils import run_all_tests_qap_workflows_utils
    from test_spatial_qc import run_all_tests_spatial_qc
    from test_temporal_qc import run_all_tests_temporal_qc

    run_all_tests_anatomical_preproc()
    run_all_tests_dvars()
    run_all_tests_functional_preproc()
    run_all_tests_qap_workflows()
    run_all_tests_qap_workflows_utils()
    run_all_tests_spatial_qc()
    run_all_tests_temporal_qc()

