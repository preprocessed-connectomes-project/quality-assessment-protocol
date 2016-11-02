"""
This Python file contains the template for the QuickPACS-style workflow 
builder functions and function runners.

The runners are designed to make it easy to run each modular workflow function 
without having it connected to a larger workflow - this allows the user to 
generate only one needed output if necessary, or allows the developer to 
easily run/test the workflow function and build unit tests for it.

Each workflow function has one associated runner function.
"""

def template_workflow(workflow, resource_pool, config, name="_"):
    """Build a Nipype workflow to...

    Keyword arguments:
      workflow -- [Nipype workflow] a Nipype workflow object which can already
                  contain other connected nodes; this function will insert the
                  following workflow into this one provided
      resource_pool -- [Python dictionary] a dictionary defining input files 
                       and pointers to Nipype node outputs / workflow 
                       connections; the keys are the resource names
      config -- [Python dictionary] a dictionary defining the configuration 
                settings for the workflow, such as directory paths or toggled 
                options
      name -- [string] (default: "_") a string to append to the end of each 
              node name

    Returns:
      workflow -- [Nipype workflow] the Nipype workflow originally provided, 
                  but with the following sub-workflow connected into it
      resource_pool -- [Python dictionary] the resource pool originally 
                       provided, but updated (if applicable) with the newest 
                       outputs and connections

    Notes:
      - If any resources/outputs required by this workflow are not in the
        resource pool, this workflow will call pre-requisite workflow builder
        functions to further populate the pipeline with workflows which will
        calculate/generate these necessary pre-requisites.

    Expected Resources in Resource Pool:
      {resource name} -- {resource description}
      {resource name} -- {resource description}
      ..

    New Resources Added to Resource Pool:
      {resource name} -- {resource description}
      {resource name} -- {resource description}
      ..

    Workflow Steps:
      1. {node 1..}
      2. {node 2..}
      3. ..
    """

    import os
    import sys
    import copy
    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    ### Check resource pool for expected resources and call pre-requisite
    ### workflows if necessary
    if "{resource name}" not in resource_pool.keys():

        from module import workflow_function
        old_rp = copy.copy(resource_pool)
        workflow, new_resource_pool = \
            workflow_function(workflow, resource_pool, config, name)

        if resource_pool == old_rp:
            return workflow, resource_pool

    ### Node declarations and configuration, node connections
    node_one = pe.Node(interface=nipype_module.Interface(),
                       name='node_one%s' % name)

    node_two = pe.Node(interface=nipype_module.Interface(),
                       name='node_two%s' % name)

    ### Check if the resource is a tuple or a string
    ### - if it's a tuple: it's a Nipype workflow connection pointer
    ### - if it's a string: it's a filepath to the resource file on disk
    if len(resource_pool["{resource name}"]) == 2:
        node, out_file = resource_pool["{resource name}"]
        workflow.connect(node, out_file, node_one, 'input_file')
    else:
        node_one.inputs.in_file = resource_pool["{resource name}"]

    workflow.connect(node_one, 'output_file', node_two, 'input_file')

    ### "Send" final output to resource pool in the form of a Nipype workflow
    ### connection pointer (the tuple)
    resource_pool["{resource name}"] = (node_two, 'output_file')

    return workflow, resource_pool


def run_template_workflow(input_resource, out_dir=None, run=True):
    """Run the 'template_workflow' function to execute the modular workflow
    with the provided inputs.

    Keyword Arguments:
      input_resource -- [string] the filepath of the { input resource }
      out_dir -- [string] (default: None) the output directory to write the 
                 results to; if left as None, will write to the current 
                 directory
      run -- [boolean] (default: True) will run the workflow; if set to False,
             will connect the Nipype workflow and return the workflow object 
             instead

    Returns:
      outpath -- [string] (if run=True) the filepath of the generated output
                 file
      workflow -- [Nipype workflow] (if run=False) the Nipype workflow object
      workflow.base_dir -- [string] (if run=False) the base directory of the 
                           workflow if it were to be run
    """

    import os
    import sys

    import glob

    import nipype.interfaces.io as nio
    import nipype.pipeline.engine as pe

    output = "{ output resource name }"

    workflow = pe.Workflow(name='workflow_name')

    if not out_dir:
        out_dir = os.getcwd()

    workflow_dir = os.path.join(out_dir, "workflow_output", output)
    workflow.base_dir = workflow_dir

    resource_pool = {}
    config = {}
    num_cores_per_subject = 1

    resource_pool["input_resource"] = input_resource

    workflow, resource_pool = \
            template_workflow(workflow, resource_pool, config)

    ds = pe.Node(nio.DataSink(), name='datasink_workflow_name')
    ds.inputs.base_directory = workflow_dir

    node, out_file = resource_pool[output]

    workflow.connect(node, out_file, ds, output)

    if run == True:
        workflow.run(plugin='MultiProc', plugin_args= \
                         {'n_procs': num_cores_per_subject})
        outpath = glob.glob(os.path.join(workflow_dir, output, "*"))[0]
        return outpath
    else:
        return workflow, workflow.base_dir