import os
import glob
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe


def build_and_run_workflow(workflow_builder, resource_pool, workflow_name, output_keys, config=None):
    """Run a workflow builder function to execute the modular
    workflow with the inputs provided in resource_pool.

    :param workflow_builder: a workflow builder function

    :type resource_pool: dict
    :param resource_pool: dictionary of resources, files from other steps of the pipeline
    :type workflow_name: string
    :param workflow_name: name to use for the pipeline and output and working directories
    :type output_keys: list
    :param output_keys: list of outputs from the workflow that should be saved to the output directory
    :type config: dict
    :param config: dictionary of extra arguments passed through to workflow builder


    :rtype: Nipype workflow object
    :return: the object returned by nipype after the workflow executed
    """

    workflow = pe.Workflow(name='{0}_workflow'.format(workflow_name))
    workflow.base_dir = os.path.join(os.getcwd(), "workflow_output", workflow_name)

    workflow, resource_pool = workflow_builder(workflow, resource_pool, config)

    print("resource pool: {0}".format(resource_pool))

    for output_key in output_keys:
        ds = pe.Node(nio.DataSink(), name='datasink_{0}_{1}'.format(workflow_name, output_key))
        ds.inputs.base_directory = workflow.base_dir
        node, out_file = resource_pool[output_key]
        workflow.connect(node, out_file, ds, output_key)

    num_cores = 2
    workflow.run(plugin='MultiProc', plugin_args={'n_procs': num_cores})

    outpaths = []
    for output_key in output_keys:
        outpaths += glob.glob(os.path.join(workflow.base_dir, output_key, "*"))

    return outpaths


def get_test_dir(key):

    import datetime
    import os

    time_stamp_string = datetime.datetime.strftime(datetime.datetime.now(), "%Y%m%d%H%M%S")
    working_dir = os.path.join('/tmp', '{0}_{1}'.format(key, time_stamp_string))
    os.makedirs(working_dir)

    return working_dir
