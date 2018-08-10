
import unittest
import pytest
import os
import glob
import pkg_resources

import nipype.interfaces.io as nipype_io
import nipype.pipeline.engine as nipype_pipe_engine

import qap_workflows

test_sub_dir = "test_data"


def setup_base_workflow(out_dir):

    import os
    import nipype.pipeline.engine as pe
    import nipype.interfaces.utility as util

    workflow = pe.Workflow(name='workflow_name')

    workflow_dir = os.path.join(out_dir, "workflow_output")
    workflow.base_dir = workflow_dir

    wf_resource_pool = {}

    for root, dirs, files in os.walk(out_dir, topdown=False):
        if files:
            for filename in files:
                if filename.endswith(".nii") or filename.endswith(".nii.gz"):
                    file_path = os.path.join(root, filename)  # type: str
                    resource_name = filename.split("_")[-1]
                    if "." in resource_name:
                        resource_name = resource_name.split(".")[0]
                    resource_name = resource_name.replace("-", "_")
                    wf_resource_pool[resource_name] = file_path

    # push them all into identity interfaces
    for resource in wf_resource_pool:
        inputnode = pe.Node(util.IdentityInterface(
                                    fields=['input'],
                                    mandatory_inputs=True),
                            name='inputnode_%s' % resource)
        inputnode.inputs.input = wf_resource_pool[resource]
        workflow.add_nodes([inputnode])
        wf_resource_pool[resource] = (inputnode, 'input')

    return workflow, wf_resource_pool


@pytest.mark.skip()
class TestQAPGatherHeaderInfo(unittest.TestCase):

    def setUp(self):
        import os
        import pkg_resources as p
        from qap.qap_workflows import qap_gather_header_info

        self.test_data = \
            p.resource_filename("qap", os.path.join("test_data", "bids_data"))

        self.wf, self.rp = setup_base_workflow(self.test_data)
        self.qghi = qap_gather_header_info

        self.config = {"subject_id": "sub01", "session_id": "ses01",
                       "scan_id": "func01", "site_name": "Caltech", "run_name": "unit_test",
                       "output_directory": self.test_data}

        self.name = "_".join(["", self.config["site_name"], self.config["subject_id"], self.config["session_id"],
                              self.config["scan_id"]])

    def test_header_extraction(self):
        print("\nresource pool: {0}".format(self.rp))
        wf, rp = self.qghi(self.wf, self.rp, self.config, name=self.name)
        wf.run()


def test_qap_mask_wf():

    anatomical_reorient = \
        pkg_resources.resource_filename("qap",
                                        os.path.join("test_data",
                                                     "bids_data",
                                                     "site-1_sub-1_ses-1",
                                                     "anat",
                                                     "site-1_sub-1_ses-1_run-1_anatomical-reorient.nii.gz"))

    allineate_linear_xfm = \
        pkg_resources.resource_filename("qap",
                                        os.path.join("test_data", "bids_data", "site-1_sub-1_ses-1", "anat",
                                                     "site-1_sub-1_ses-1_run-1_allineate-linear-xfm.aff12.1D"))

    resource_pool = {"anatomical_reorient": anatomical_reorient,
                     "allineate_linear_xfm": allineate_linear_xfm}

    base_dir = os.path.join(os.getcwd(), 'qap_mask_workflow')

    wf = nipype_pipe_engine.Workflow(name='nipype_workflow')
    wf.base_dir = os.path.join(base_dir, "working_dir")

    wf, resource_pool = qap_workflows.qap_mask_workflow(wf, resource_pool, config={})

    for key in resource_pool.keys():
        ds = nipype_pipe_engine.Node(nipype_io.DataSink(), name='datasink_{0}'.format(key))
        ds.inputs.base_directory = os.path.join(base_dir, "output")
        if isinstance(resource_pool[key], tuple):
            node, out_file = resource_pool[key]
            wf.connect(node, out_file, ds, key)

    config_dict = {"num_processors":2}

    if "num_processors" not in config_dict.keys():
        config_dict["num_processors"] = 1
    retVal = wf.run(plugin='MultiProc',
                    plugin_args={'n_procs': config_dict["num_processors"]})
    out_paths = glob.glob(os.path.join(base_dir, "output", "*"))


    assert(len(out_paths) > 0)
