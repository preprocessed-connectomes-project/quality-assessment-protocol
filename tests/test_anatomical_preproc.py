import os
import pytest
import unittest
import shutil
import qap.anatomical_preproc as anatomical_preproc

from tests.test_utils import build_and_run_workflow


@pytest.mark.long()
class TestAnatomicalPreproc(unittest.TestCase):

    def setUp(self):

        self.local_dir = os.path.dirname(os.path.abspath(__file__))
        self.working_dir = os.path.join(self.local_dir, 'workflow_output')

        self.anatomical_scan = os.path.join(
            self.local_dir, "test_data/anatomical_scan.nii.gz")
        self.anat_brain = os.path.join(
            self.local_dir, "test_data/anat_brain.nii.gz")
        self.anat_reorient = os.path.join(
            self.local_dir, "test_data/anat_reorient.nii.gz")

        self.config = {
            "skull_on_registration": True,
            "anatomical_template": os.path.join(self.local_dir, "test_data/MNI152_T1_3mm.nii.gz"),
            "execution": {
                "plugin": "Linear",
                "plugin_args": {},
                "working_dir": self.working_dir
            }
        }

        # remove any remnants of previous tests
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    def test_anatomical_reorient(self):
        """

        test the anatomical reorient workflow

        :return: asserts on error
        """

        resource_pool = {'anatomical_scan': self.anatomical_scan}
        out_paths = build_and_run_workflow(anatomical_preproc.anatomical_reorient_workflow,
                                           resource_pool, 'anatomical_reorient',
                                           ['anat_reorient'], config=self.config)

        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_afni_segmentation(self):
        """Run the 'afni_segmentation_workflow' function to execute the modular
        workflow with the provided inputs.
        """

        resource_pool = {'anat_brain': self.anat_brain}
        out_paths = build_and_run_workflow(anatomical_preproc.afni_segmentation_workflow,
                                           resource_pool, 'afni_segmentation',
                                           ['anat_csf_mask', 'anat_gm_mask',
                                               'anat_wm_mask'],
                                           config=self.config)

        assert(out_paths != [] and len(out_paths) == 3)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_afni_segmentation_buildup(self):
        """Run the 'afni_segmentation_workflow' function to execute the modular
        workflow with the provided inputs.
        """

        resource_pool = {'anatomical_scan': self.anatomical_scan}
        out_paths = build_and_run_workflow(anatomical_preproc.afni_segmentation_workflow,
                                           resource_pool, 'afni_segmentation',
                                           ['anat_csf_mask', 'anat_gm_mask',
                                               'anat_wm_mask'],
                                           config=self.config)

        assert(out_paths != [] and len(out_paths) == 3)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_afni_anatomical_linear_registration(self):
        """Run the 'afni_anatomical_linear_registration' function to execute the
        modular workflow with the provided inputs."""
        resource_pool = {'anat_reorient': self.anat_reorient}
        out_paths = build_and_run_workflow(anatomical_preproc.afni_anatomical_linear_registration,
                                           resource_pool, 'afni_linear_registration',
                                           ['anat_linear_xfm',
                                               'anat_linear_warped_anat_reorient'],
                                           config=self.config)

        assert(out_paths != [] and len(out_paths) == 2)
        assert(out_paths[0].endswith(".aff12.1D")
               and out_paths[1].endswith(".nii.gz"))

    def test_anatomical_skullstrip(self):
        """Run the 'anatomical_skullstrip_workflow' function to execute the
        modular workflow with the provided inputs."""
        resource_pool = {'anat_reorient': self.anat_reorient}
        out_paths = build_and_run_workflow(anatomical_preproc.anatomical_skullstrip_workflow,
                                           resource_pool, 'anatomical_skullstrip',
                                           ['anat_brain'],
                                           config=self.config)

        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_afni_segmentation_buildup_noresource(self):
        """Run the 'afni_segmentation_workflow' function to execute the modular
        workflow with the provided inputs.
        """

        out_paths = None
        resource_pool = {}
        test_failed = True
        try:
            out_paths = build_and_run_workflow(anatomical_preproc.afni_segmentation_workflow,
                                               resource_pool, 'afni_segmentation',
                                               ['anat_csf_mask', 'anat_gm_mask',
                                                   'anat_wm_mask'],
                                               config=self.config)

        except ValueError:
            test_failed = False

        assert(test_failed is False)
        assert(out_paths is None)
