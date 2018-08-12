import os
import pytest
import unittest
import nibabel
import shutil

import qap.functional_preproc as functional_preproc

from tests.test_utils import build_and_run_workflow


@pytest.mark.long()
class TestFunctionalPreproc(unittest.TestCase):

    def setUp(self):

        self.local_dir = os.path.dirname(os.path.abspath(__file__))
        self.working_dir = os.path.join(self.local_dir, 'workflow_output')

        self.functional_scan = os.path.join(self.local_dir, "test_data/functional_scan.nii.gz")
        self.func_reorient = os.path.join(self.local_dir, "test_data/func_reorient.nii.gz")
        self.functional_brain_mask = os.path.join(self.local_dir,
                                                  "test_data/functional_brain_mask.nii.gz")

        self.config = {"start_idx": 4,
                       "stop_idx": "End"}

        # remove any remnants of previous tests
        if os.path.exists(self.working_dir):
            shutil.rmtree(self.working_dir)

    @pytest.mark.quick
    def test_get_idx_whole_timeseries(self):

        idx_tuple = functional_preproc.get_idx(self.functional_scan, "End", 0)

        assert idx_tuple == (123, 0)

    @pytest.mark.quick
    def test_get_idx_partial_timeseries(self):

        idx_tuple = functional_preproc.get_idx(self.functional_scan, 100, 20)

        assert idx_tuple == (100, 20)

    @pytest.mark.quick
    def test_get_idx_partial_timeseries_overshoot(self):

        idx_tuple = functional_preproc.get_idx(self.functional_scan, 250, 20)

        assert idx_tuple == (123, 20)

    def test_func_preproc(self):
        """Run the 'func_preproc_workflow' function to execute the modular
        workflow with the provided inputs.
        """

        resource_pool = {'functional_scan': self.functional_scan}
        out_paths = build_and_run_workflow(functional_preproc.func_preproc_workflow,
                                           resource_pool, 'func_preproc',
                                           ['func_reorient'], config=self.config)

        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

        original_img = nibabel.load(self.functional_scan)
        preproc_img = nibabel.load(out_paths[0])
        assert(original_img.shape[3] - preproc_img.shape[3] == 4)

    def test_func_motion_correct(self):
        """Run the 'func_motion_correct_workflow' function to execute the modular
        workflow with the provided inputs.
        """

        resource_pool = {'func_reorient': self.func_reorient}
        out_paths = build_and_run_workflow(functional_preproc.func_motion_correct_workflow,
                                           resource_pool, 'func_motion_correct',
                                           ['func_motion_correct', 'func_coordinate_transformation',
                                            'func_motion_estimates'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 3)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_func_brain_mask(self):
        """ Run the 'functional_brain_mask_workflow' function to execute the modular workflow with the
            provided inputs.
        """

        resource_pool = {'func_reorient': self.func_reorient}
        out_paths = build_and_run_workflow(functional_preproc.functional_brain_mask_workflow,
                                           resource_pool, 'functional_brain_mask',
                                           ['func_brain_mask'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_invert_func_brain_mask(self):
        """ Run the 'invert_functional_brain_mask_workflow' function to execute the
            modular workflow with the provided inputs.
        """

        resource_pool = {'func_brain_mask': self.functional_brain_mask}
        out_paths = build_and_run_workflow(functional_preproc.invert_functional_brain_mask_workflow,
                                           resource_pool, 'inverted_functional_brain_mask',
                                           ['func_inverted_brain_mask'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_mean_functional(self):
        """ Run the 'mean_functional_workflow' function to execute the modular
            workflow with the provided inputs.
        """

        resource_pool = {'func_reorient': self.func_reorient}
        out_paths = build_and_run_workflow(functional_preproc.mean_functional_workflow,
                                           resource_pool, 'mean_functional',
                                           ['func_mean'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_tstd_functional(self):
        """ Run the 'mean_functional_workflow' function to execute the modular
            workflow with the provided inputs.
        """

        resource_pool = {'func_reorient': self.func_reorient}
        out_paths = build_and_run_workflow(functional_preproc.tstd_functional_workflow,
                                           resource_pool, 'tstd_functional',
                                           ['func_tstd'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_invert_mask_buildup(self):
        """
        If the invert mask cannot find its required input it should trigger the addition of a workflow
        to generate the input. This will recurse back until it no longer knows how to create the missing
        resources.
        :return:
          multiple asserts
        """

        resource_pool = {'functional_scan': self.functional_scan}
        out_paths = build_and_run_workflow(functional_preproc.invert_functional_brain_mask_workflow,
                                           resource_pool, 'inverted_functional_brain_mask',
                                           ['func_inverted_brain_mask'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 1)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_func_motion_correct_buildup(self):
        """
        If the motion correction workflow cannot find its required input it should trigger the addition of a workflow
        to generate the input. This will recurse back until it no longer knows how to create the missing
        resources.
        """

        resource_pool = {'functional_scan': self.functional_scan}
        out_paths = build_and_run_workflow(functional_preproc.func_motion_correct_workflow,
                                           resource_pool, 'func_motion_correct',
                                           ['func_motion_correct', 'func_coordinate_transformation',
                                            'func_motion_estimates'], config=self.config)

        print(out_paths)
        assert(out_paths != [] and len(out_paths) == 3)
        assert(out_paths[0].endswith(".nii.gz"))

    def test_invert_mask_buildup_empty_resource_pool(self):
        """
        If the invert mask cannot find its required input it should trigger the addition of a workflow
        to generate the input. This will recurse back until it no longer knows how to create the missing
        resources.

        This workflow builder should fail to find the necessary resources, and since it cant find them
        :return:
          multiple asserts
        """

        resource_pool = {}

        out_paths = None
        test_failed = True
        try:
            out_paths = build_and_run_workflow(functional_preproc.invert_functional_brain_mask_workflow,
                                               resource_pool, 'inverted_functional_brain_mask',
                                               ['func_inverted_brain_mask'], config=self.config)
        except ValueError:
            test_failed = False

        assert(test_failed is False)
        assert(out_paths is None)
