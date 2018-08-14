
import os
import unittest
import pytest

import qap.bids_utils as bids_utils


@pytest.mark.long()
class TestBidsUtils(unittest.TestCase):

    def setUp(self):
        self.local_dir = os.path.dirname(os.path.abspath(__file__))
        self.local_bids_dir = \
            os.path.join(self.local_dir, 'test_data/bids_data')
        self.s3_bids_dir = "s3://fcp-indi/data/Projects/CORR/RawDataBIDS/IBA_TRT"

    def test_bids_decode_filename(self):

        decoded = bids_utils.bids_decode_filename(
            'BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz')

        assert decoded['site'] == "BMB1"
        assert decoded['sub'] == "0003001"
        assert decoded['ses'] == "1"
        assert decoded['run'] == "2"
        assert decoded['task'] == "rest"
        assert decoded['scantype'] == "bold"

    def test_bids_retrieve_parameters(self):
        from_filename = bids_utils.bids_decode_filename(
            'BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz')

        parameters = bids_utils.bids_retrieve_parameters({
            "scantype-bold": {
                "site-BMB1": {
                    "sub-0003001": {
                        "ses-1": {
                            "task-rest": {
                                "acq-none": {
                                    "rec-none": {
                                        "run-2": {
                                            "RepetitionTime": 2
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }, from_filename)

        assert parameters["RepetitionTime"] == 2

    def test_bids_decode_s3(self):
        bids_dir = self.s3_bids_dir
        bids_files, derivative_files, _ = bids_utils.collect_bids_files_configs(
            bids_dir, [])
        bids_files = bids_files + derivative_files
        assert bids_files is not None

    def test_collect_bids_files_configs(self):
        # just make sure it runs
        from qap.bids_utils import collect_bids_files_configs
        s3_bids_dir = "s3://fcp-indi/data/Projects/CORR/RawDataBIDS"
        file_paths, derivative_files, config_dict = collect_bids_files_configs(
            s3_bids_dir)
        assert file_paths is not None
        assert derivative_files == []
        assert config_dict is not None


class TestExtractBidsData(unittest.TestCase):

    def setUp(self):
        # inputs
        self.input_s3_bids_dir = "s3://fcp-indi/data/Projects/CORR/RawDataBIDS"
        self.input_BIDS_CoRR = [
            'BMB_1/T1w.json',
            'BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz',
            'BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz',
            'BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz',
            'BMB_1/sub-0003001/sub-0003001_sessions.tsv',
            'BMB_1/sub-0003002/ses-1/anat/sub-0003002_ses-1_run-1_T1w.nii.gz',
            'BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz',
            'BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz',
            'BMB_1/sub-0003002/sub-0003002_sessions.tsv'
        ]

        # expected outputs
        self.corr_ref_s3_dict = {
            'site-BMB1': {
                'sub-0003001': {
                    'ses-1': {
                        'run-1_T1w': {
                            'anatomical_scan':
                                's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/anat'
                                '/sub-0003001_ses-1_run-1_T1w.nii.gz'},
                        'task-rest_run-1_bold': {
                            'functional_scan':
                                's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func'
                                '/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz'},
                        'task-rest_run-2_bold': {
                            'functional_scan':
                                's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func'
                                '/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz'}}},
                'sub-0003002': {
                    'ses-1': {
                        'run-1_T1w': {
                            'anatomical_scan':
                                's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/anat'
                                '/sub-0003002_ses-1_run-1_T1w.nii.gz'},
                        'task-rest_run-1_bold': {
                            'functional_scan':
                                's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func'
                                '/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz'},
                        'task-rest_run-2_bold': {
                            'functional_scan':
                                's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func'
                                '/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz'
                        }
                    }
                }
            }
        }

    def test_s3_paths(self):
        test_dict = bids_utils.bids_generate_qap_data_configuration(
            self.input_s3_bids_dir, self.input_BIDS_CoRR)
        self.assertDictEqual(self.corr_ref_s3_dict, test_dict)
