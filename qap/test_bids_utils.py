
import pytest
import unittest


class TestExtractBidsData(unittest.TestCase):

    def setUp(self):
        # setup
        from qap.bids_utils import bids_gen_qap_sublist
        self.bids_gen_qap_sublist = bids_gen_qap_sublist
        self.maxDiff = None

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
            'BMB_1/sub-0003002/sub-0003002_sessions.tsv']

        # outputs
        self.corr_ref_s3_dict = {
            'sub-0003001': {
                'ses-1': {
                    'anatomical_scan': {
                        'run-1_T1w': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz'},
                    'functional_scan': {
                        'task-rest_run-1': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz',
                        'task-rest_run-2': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz'},
                    'creds_path': '', 'site_name': 'site-BMB1'}},
            'sub-0003002': {
                'ses-1': {
                    'anatomical_scan': {
                        'run-1_T1w': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/anat/sub-0003002_ses-1_run-1_T1w.nii.gz'},
                    'functional_scan': {
                        'task-rest_run-1': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz',
                        'task-rest_run-2': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz'},
                    'creds_path': '', 'site_name': 'site-BMB1'}}}

    def test_s3_paths(self):
        test_subdict = self.bids_gen_qap_sublist(self.input_s3_bids_dir,
                                                 self.input_BIDS_CoRR)
        self.assertDictEqual(self.corr_ref_s3_dict, test_subdict)

    def test_s3_paths_creds(self):
        # TODO
        pass

