
import pytest
import unittest


class TestGatherFilepathList(unittest.TestCase):

    def setUp(self):
        # setup
        import os
        import pkg_resources as p
        from qap.script_utils import gather_filepath_list
        self.gather_filepath_list = gather_filepath_list

        # inputs
        self.data_folder = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "data_folder"))

        # outputs
        self.ref_path_list = [
            "site_1/sub_01/ses_01/anat_1/anatomical_scan.nii.gz",
            "site_1/sub_01/ses_01/rest_1/functional_scan.nii.gz"]

    def test_custom_filepaths(self):
        test_path_list = self.gather_filepath_list(self.data_folder)
        self.assertListEqual(self.ref_path_list, test_path_list)


@pytest.mark.long
class TestPullS3Sublist(unittest.TestCase):

    # will fail if no internet connection
    # use this test fixture periodically

    def setUp(self):
        # setup
        from qap.script_utils import pull_s3_sublist
        self.pull_s3_sublist = pull_s3_sublist

        # inputs
        self.bids_path = "s3://fcp-indi/data/Projects/CORR/RawDataBIDS"
        self.custom_path = "s3://fcp-indi/data/Projects/CORR/RawData"
        self.invalid_bucket_path = "s3://fcp--indi/data/Projects/CORR/RawDataBIDS"
        self.invalid_dir_path = "s3://fcp-indi/data/Projects/corr/RawDataBIDS"

        # outputs
        self.bids_s3_list = [
            'BMB_1/T1w.json',
            'BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz',
            'BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz',
            'BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz',
            'BMB_1/sub-0003001/sub-0003001_sessions.tsv']
        self.custom_s3_list = [
            'BMB_1/0003001/session_1/anat_1/anat.nii.gz',
            'BMB_1/0003001/session_1/rest_1/rest.nii.gz',
            'BMB_1/0003001/session_1/rest_2/rest.nii.gz',
            'BMB_1/0003002/session_1/anat_1/anat.nii.gz',
            'BMB_1/0003002/session_1/rest_1/rest.nii.gz']

    def test_BIDS(self):
        test_bids_s3_list = self.pull_s3_sublist(self.bids_path)
        self.assertListEqual(self.bids_s3_list, test_bids_s3_list[0][0:5])

    def test_custom(self):
        test_custom_s3_list = self.pull_s3_sublist(self.custom_path)
        self.assertListEqual(self.custom_s3_list, test_custom_s3_list[0][0:5])

    def test_invalid_bucket_name(self):
        with self.assertRaises(Exception):
            self.pull_s3_sublist(self.invalid_bucket_path)

    def test_wrong_dirpath(self):
        test_wrong_list = self.pull_s3_sublist(self.invalid_dir_path)
        self.assertEquals(0, len(test_wrong_list[0]))

    def test_invalid_creds_path(self):
        with self.assertRaises(Exception):
            self.pull_s3_sublist(self.bids_path, "/path/to/nothing.csv")


class TestParseRawDataList(unittest.TestCase):

    # for non-BIDS data directory formats

    def setUp(self):
        # setup
        from qap.script_utils import parse_raw_data_list
        self.parse_raw_data_list = parse_raw_data_list

        # inputs
        self.local_data_folder = "/data/dir"
        self.local_file_list = ["site_1/sub_01/ses_01/anat_1/mprage.nii.gz",
                                "site_1/sub_02/ses_01/func_1/rest.nii.gz"]
        self.wrong_file_list = ["site_1/sub_01/anat_1/mprage.nii.gz"]
        self.s3_data_folder = "s3://data/Projects/RawData"
        self.s3_file_list = ["site_1/sub_01/ses_01/anat_1/mprage.nii.gz",
                             "site_1/sub_02/ses_01/func_1/rest.nii.gz"]

        # outputs
        self.ref_local_subdict = {
            'sub_01': {
                'ses_01': {
                    'anatomical_scan': {
                        'anat_1': '/data/dir/site_1/sub_01/ses_01/anat_1/mprage.nii.gz'},
                    'site_name': 'site_1'}},
            'sub_02': {
                'ses_01': {
                    'functional_scan': {
                        'func_1': '/data/dir/site_1/sub_02/ses_01/func_1/rest.nii.gz'},
                    'site_name': 'site_1'}}}

    def test_local_filepaths(self):
        test_local = self.parse_raw_data_list(self.local_file_list,
                                              self.local_data_folder)
        self.assertDictEqual(self.ref_local_subdict, test_local)

    def test_s3_filepaths(self):
        # TODO
        pass

    def test_inclusion(self):
        # TODO
        pass

    def test_wrong_dir_format(self):
        with self.assertRaises(IndexError):
            self.parse_raw_data_list(self.wrong_file_list,
                                     self.local_data_folder)


@pytest.mark.quick
def test_gather_custom_raw_data():

    from qap.script_utils import gather_custom_raw_data

    # we are starting in the directory containing the site folders!
    site_folder = "/home/data"

    format = "/{site}/{participant}/{session}/{series}"

    anatomical_keywords = "mprage"
    functional_keywords = "rest func"

    filepath_list = [
      "/home/data/site01/sub01/sess01/anat_1/mprage.nii.gz",
      "/home/data/site01/sub01/sess02/anat_1/mprage.nii.gz",
      "/home/data/site01/sub02/sess01/anat_1/mprage.nii.gz",
      "/home/data/site01/sub02/sess02/anat_1/mprage.nii.gz",
      "/home/data/site01/sub01/sess01/rest_1/rest.nii.gz",
      "/home/data/site01/sub01/sess01/rest_2/rest.nii.gz",
      "/home/data/site01/sub01/sess02/rest_1/func.nii.gz",
      "/home/data/site01/sub02/sess01/rest_1/rest.nii.gz",
      "/home/data/site01/sub02/sess01/rest_2/rest.nii.gz",
      "/home/data/site01/sub02/sess02/rest_1/func.nii.gz",
    ]

    # include sites
    ref_sub_dict = {
      'sub01': {'sess01': {'anatomical_scan': {'anat_1': '/home/data/site01/sub01/sess01/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub01/sess01/rest_1/rest.nii.gz',
                                               'rest_2': '/home/data/site01/sub01/sess01/rest_2/rest.nii.gz'},
                           'site_name': 'site01'},
                'sess02': {'anatomical_scan': {'anat_1': '/home/data/site01/sub01/sess02/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub01/sess02/rest_1/func.nii.gz'},
                           'site_name': 'site01'}},
      'sub02': {'sess01': {'anatomical_scan': {'anat_1': '/home/data/site01/sub02/sess01/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub02/sess01/rest_1/rest.nii.gz',
                                               'rest_2': '/home/data/site01/sub02/sess01/rest_2/rest.nii.gz'},
                           'site_name': 'site01'},
                'sess02': {'anatomical_scan': {'anat_1': '/home/data/site01/sub02/sess02/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub02/sess02/rest_1/func.nii.gz'},
                           'site_name': 'site01'}}}

    sub_dict = gather_custom_raw_data(filepath_list, site_folder, format, 
        anatomical_keywords, functional_keywords)

    assert ref_sub_dict == sub_dict


@pytest.mark.quick
def test_gather_custom_raw_data_scans_folder():

    from qap.script_utils import gather_custom_raw_data

    # we are starting in the directory containing the site folders!
    site_folder = "/home/data"

    format = "/{site}/{participant}/{session}/scans/{series}"

    anatomical_keywords = "mprage"
    functional_keywords = "rest func"

    # inclusion of a "scans" folder in between the session and scan folders
    filepath_list = [
      "/home/data/site01/sub01/sess01/scans/anat_1/mprage.nii.gz",
      "/home/data/site01/sub01/sess02/scans/anat_1/mprage.nii.gz",
      "/home/data/site01/sub02/sess01/scans/anat_1/mprage.nii.gz",
      "/home/data/site01/sub02/sess02/scans/anat_1/mprage.nii.gz",
      "/home/data/site01/sub01/sess01/scans/rest_1/rest.nii.gz",
      "/home/data/site01/sub01/sess01/scans/rest_2/rest.nii.gz",
      "/home/data/site01/sub01/sess02/scans/rest_1/func.nii.gz",
      "/home/data/site01/sub02/sess01/scans/rest_1/rest.nii.gz",
      "/home/data/site01/sub02/sess01/scans/rest_2/rest.nii.gz",
      "/home/data/site01/sub02/sess02/scans/rest_1/func.nii.gz",
    ]

    # include sites
    ref_sub_dict = {
      'sub01': {'sess01': {'anatomical_scan': {'anat_1': '/home/data/site01/sub01/sess01/scans/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub01/sess01/scans/rest_1/rest.nii.gz',
                                               'rest_2': '/home/data/site01/sub01/sess01/scans/rest_2/rest.nii.gz'},
                           'site_name': 'site01'},
                'sess02': {'anatomical_scan': {'anat_1': '/home/data/site01/sub01/sess02/scans/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub01/sess02/scans/rest_1/func.nii.gz'},
                           'site_name': 'site01'}},
      'sub02': {'sess01': {'anatomical_scan': {'anat_1': '/home/data/site01/sub02/sess01/scans/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub02/sess01/scans/rest_1/rest.nii.gz',
                                               'rest_2': '/home/data/site01/sub02/sess01/scans/rest_2/rest.nii.gz'},
                           'site_name': 'site01'},
                'sess02': {'anatomical_scan': {'anat_1': '/home/data/site01/sub02/sess02/scans/anat_1/mprage.nii.gz'},
                           'functional_scan': {'rest_1': '/home/data/site01/sub02/sess02/scans/rest_1/func.nii.gz'},
                           'site_name': 'site01'}}}

    sub_dict = gather_custom_raw_data(filepath_list, site_folder, format, 
        anatomical_keywords, functional_keywords)

    assert ref_sub_dict == sub_dict
