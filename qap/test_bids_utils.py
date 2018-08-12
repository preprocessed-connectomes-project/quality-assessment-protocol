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
        self.local_bids_dir = os.path.join(self.local_dir, 'test_data/bids_data')


def test_bids_decode(bids_dir, dbg=True):
    (bids_files, derivative_files, config) = collect_bids_files_configs(bids_dir, [])
    if dbg:
        print("Found %d config files for %d image files and %d derivatives" % (len(config), len(bids_files),
                                                                               len(derivative_files)))
    bids_files = bids_files + derivative_files

    for bids_file in bids_files:
        print("{0} :: {1}".format(bids_file, bids_decode_filename(bids_file)))

    return


def test_gen_bids_sublist_qap(bids_dir, test_yml, aws_input_credentials=None, debug=False, cfg=False):
    (bids_files, derivative_files, config) = collect_bids_files_configs(bids_dir,
                                                                        aws_input_credentials=aws_input_credentials)
    if debug:
        num_json = 0
        for df in derivative_files:
            if df.endswith("json"):
                num_json += 1
        print("Found {0} config files for {1} image files and {2} derivatives ({3} json derivatives)".format(
            len(config), len(bids_files), len(derivative_files), num_json))

    bids_files = bids_files + derivative_files

    if cfg:
        data_configuration = bids_generate_qap_data_configuration(bids_dir, bids_files, config,
                                                                  credentials_path=aws_input_credentials, debug=debug)
    else:
        data_configuration = bids_generate_qap_data_configuration(bids_dir, bids_files,
                                                                  credentials_path=aws_input_credentials, debug=debug)

    with open(test_yml, "w") as ofd:
        yaml.dump(data_configuration, ofd, encoding='utf-8')

    assert data_configuration


if __name__ == '__main__':
    test_gen_bids_sublist_qap("/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/bids_dir",
                              "/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/qap_test.yml",
                              debug=True)
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/data/ADHD200/RawDataBIDS/",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/test/rs_subject_list.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)

    # test_bids_decode("/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/bids_dir")

    #
    # test_gen_bids_sublist(
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/data/ADHD200/RawDataBIDS/",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/test/rs_subject_list.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
    # test_gen_bids_sublist(
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/data/ADHD200/RawDataBIDS/Peking_3",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC"
    #     "/test/rs_subject_list_pk3.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
    test_gen_bids_sublist_qap(
        "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/",
        "/Users/cameron.craddock/workspace/git_temp/qap_v1.9/qap_test_data/rs_subject_list_s3.yml",
        aws_input_credentials="/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
        debug=False)
    #
    # test_gen_bids_sublist(
    #     "s3://fcp-indi/data/Projects/ADHD200/RawDataBIDS/Peking_3",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC/test/"
    #        "rs_subject_list_pk3_s3.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
    # test_gen_bids_sublist(
    #     "s3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1",
    #     "/Users/cameron.craddock/workspace/git_temp/CPAC/test/"
    #        "rs_subject_list_corr_bmb3_s3.yml",
    #     "/Users/cameron.craddock/AWS/ccraddock-fcp-indi-keys2.csv",
    #     dbg=False)
    #
