
import os
import unittest

import pytest

import qap.qap_cfg as qap_cfg


@pytest.mark.long()
class TestQapCfg(unittest.TestCase):

    def setUp(self):
        self.local_dir = os.path.dirname(os.path.abspath(__file__))
        self.local_bids_dir = os.path.join(self.local_dir, 'test_data/bids_data')
        self.s3_bids_dir = "s3://fcp-indi/data/Projects/CORR/RawDataBIDS/IBA_TRT"

    def test_qap_cfg(self):
        qap_cfg.write_pipeline_configuration("test_config.yml", qap_cfg.default_pipeline_configuration)

        print(qap_cfg.configuration_output_string.format(**qap_cfg.default_pipeline_configuration))

        import yaml

        test_configuration = yaml.load(open("test_config.yml", "r"))
        print(qap_cfg.configuration_output_string.format(**test_configuration))

        qap_cfg.validate_pipeline_configuration(test_configuration)
        qap_cfg.validate_pipeline_configuration(qap_cfg.default_pipeline_configuration)

        assert (test_configuration != qap_cfg.default_pipeline_configuration)