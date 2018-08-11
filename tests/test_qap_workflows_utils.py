
import pytest
import unittest


class TestQapFunctionalSpatial(unittest.TestCase):

    def setUp(self):
        # init
        import os
        import pickle
        import pkg_resources as p
        from qap.qap_workflows_utils import qap_functional_spatial as qfs
        self.qfs = qfs

        # inputs
        self.func_mean = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "mean_functional.nii.gz"))
        self.func_mask = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "functional_brain_mask.nii.gz"))

    def test_qap_measures_dict(self):
        qap_dict = self.qfs(self.func_mean, self.func_mask, "y",
                            "sub1", "ses1", "scan1", "site1")
        partic = qap_dict['sub1 ses1 scan1']['Participant']
        site = qap_dict['sub1 ses1 scan1']['Site']
        ghost_y = qap_dict['sub1 ses1 scan1']['functional_spatial']['Ghost_y']
        self.assertEqual('sub1', partic)
        self.assertEqual('site1', site)
        #self.assertAlmostEqual(0.03716297338, float(qual_mean))
        print ghost_y


class TestQapFunctionalTemporal(unittest.TestCase):

    def setUp(self):
        # init
        import os
        import pickle
        import pkg_resources as p
        from qap.qap_workflows_utils import qap_functional_temporal as qft
        self.qft = qft

        # inputs
        self.func_ts = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "func_reorient.nii.gz"))
        self.func_mask = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "functional_brain_mask.nii.gz"))
        self.func_bg_mask = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "inverted_functional_brain_mask.nii.gz"))
        self.rmsd_file = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "meanFD.1D"))

        # outputs
        self.ref_motion_param_ts = None

        with open(self.rmsd_file, "r") as f:
            fd_lines = f.readlines()
        self.fd_floats = [float(x) for x in fd_lines]

    def test_qap_measures_dict(self):
        qap_dict, qa_dict = self.qft(self.func_ts, self.func_mask,
                                     self.func_bg_mask, self.rmsd_file,
                                     "sub1", "ses1", "scan1", "site1")
        partic = qap_dict['sub1 ses1 scan1']['Participant']
        site = qap_dict['sub1 ses1 scan1']['Site']
        qual_mean = qap_dict['sub1 ses1 scan1']['functional_temporal']['Quality (Mean)']
        self.assertEqual('sub1', partic)
        self.assertEqual('site1', site)
        self.assertAlmostEqual(0.03716297338, float(qual_mean))

    def test_qa_dict(self):
        qap_dict, qa_dict = self.qft(self.func_ts, self.func_mask,
                                     self.func_bg_mask, self.rmsd_file,
                                     "sub1", "ses1", "scan1", "site1")
        fd = qa_dict['sub1 ses1 scan1']['RMSD']
        self.assertListEqual(self.fd_floats, fd)
