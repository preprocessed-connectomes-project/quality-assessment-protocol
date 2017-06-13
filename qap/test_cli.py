
import os
import pytest
import unittest
import shutil


@pytest.mark.long()
class TestFullPipelines(unittest.TestCase):

    def setUp(self):
        # outdir stuff
        self.outdir = os.path.join(os.getcwd(), "test_data", "test_out_dir")
        self.workdir = os.path.join(self.outdir, "work_dir")
        if os.path.isdir(self.outdir):
            try:
                shutil.rmtree(self.outdir)
            except OSError:
                pass
        os.makedirs(self.outdir)
        os.makedirs(self.workdir)

        import pkg_resources as p

        self.local_bids = p.resource_filename("qap", os.path.join("test_data",
                                              "bids_data"))

        self.anat_template = p.resource_filename("qap",
                                                 os.path.join("test_data",
                                                              "MNI152_T1_3mm.nii.gz"))

    def tearDown(self):
        pass
        '''
        try:
            shutil.rmtree(self.outdir)
        except OSError:
            pass
        try:
            shutil.rmtree(self.outdir)
        except OSError:
            pass
        '''

    def test_qap_local(self):
        from qap.cli import process_args
        # just making sure it runs on the most basic input
        process_args(self.local_bids, self.outdir, "participant",
                     working_dir=self.workdir, save_working_dir=False,
                     mni_template=self.anat_template)

    @pytest.mark.skip()
    def test_qap_s3(self):
        from qap.cli import process_args
        # just making sure it runs on the most basic input
        bids_dir = "s3://fcp-indi/data/Projects/CORR/RawDataBIDS"
        process_args(bids_dir, self.outdir, "participant",
                     working_dir=self.workdir, save_working_dir=True,
                     mni_template=self.anat_template)
