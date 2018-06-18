
import pytest
import unittest


class TestPlotAll(unittest.TestCase):

    def setUp(self):
        import os
        import pandas as pd
        import pkg_resources as p
        from qap.viz.plotting import plot_all
        self.plot_all = plot_all

        anat_spat_csv = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "qap_anatomical_spatial_5rows.csv"))
        func_spat_csv = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "qap_functional_spatial_5rows.csv"))
        func_temp_csv = \
            p.resource_filename("qap", os.path.join("test_data",
                                                    "qap_functional_temporal_5rows.csv"))

        self.anat_spat_df = pd.read_csv(anat_spat_csv)
        self.func_spat_df = pd.read_csv(func_spat_csv)
        self.func_temp_df = pd.read_csv(func_temp_csv)

        self.anat_spat_sessions = \
            sorted(pd.unique(self.anat_spat_df.Session.ravel()))
        self.func_spat_sessions = \
            sorted(pd.unique(self.func_spat_df.Session.ravel()))
        self.func_temp_sessions = \
            sorted(pd.unique(self.func_temp_df.Session.ravel()))

        self.anat_spat_groups = [['CNR'],
                                ['Cortical Contrast'],
                                ['EFC'],
                                ['FBER'],
                                ['FWHM', 'FWHM_x', 'FWHM_y', 'FWHM_z'],
                                ['Qi1'],
                                ['SNR']]

        self.func_spat_groups = [['EFC'],
                                ['FBER'],
                                ['FWHM', 'FWHM_x', 'FWHM_y', 'FWHM_z'],
                                ['Ghost_%s' % a for a in ['x', 'y', 'z']],
                                ['SNR']]

        self.func_temp_groups = [['Fraction of Outliers (Mean)',
                                  'Fraction of Outliers (Median)',
                                  'Fraction of Outliers (Std Dev)',
                                  'Fraction of Outliers IQR'],
                                 ['GCOR'],
                                 ['Quality (Mean)', 'Quality (Median)',
                                  'Quality (Std Dev)', 'Quality IQR',
                                  'Quality percent outliers'],
                                 ['RMSD (Mean)', 'RMSD (Median)',
                                  'RMSD (Std Dev)', 'RMSD IQR'],
                                 ['Std. DVARS (Mean)', 'Std. DVARS (Median)',
                                  'Std. DVARS percent outliers',
                                  'Std. DVARs IQR']]

    def test_plot_all_anat_spat(self):
        df = self.anat_spat_df.copy()
        for ss in self.anat_spat_sessions:
            sesdf = df.copy().loc[df['Session'] == ss]
            subtitle = '(%s)' % (ss)
            sub_id = None
            fig = self.plot_all(sesdf, self.anat_spat_groups, subject=sub_id,
                                title='QC measures ' + subtitle)

    def test_plot_all_func_spat(self):
        df = self.func_spat_df.copy()
        for ss in self.func_spat_sessions:
            sesdf = df.copy().loc[df['Session'] == ss]
            subtitle = '(%s)' % (ss)
            sub_id = None
            fig = self.plot_all(sesdf, self.func_spat_groups, subject=sub_id,
                                title='QC measures ' + subtitle)

    def test_plot_all_func_temp(self):
        df = self.func_temp_df.copy()
        for ss in self.func_temp_sessions:
            sesdf = df.copy().loc[df['Session'] == ss]
            subtitle = '(%s)' % (ss)
            sub_id = None
            fig = self.plot_all(sesdf, self.func_temp_groups, subject=sub_id,
                                title='QC measures ' + subtitle)

    def test_plot_all_anat_spat_onesub(self):
        df = self.anat_spat_df.copy()
        for ss in self.anat_spat_sessions:
            sesdf = df.copy().loc[df['Session'] == ss]
            subtitle = '(Participant %s_%s)' % ("0003001", ss)
            sub_id = "0003001"
            fig = self.plot_all(sesdf, self.anat_spat_groups, subject=sub_id,
                                title='QC measures ' + subtitle)


class TestCalcFD(unittest.TestCase):

    def setUp(self):
        import os
        import pkg_resources as p
        from qap.viz.plotting import _calc_fd
        self.calc_fd = _calc_fd

        self.fd_file = p.resource_filename("qap", os.path.join("test_data",
                                                               "meanFD.1D"))

    def test_fd_calculation(self):
        ref_fd = [0., 0.15495802, 0.01828499, 0.0494406, 0.00662847]
        fd_power = self.calc_fd(self.fd_file)
        for idx in range(0, len(ref_fd)):
            self.assertAlmostEqual(ref_fd[idx], fd_power[idx])
