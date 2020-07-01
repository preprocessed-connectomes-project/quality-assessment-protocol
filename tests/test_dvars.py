import unittest
test_sub_dir = "test_data"


class TestStandardizedDVARS(unittest.TestCase):

    def setUp(self):
        import os
        import pkg_resources as p

        self.func_reorient = \
            p.resource_filename("tests", os.path.join(test_sub_dir, "func_reorient.nii.gz"))
        self.func_mask = \
            p.resource_filename("tests", os.path.join(test_sub_dir, "fsl_bet_mask.nii.gz"))

        self.rdvars_reference = p.resource_filename("tests", os.path.join(test_sub_dir, "DVARSout_nomask_noz"))
        self.dse_reference = p.resource_filename("tests", os.path.join(test_sub_dir, "dvars_matlab_result.json"))

    def test_calc_dvars(self):
        import json
        import numpy as np
        import numpy.testing as nt
        from qap.dvars import calc_dse

        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        dse = calc_dse(self.func_reorient, self.func_mask)

        with open("/tmp/dse_test_out.json", 'w') as ofd:
            json.dump(dse, ofd, indent=2, cls=NumpyEncoder)

        ref_dvars = np.loadtxt(self.rdvars_reference)
        with open(self.dse_reference, 'r') as infd:
            ref_dse = json.load(infd)

        nt.assert_array_almost_equal(dse['DSE_Relative_DVARS'], ref_dvars, decimal=5)

        for python_name, matlab_name in [('DSE_Total_Variance', 'Avar'),
                                         ('DSE_Delta_Percent_DVAR', 'DeltapDvar'),
                                         ('DSE_P_Values', 'pvals')]:
            nt.assert_array_almost_equal(dse[python_name], ref_dse[matlab_name])
