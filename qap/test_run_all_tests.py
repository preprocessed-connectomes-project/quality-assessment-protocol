
def run_all_qap_tests():

    from qap import test_anatomical_preproc
    from qap import test_functional_preproc
    from qap import test_qap_workflows
    from qap import test_spatial_qc
    from qap import test_temporal_qc
    from qap import test_dvars

    test_anatomical_preproc.run_all_tests_anatomical_preproc()
    test_functional_preproc.run_all_tests_functional_preproc()
    test_qap_workflows.run_all_tests_qap_workflows()
    test_spatial_qc.run_all_tests_spatial_qc()
    test_temporal_qc.run_all_tests_temporal_qc()
    test_dvars.run_all_tests_dvars()