
import pytest
import unittest


# un-skip this once the CLI is un-classed
@pytest.mark.skip()
class TestValidateConfigDict(unittest.TestCase):

    def setUp(self):
        from qap.cli import validate_config_dict
        self.validate_config_dict = validate_config_dict
        self.good_config_dict = {"num_processors": 4,
                                 "output_directory": "/path/to/output"}
        self.bad_config_dict = {"num_processors": 4,
                                "output_directory": "/path/to/output",
                                "num_bundles_at_once": 2}

    def test_no_obsolete_keys(self):
        ret = self.validate_config_dict(self.good_config_dict)
        self.assertEquals(0, ret)

    def test_obsolete_keys(self):
        with self.assertRaises(Exception):
            self.validate_config_dict(self.bad_config_dict)


# un-skip this once the CLI is un-classed
@pytest.mark.skip()
class TestCreateSessionDict(unittest.TestCase):

    def setUp(self):
        # setup
        from qap.cli import create_session_dict
        self.create_session_dict = create_session_dict
        self.maxDiff = None

        # inputs
        self.input_subdict = {
            'sub_001': {
                'session_01': {
                    'anatomical_scan': {
                        'anat_1': '/file/site_01/sub_001/session_01/anat_1/mprage.nii.gz'},
                    'functional_scan': {
                        'rest_1': '/file/site_01/sub_001/session_01/rest_1/rest.nii.gz',
                        'rest_2': '/file/site_01/sub_001/session_01/rest_2/rest.nii.gz'},
                    'site_name': 'site_01'}},
            'sub_002': {
                'session_01': {
                    'anatomical_scan': {
                        'anat_1': '/file/site_01/sub_002/session_01/anat_1/mprage.nii.gz'},
                    'functional_scan': {
                        'rest_1': '/file/site_01/sub_002/session_01/rest_1/rest.nii.gz',
                        'rest_2': '/file/site_01/sub_002/session_01/rest_2/rest.nii.gz'},
                    'site_name': 'site_01'}}}

        # outputs
        self.ref_subdict = {
            ("sub_001", "session_01"): {
                "anat_1": {
                    "anatomical_scan": "/file/site_01/sub_001/session_01/anat_1/mprage.nii.gz"},
                "rest_1": {
                    "functional_scan": "/file/site_01/sub_001/session_01/rest_1/rest.nii.gz"},
                "rest_2": {
                    "functional_scan": "/file/site_01/sub_001/session_01/rest_2/rest.nii.gz"},
                "site_name": "site_01"},
            ("sub_002", "session_01"): {
                "anat_1": {
                    "anatomical_scan": "/file/site_01/sub_002/session_01/anat_1/mprage.nii.gz"},
                "rest_1": {
                    "functional_scan": "/file/site_01/sub_002/session_01/rest_1/rest.nii.gz"},
                "rest_2": {
                    "functional_scan": "/file/site_01/sub_002/session_01/rest_2/rest.nii.gz"},
                "site_name": "site_01"}}

    def test_two_subs(self):
        test_subdict = self.create_session_dict(self.input_subdict)
        self.assertDictEqual(self.ref_subdict, test_subdict)


@pytest.mark.quick
class TestCLI(unittest.TestCase):

    def setUp(self):
        # setup
        import os
        from qap import cli

        out_dir = os.path.join(os.getcwd(), "unit_tests_cli", "out")
        work_dir = os.path.join(os.getcwd(), "unit_tests_cli", "work")

        self.cli = cli.QAProtocolCLI(parse_args=False)
        self.cli._config = {}
        self.cli._config["output_directory"] = out_dir
        self.cli._config["working_directory"] = work_dir
        self.cli._config["qap_type"] = "anatomical_spatial"
        self.cli._config["template_head_for_anat"] = "/Library/Python/2.7/site-packages/qap-1.0.8-py2.7.egg/qap/test_data/MNI152_T1_3MM.nii.gz"
        self.cli._num_processors = 1
        self.cli._run_log_dir = "/custom/log/dir"
        self.cli.runargs = {'plugin': 'MultiProc'}
        self.cli._run_name = "qap_unit_test"

        # inputs/outputs
        self.input_subdict = {
            'sub_001': {
                'session_01': {
                    'anatomical_scan': {
                        'anat_1': '/file/site_01/sub_001/session_01/anat_1/mprage.nii.gz'},
                    'functional_scan': {
                        'rest_1': '/file/site_01/sub_001/session_01/rest_1/rest.nii.gz',
                        'rest_2': '/file/site_01/sub_001/session_01/rest_2/rest.nii.gz'},
                    'site_name': 'site_01'}},
            'sub_002': {
                'session_01': {
                    'anatomical_scan': {
                        'anat_1': '/file/site_01/sub_002/session_01/anat_1/mprage.nii.gz'},
                    'functional_scan': {
                        'rest_1': '/file/site_01/sub_002/session_01/rest_1/rest.nii.gz',
                        'rest_2': '/file/site_01/sub_002/session_01/rest_2/rest.nii.gz'},
                    'site_name': 'site_01'}}}

        self.BIDS_input_subdict = {
            "sub-0003001": {
                "ses-1": {
                    "anatomical_scan": {
                        "run-1_T1w": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz'},
                    "creds_path": '',
                    "functional_scan": {
                        "task-rest_run-1": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz',
                        "task-rest_run-2": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz'},
                    "site_name": "site-BMB1"}},
            "sub-0003002": {
                "ses-1": {
                    "anatomical_scan": {
                        "run-1_T1w": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/anat/sub-0003002_ses-1_run-1_T1w.nii.gz'},
                    "creds_path": '',
                    "functional_scan": {
                        "task-rest_run-1": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz',
                        "task-rest_run-2": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz'},
                    "site_name": "site-BMB1"}}}

        self.session_dict = {
            ("sub_001", "session_01"): {
                "anat_1": {
                    "anatomical_scan": "/file/site_01/sub_001/session_01/anat_1/mprage.nii.gz"},
                "rest_1": {
                    "functional_scan": "/file/site_01/sub_001/session_01/rest_1/rest.nii.gz"},
                "rest_2": {
                    "functional_scan": "/file/site_01/sub_001/session_01/rest_2/rest.nii.gz"},
                "site_name": "site_01"},
            ("sub_002", "session_01"): {
                "anat_1": {
                    "anatomical_scan": "/file/site_01/sub_002/session_01/anat_1/mprage.nii.gz"},
                "rest_1": {
                    "functional_scan": "/file/site_01/sub_002/session_01/rest_1/rest.nii.gz"},
                "rest_2": {
                    "functional_scan": "/file/site_01/sub_002/session_01/rest_2/rest.nii.gz"},
                "site_name": "site_01"}}

        self.BIDS_session_dict = {
            ("sub-0003001", "ses-1"): {
                "run-1_T1w": {
                    "anatomical_scan": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz'},
                "creds_path": '',
                "task-rest_run-1": {
                    "functional_scan": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz'},
                "task-rest_run-2": {
                    "functional_scan": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz'},
                "site_name": "site-BMB1"},
            ("sub-0003002", "ses-1"): {
                "run-1_T1w": {
                    "anatomical_scan": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/anat/sub-0003002_ses-1_run-1_T1w.nii.gz'},
                "creds_path": '',
                "task-rest_run-1": {
                    "functional_scan": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz'},
                "task-rest_run-2": {
                    "functional_scan": 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz'},
                "site_name": "site-BMB1"}}

        self.cli._bundles_list = [
            {
                ('sub_001', 'session_01', 'anat_1'): {
                    'anatomical_scan': '/file/path/sub_001/session_01/anatomical_scan/anat_1/file.nii.gz'},
                ('sub_002', 'session_01', 'anat_1'): {
                    'anatomical_scan': '/file/path/sub_002/session_01/anatomical_scan/anat_1/file.nii.gz'},
                ('sub_002', 'session_01', 'anat_2'): {
                    'anatomical_scan': '/file/path/sub_002/session_01/anatomical_scan/anat_2/file.nii.gz'}},
            {
                ('sub_001', 'session_02', 'anat_1'): {
                    'anatomical_scan': '/file/path/sub_001/session_02/anatomical_scan/anat_1/file.nii.gz'}}]

    def tearDown(self):
        import shutil
        try:
            shutil.rmtree(self.cli._config["output_directory"])
        except OSError:
            pass
        try:
            shutil.rmtree(self.cli._config["working_directory"])
        except OSError:
            pass

    def test_create_bundles_one_ses_each(self):
        self.cli._config["num_sessions_at_once"] = 1
        self.cli._sub_dict = self.session_dict
        bundles = self.cli.create_bundles()
        # how many bundles there should be
        #   2 sub-sessions, 1 per bundle = 2 bundles
        self.assertEqual(len(bundles), 2)

    def test_create_bundles_two_ses_each(self):
        self.cli._config["num_sessions_at_once"] = 2
        self.cli._sub_dict = self.session_dict
        bundles = self.cli.create_bundles()
        # how many bundles there should be
        #   2 sub-sessions, 2 per bundle = 1 bundles
        self.assertEqual(len(bundles), 1)

    def test_create_bundles_five_ses_each(self):
        self.cli._config["num_sessions_at_once"] = 5
        self.cli._sub_dict = self.session_dict
        bundles = self.cli.create_bundles()
        # how many bundles there should be
        #   2 sub-sessions, 5 per bundle = 1 bundles
        self.assertEqual(len(bundles), 1)
        # should only have 6 filepaths
        self.assertEqual(len(bundles[0]), 6)

    def test_create_session_dict_nonBIDS(self):
        test_ses_dict = self.cli.create_session_dict(self.input_subdict)
        self.assertDictEqual(self.session_dict, test_ses_dict)

    def test_create_session_dict_BIDS(self):
        test_ses_dict = self.cli.create_session_dict(self.BIDS_input_subdict)
        self.assertDictEqual(self.BIDS_session_dict, test_ses_dict)

    def test_create_bundles_BIDS_session_dict(self):
        self.cli._config["num_sessions_at_once"] = 1
        self.cli._sub_dict = self.BIDS_session_dict
        bundles = self.cli.create_bundles()
        # how many bundles there should be
        #   2 sub-sessions, 1 per bundle = 2 bundles
        self.assertEqual(len(bundles), 2)

    def test_run_one_bundle(self):
        # make sure the number of sessions in the bundle being run matches
        # runs it off self.cli._bundles_list created in setUp
        # inputs 1 for bundle_idx parameter, first bundle in bundles_list, has
        # 3 sessions in it
        wfargs = self.cli.run_one_bundle(1, run=False)
        self.assertEqual(len(wfargs[1]), 3)


@pytest.mark.skip()
@pytest.mark.quick
def init_cli_obj():
    # type: () -> object

    import os
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli", "out")
    work_dir = os.path.join(os.getcwd(), "unit_tests_cli", "work")

    cli_obj = cli.QAProtocolCLI(parse_args=False)
    cli_obj._config = {}
    cli_obj._config["output_directory"] = out_dir
    cli_obj._config["working_directory"] = work_dir
    cli_obj._config["qap_type"] = "anatomical_spatial"
    cli_obj._config["template_head_for_anat"] = "/Library/Python/2.7/site-packages/qap-1.0.8-py2.7.egg/qap/test_data/MNI152_T1_3MM.nii.gz"
    cli_obj._num_processors = 4

    return cli_obj


@pytest.mark.skip()
@pytest.mark.quick
def init_flat_sub_dict_dict():

    ref_flat_dict = {}
    ref_flat_dict[("sub_001","session_01","anat_1")] = {}
    ref_flat_dict[("sub_001","session_01","anat_1")]["anatomical_scan"] = \
        "/file/path/sub_001/session_01/anatomical_scan/anat_1/file.nii.gz"
    ref_flat_dict[("sub_001","session_02","anat_1")] = {}
    ref_flat_dict[("sub_001","session_02","anat_1")]["anatomical_scan"] = \
        "/file/path/sub_001/session_02/anatomical_scan/anat_1/file.nii.gz"
    ref_flat_dict[("sub_002","session_01","anat_1")] = {}
    ref_flat_dict[("sub_002","session_01","anat_1")]["anatomical_scan"] = \
        "/file/path/sub_002/session_01/anatomical_scan/anat_1/file.nii.gz"
    ref_flat_dict[("sub_002","session_01","anat_2")] = {}
    ref_flat_dict[("sub_002","session_01","anat_2")]["anatomical_scan"] = \
        "/file/path/sub_002/session_01/anatomical_scan/anat_2/file.nii.gz"

    return ref_flat_dict


@pytest.mark.skip()
@pytest.mark.quick
def test_submit_cluster_batch_file_for_SGE_s3_paths():

    import os
    import shutil
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli")

    run_name = "run_1"
    num_bundles = 5

    cli_obj = init_cli_obj()

    cli_obj._config["subject_list"] = os.path.join(out_dir, "s3_sublist.yml")
    cli_obj._config["pipeline_config_yaml"] = os.path.join(out_dir,
                                                           "pipeline_config.yml")
    cli_obj._platform = "SGE"

    file_contents, file_path, exec_cmd, confirm_str, cluster_files_dir = \
        cli_obj._prepare_cluster_batch_file(run_name, num_bundles)

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert "SGE" in file_contents.split("\n")[1]
    assert file_contents.split("\n")[3] == "#$ -N run_1"
    assert file_contents.split("\n")[4] == "#$ -t 1-5"
    assert file_contents.split("\n")[6] == "#$ -pe mpi_smp 4"
    assert file_contents.split("\n")[14] == "qap_measures_pipeline.py " \
        "--bundle_idx $SGE_TASK_ID %s %s" \
        % (os.path.join(out_dir, "s3_sublist.yml"),
           os.path.join(out_dir, "pipeline_config.yml"))
    assert exec_cmd == "qsub"
    assert confirm_str == "(?<=Your job-array )\d+"
    assert cluster_files_dir == os.path.join(out_dir, "cluster_files")


@pytest.mark.skip()
@pytest.mark.quick
def test_prepare_cluster_batch_file_for_SGE_with_sublist():

    import os
    import shutil
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli")

    run_name = "run_1"
    num_bundles = 5

    cli_obj = init_cli_obj()

    cli_obj._config["pipeline_config_yaml"] = os.path.join(out_dir,
                                                           "pipeline_config.yml")
    cli_obj._config["subject_list"] = \
        os.path.join(out_dir, "participant_list.yml")
    cli_obj._s3_dict_yml = None
    cli_obj._platform = "SGE"

    file_contents, file_path, exec_cmd, confirm_str, cluster_files_dir = \
        cli_obj._prepare_cluster_batch_file(run_name, num_bundles)

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert "SGE" in file_contents.split("\n")[1]
    assert file_contents.split("\n")[3] == "#$ -N run_1"
    assert file_contents.split("\n")[4] == "#$ -t 1-5"
    assert file_contents.split("\n")[6] == "#$ -pe mpi_smp 4"
    assert file_contents.split("\n")[14] == "qap_anatomical_spatial.py " \
        "--sublist %s --bundle_idx $SGE_TASK_ID %s" \
        % (os.path.join(out_dir, "participant_list.yml"), \
           os.path.join(out_dir, "pipeline_config.yml"))
    assert exec_cmd == "qsub"
    assert confirm_str == "(?<=Your job-array )\d+"
    assert cluster_files_dir == os.path.join(out_dir, "cluster_files")


@pytest.mark.skip()
@pytest.mark.quick
def test_prepare_cluster_batch_file_for_PBS_with_sublist():

    import os
    import shutil
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli")

    run_name = "run_1"
    num_bundles = 5

    cli_obj = init_cli_obj()

    cli_obj._config["pipeline_config_yaml"] = os.path.join(out_dir, \
    	"pipeline_config.yml")
    cli_obj._config["subject_list"] = \
        os.path.join(out_dir, "participant_list.yml")
    cli_obj._s3_dict_yml = None
    cli_obj._platform = "PBS"

    file_contents, file_path, exec_cmd, confirm_str, cluster_files_dir = \
        cli_obj._prepare_cluster_batch_file(run_name, num_bundles)

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert "PBS" in file_contents.split("\n")[1]
    assert file_contents.split("\n")[3] == "#PBS -N run_1"
    assert file_contents.split("\n")[4] == "#PBS -t 1-5"
    assert file_contents.split("\n")[6] == "#PBS -l nodes=1:ppn=4"
    assert file_contents.split("\n")[14] == "qap_anatomical_spatial.py " \
        "--sublist %s --bundle_idx $PBS_ARRAYID %s" \
        % (os.path.join(out_dir, "participant_list.yml"), \
           os.path.join(out_dir, "pipeline_config.yml"))
    assert exec_cmd == "qsub"
    assert confirm_str == "(?<=Your job-array )\d+"
    assert cluster_files_dir == os.path.join(out_dir, "cluster_files")


@pytest.mark.skip()
@pytest.mark.quick
def test_prepare_cluster_batch_file_for_SLURM_with_sublist():

    import os
    import shutil
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli")

    run_name = "run_1"
    num_bundles = 5

    cli_obj = init_cli_obj()

    cli_obj._config["pipeline_config_yaml"] = os.path.join(out_dir, \
    	"pipeline_config.yml")
    cli_obj._config["subject_list"] = \
        os.path.join(out_dir, "participant_list.yml")
    cli_obj._s3_dict_yml = None
    cli_obj._platform = "SLURM"

    file_contents, file_path, exec_cmd, confirm_str, cluster_files_dir = \
        cli_obj._prepare_cluster_batch_file(run_name, num_bundles)

    try:
        shutil.rmtree(out_dir)
    except:
        pass

    assert "SLURM" in file_contents.split("\n")[1]
    assert file_contents.split("\n")[2] == "#SBATCH --job-name=run_1"
    assert file_contents.split("\n")[3] == "#SBATCH --array=1-5"
    assert file_contents.split("\n")[4] == "#SBATCH --cpus-per-task=4"
    assert file_contents.split("\n")[13] == "qap_anatomical_spatial.py " \
        "--sublist %s --bundle_idx $SLURM_ARRAY_TASK_ID %s" \
        % (os.path.join(out_dir, "participant_list.yml"), \
           os.path.join(out_dir, "pipeline_config.yml"))
    assert exec_cmd == "sbatch"
    assert confirm_str == "(?<=Submitted batch job )\d+"
    assert cluster_files_dir == os.path.join(out_dir, "cluster_files")


@pytest.mark.skip()
@pytest.mark.quick
def test_create_bundles_one_sub_per_bundle():

    from qap import cli

    cli_obj = init_cli_obj()
    cli_obj._num_participants_at_once = 1

    flat_sub_dict_dict = init_flat_sub_dict_dict()

    bundles = cli_obj._create_bundles(flat_sub_dict_dict)

    # how many bundles there should be
    #   4 sub-session-scans, 1 per bundle = 4 bundles
    assert len(bundles) == 4

    # how many sub-session-scans should be in each bundle
    assert len(bundles[0]) == cli_obj._num_participants_at_once


@pytest.mark.skip()
@pytest.mark.quick
def test_create_bundles_two_subs_per_bundle():

    cli_obj = init_cli_obj()
    cli_obj._num_participants_at_once = 2

    flat_sub_dict_dict = init_flat_sub_dict_dict()

    bundles = cli_obj._create_bundles(flat_sub_dict_dict)

    # how many bundles there should be
    #   4 sub-session-scans, 2 per bundle = 2 bundles
    assert len(bundles) == 2

    # how many sub-session-scans should be in each bundle
    assert len(bundles[0]) == cli_obj._num_participants_at_once


@pytest.mark.skip()
@pytest.mark.quick
def test_create_bundles_six_subs_per_bundle():

    cli_obj = init_cli_obj()
    cli_obj._num_participants_at_once = 6

    flat_sub_dict_dict = init_flat_sub_dict_dict()

    bundles = cli_obj._create_bundles(flat_sub_dict_dict)

    # how many bundles there should be
    #   4 sub-session-scans, 6 per bundle = 1 bundle
    assert len(bundles) == 1

    # how many sub-session-scans should be in each bundle
    #   number of sub-session-scans less than "num subs per bundle", so should
    #   just be the total number of sub-session-scans
    assert len(bundles[0]) == 4


@pytest.mark.skip()
@pytest.mark.quick
def test_get_num_bundles():

    cli_obj = init_cli_obj()
    cli_obj._config["num_participants_at_once"] = 2
    size_of_sublist = 5

    num_bundles = cli_obj._get_num_bundles(size_of_sublist)
    assert num_bundles == 3

    cli_obj._config["num_participants_at_once"] = 7
    size_of_sublist = 5

    num_bundles = cli_obj._get_num_bundles(size_of_sublist)
    assert num_bundles == 1


@pytest.mark.skip()
@pytest.mark.quick
def test_run_workflow_build_only():

    import os
    from qap import cli
    import pkg_resources as p

    anat_scan_1 = p.resource_filename("qap", os.path.join("test_data", "input_data", "site_1",
                                                        "sub_001", "session_01", "anatomical_scan",
                                                        "anat_1", "anat.nii.gz"))
    anat_scan_2 = p.resource_filename("qap", os.path.join("test_data", "input_data", "site_1",
                                                        "sub_001", "session_02", "anatomical_scan",
                                                        "anat_1", "anat.nii.gz"))

    one_bundle = {('sub_001', 'session_01', 'anat_1'): {'anatomical_scan': anat_scan_1},
                  ('sub_001', 'session_02', 'anat_1'): {'anatomical_scan': anat_scan_2}}

    run_name = "test_run_workflow"
    run_args = {'plugin': 'MultiProc'}
    bundle_idx = 1

    cli_obj = init_cli_obj()

    args = (one_bundle, one_bundle.keys(), cli_obj._config, run_name, run_args, bundle_idx)

    workflow = cli._run_workflow(args, False)

    node_names = workflow.list_node_names()

    terminal_nodes = [
        'gather_header_info_sub_001_session_01_anat_1',
        'gather_header_info_sub_001_session_02_anat_1',
        'qap_anatomical_spatial_sub_001_session_01_anat_1',
        'qap_anatomical_spatial_sub_001_session_02_anat_1',
        'qap_anatomical_spatial_to_json_sub_001_session_01_anat_1',
        'qap_anatomical_spatial_to_json_sub_001_session_02_anat_1',
        'qap_header_to_json_sub_001_session_01_anat_1',
        'qap_header_to_json_sub_001_session_02_anat_1']

    for node in terminal_nodes:
        assert node in node_names
