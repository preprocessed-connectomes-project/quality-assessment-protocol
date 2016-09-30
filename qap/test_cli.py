
import pytest


@pytest.mark.quick
def test_cli():

    err = "some minor refactoring needed! no proper error message when you send in an S3 list as a subject list. flatten sublist, run_here, etc. seems convoluted"
    pass


@pytest.mark.quick
def init_cli_obj():

    import os
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli")

    cli_obj = cli.QAProtocolCLI(parse_args=False)
    cli_obj._config = {}
    cli_obj._config["output_directory"] = out_dir
    cli_obj._config["qap_type"] = "anatomical_spatial"
    cli_obj._num_subjects_per_bundle = 4

    return cli_obj


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


@pytest.mark.quick
def test_prepare_cluster_batch_file_for_SGE_with_s3_dict():

    import os
    import shutil
    from qap import cli

    out_dir = os.path.join(os.getcwd(), "unit_tests_cli")

    run_name = "run_1"
    num_bundles = 5

    cli_obj = init_cli_obj()

    cli_obj._s3_dict_yml = os.path.join(out_dir, "s3_dict_yml.yml")
    cli_obj._config["pipeline_config_yaml"] = os.path.join(out_dir, \
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
    assert file_contents.split("\n")[14] == "qap_anatomical_spatial.py " \
        "--s3_dict_yml %s --bundle_idx $SGE_TASK_ID %s" \
        % (os.path.join(out_dir, "s3_dict_yml.yml"), \
           os.path.join(out_dir, "pipeline_config.yml"))
    assert exec_cmd == "qsub"
    assert confirm_str == "(?<=Your job-array )\d+"
    assert cluster_files_dir == os.path.join(out_dir, "cluster_files")


@pytest.mark.quick
def test_prepare_cluster_batch_file_for_SGE_with_sublist():

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


@pytest.mark.quick
def test_create_flat_sub_dict_dict():

    from qap import cli

    cli_obj = init_cli_obj()

    subdict = {}
    subdict["sub_001"] = {}
    subdict["sub_001"]["session_01"] = {}
    subdict["sub_001"]["session_01"]["anatomical_scan"] = {}
    subdict["sub_001"]["session_01"]["anatomical_scan"]["anat_1"] = \
	    "/file/path/sub_001/session_01/anatomical_scan/anat_1/file.nii.gz"

    subdict["sub_001"]["session_02"] = {}
    subdict["sub_001"]["session_02"]["anatomical_scan"] = {}
    subdict["sub_001"]["session_02"]["anatomical_scan"]["anat_1"] = \
	    "/file/path/sub_001/session_02/anatomical_scan/anat_1/file.nii.gz"

    subdict["sub_002"] = {}
    subdict["sub_002"]["session_01"] = {}
    subdict["sub_002"]["session_01"]["anatomical_scan"] = {}
    subdict["sub_002"]["session_01"]["anatomical_scan"]["anat_1"] = \
	    "/file/path/sub_002/session_01/anatomical_scan/anat_1/file.nii.gz"

    subdict["sub_002"]["session_01"]["anatomical_scan"]["anat_2"] = \
	    "/file/path/sub_002/session_01/anatomical_scan/anat_2/file.nii.gz"

    flat_sub_dict_dict = cli_obj.create_flat_sub_dict_dict(subdict)

    ref_flat_dict = init_flat_sub_dict_dict()

    assert flat_sub_dict_dict == ref_flat_dict


@pytest.mark.quick
def test_create_bundles_one_sub_per_bundle():

    from qap import cli

    cli_obj = init_cli_obj()
    cli_obj._num_subjects_per_bundle = 1

    flat_sub_dict_dict = init_flat_sub_dict_dict()

    bundles = cli_obj._create_bundles(flat_sub_dict_dict)

    # how many bundles there should be
    #   4 sub-session-scans, 1 per bundle = 4 bundles
    assert len(bundles) == 4

    # how many sub-session-scans should be in each bundle
    assert len(bundles[0]) == cli_obj._num_subjects_per_bundle


@pytest.mark.quick
def test_create_bundles_two_subs_per_bundle():

    from qap import cli

    cli_obj = init_cli_obj()
    cli_obj._num_subjects_per_bundle = 2

    flat_sub_dict_dict = init_flat_sub_dict_dict()

    bundles = cli_obj._create_bundles(flat_sub_dict_dict)

    # how many bundles there should be
    #   4 sub-session-scans, 2 per bundle = 2 bundles
    assert len(bundles) == 2

    # how many sub-session-scans should be in each bundle
    assert len(bundles[0]) == cli_obj._num_subjects_per_bundle


@pytest.mark.quick
def test_create_bundles_six_subs_per_bundle():

    from qap import cli

    cli_obj = init_cli_obj()
    cli_obj._num_subjects_per_bundle = 6

    flat_sub_dict_dict = init_flat_sub_dict_dict()

    bundles = cli_obj._create_bundles(flat_sub_dict_dict)

    # how many bundles there should be
    #   4 sub-session-scans, 6 per bundle = 1 bundle
    assert len(bundles) == 1

    # how many sub-session-scans should be in each bundle
    #   number of sub-session-scans less than "num subs per bundle", so should
    #   just be the total number of sub-session-scans
    assert len(bundles[0]) == 4
