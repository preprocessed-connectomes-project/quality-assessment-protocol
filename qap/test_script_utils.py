
import pytest

# site_1
# 0026005 - 3 sessions, 1 anat and 2 func each
# 0026006 - 3 sessions, first two have 1 anat and 2 func each, third has
#                       1 anat and no func
# site_2
# 0033005 - 2 sessions, 1 anat and 2 func each
nonBIDS_s3_list = [
 'data/RawData/site_1/0026005/session_1/anat_1/anat.nii.gz',
 'data/RawData/site_1/0026005/session_1/dti_1/dti.bval',
 'data/RawData/site_1/0026005/session_1/dti_1/dti.bvec',
 'data/RawData/site_1/0026005/session_1/dti_1/dti.nii.gz',
 'data/RawData/site_1/0026005/session_1/rest_1/rest.nii.gz',
 'data/RawData/site_1/0026005/session_1/rest_2/rest.nii.gz',
 'data/RawData/site_1/0026005/session_2/anat_1/anat.nii.gz',
 'data/RawData/site_1/0026005/session_2/dti_1/dti.bval',
 'data/RawData/site_1/0026005/session_2/dti_1/dti.bvec',
 'data/RawData/site_1/0026005/session_2/dti_1/dti.nii.gz',
 'data/RawData/site_1/0026005/session_2/rest_1/rest.nii.gz',
 'data/RawData/site_1/0026005/session_2/rest_2/rest.nii.gz',
 'data/RawData/site_1/0026005/session_3/anat_1/anat.nii.gz',
 'data/RawData/site_1/0026005/session_3/dti_1/dti.bval',
 'data/RawData/site_1/0026005/session_3/dti_1/dti.bvec',
 'data/RawData/site_1/0026005/session_3/dti_1/dti.nii.gz',
 'data/RawData/site_1/0026005/session_3/rest_1/rest.nii.gz',
 'data/RawData/site_1/0026005/session_3/rest_2/rest.nii.gz',
 'data/RawData/site_1/0026006/session_1/anat_1/anat.nii.gz',
 'data/RawData/site_1/0026006/session_1/dti_1/dti.bval',
 'data/RawData/site_1/0026006/session_1/dti_1/dti.bvec',
 'data/RawData/site_1/0026006/session_1/dti_1/dti.nii.gz',
 'data/RawData/site_1/0026006/session_1/rest_1/rest.nii.gz',
 'data/RawData/site_1/0026006/session_1/rest_2/rest.nii.gz',
 'data/RawData/site_1/0026006/session_2/anat_1/anat.nii.gz',
 'data/RawData/site_1/0026006/session_2/dti_1/dti.bval',
 'data/RawData/site_1/0026006/session_2/dti_1/dti.bvec',
 'data/RawData/site_1/0026006/session_2/dti_1/dti.nii.gz',
 'data/RawData/site_1/0026006/session_2/rest_1/rest.nii.gz',
 'data/RawData/site_1/0026006/session_2/rest_2/rest.nii.gz',
 'data/RawData/site_1/0026006/session_3/anat_1/anat.nii.gz',
 'data/RawData/site_2/0033005/session_1/anat_1/anat.nii.gz',
 'data/RawData/site_2/0033005/session_1/dti_1/dti.bval',
 'data/RawData/site_2/0033005/session_1/dti_1/dti.bvec',
 'data/RawData/site_2/0033005/session_1/dti_1/dti.nii.gz',
 'data/RawData/site_2/0033005/session_1/rest_1/rest.nii.gz',
 'data/RawData/site_2/0033005/session_1/rest_2/rest.nii.gz',
 'data/RawData/site_2/0033005/session_2/anat_1/anat.nii.gz',
 'data/RawData/site_2/0033005/session_2/dti_1/dti.bval',
 'data/RawData/site_2/0033005/session_2/dti_1/dti.bvec',
 'data/RawData/site_2/0033005/session_2/dti_1/dti.nii.gz',
 'data/RawData/site_2/0033005/session_2/rest_1/rest.nii.gz',
 'data/RawData/site_2/0033005/session_2/rest_2/rest.nii.gz']


@pytest.mark.quick
def test_parse_raw_data_list():

    from qap.script_utils import parse_raw_data_list

    # we are starting in the directory containing the site folders!
    site_folder = "/home/data"

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

    sub_dict = parse_raw_data_list(filepath_list, site_folder)

    assert ref_sub_dict == sub_dict


@pytest.mark.quick
def test_parse_raw_data_list_S3():

    from qap.script_utils import parse_raw_data_list

    ref_subdict = {'0026005': {'session_1': {'anatomical_scan': {'anat_1': 'data/RawData/site_1/0026005/session_1/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_1/0026005/session_1/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_1/0026005/session_1/rest_2/rest.nii.gz'},
   'site_name': 'site_1'},
  'session_2': {'anatomical_scan': {'anat_1': 'data/RawData/site_1/0026005/session_2/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_1/0026005/session_2/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_1/0026005/session_2/rest_2/rest.nii.gz'},
   'site_name': 'site_1'},
  'session_3': {'anatomical_scan': {'anat_1': 'data/RawData/site_1/0026005/session_3/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_1/0026005/session_3/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_1/0026005/session_3/rest_2/rest.nii.gz'},
   'site_name': 'site_1'}},
 '0026006': {'session_1': {'anatomical_scan': {'anat_1': 'data/RawData/site_1/0026006/session_1/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_1/0026006/session_1/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_1/0026006/session_1/rest_2/rest.nii.gz'},
   'site_name': 'site_1'},
  'session_2': {'anatomical_scan': {'anat_1': 'data/RawData/site_1/0026006/session_2/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_1/0026006/session_2/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_1/0026006/session_2/rest_2/rest.nii.gz'},
   'site_name': 'site_1'},
  'session_3': {'anatomical_scan': {'anat_1': 'data/RawData/site_1/0026006/session_3/anat_1/anat.nii.gz'},
   'site_name': 'site_1'}},
 '0033005': {'session_1': {'anatomical_scan': {'anat_1': 'data/RawData/site_2/0033005/session_1/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_2/0033005/session_1/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_2/0033005/session_1/rest_2/rest.nii.gz'},
   'site_name': 'site_2'},
  'session_2': {'anatomical_scan': {'anat_1': 'data/RawData/site_2/0033005/session_2/anat_1/anat.nii.gz'},
   'functional_scan': {'rest_1': 'data/RawData/site_2/0033005/session_2/rest_1/rest.nii.gz',
    'rest_2': 'data/RawData/site_2/0033005/session_2/rest_2/rest.nii.gz'},
   'site_name': 'site_2'}}}

    subdict = parse_raw_data_list(nonBIDS_s3_list, "data/RawData",
                                  s3_bucket=True)

    assert ref_subdict == subdict


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
