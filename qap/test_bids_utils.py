
import pytest

BIDS_s3_list = [
 'data/RawDataBIDS/sub-A0200/ses-A01/anat/sub-A0200_ses-A01_T1w.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/anat/sub-A0200_ses-A01_T1w_brainmask.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/anat/sub-A0200_ses-A01_T1w_skullstripped.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/fmap/sub-A0200_ses-A01_magnitude1.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/fmap/sub-A0200_ses-A01_phasediff.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY1_bold.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY1_events.tsv',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY1_physio.json',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY1_physio.tsv.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY2_bold.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY2_physio.json',
 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY2_physio.tsv.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/anat/sub-A0200_ses-A02_T1w.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/anat/sub-A0200_ses-A02_T1w_brainmask.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/anat/sub-A0200_ses-A02_T1w_skullstripped.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/fmap/sub-A0200_ses-A02_magnitude1.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/fmap/sub-A0200_ses-A02_phasediff.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY1_bold.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY1_events.tsv',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY1_physio.json',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY1_physio.tsv.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY2_bold.nii.gz',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY2_physio.json',
 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY2_physio.tsv.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/anat/sub-A0300_ses-A01_T1w.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/anat/sub-A0300_ses-A01_T1w_brainmask.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/anat/sub-A0300_ses-A01_T1w_skullstripped.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/fmap/sub-A0300_ses-A01_magnitude1.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/fmap/sub-A0300_ses-A01_phasediff.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY1_bold.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY1_events.tsv',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY1_physio.json',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY1_physio.tsv.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY2_bold.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY2_physio.json',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY2_physio.tsv.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY3_bold.nii.gz',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY3_physio.json',
 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY3_physio.tsv.gz'
]


BIDS_CoRR = [
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/T1w.json',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/sub-0003001_sessions.tsv',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/anat/sub-0003002_ses-1_run-1_T1w.nii.gz',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz',
     's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/sub-0003002_sessions.tsv'
]


def test_extract_bids_data_S3():

    from qap.bids_utils import extract_bids_data

    ref_subdict = {
        'sub-A0200': {
            'ses-A01': {
                'anatomical_scan': {
                    'T1w': 'data/RawDataBIDS/sub-A0200/ses-A01/anat/sub-A0200_ses-A01_T1w.nii.gz'},
                'functional_scan': {
                    'bold_task-STUDY1': 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY1_bold.nii.gz',
                    'bold_task-STUDY2': 'data/RawDataBIDS/sub-A0200/ses-A01/func/sub-A0200_ses-A01_task-STUDY2_bold.nii.gz'}},
            'ses-A02': {
                'anatomical_scan': {
                    'T1w': 'data/RawDataBIDS/sub-A0200/ses-A02/anat/sub-A0200_ses-A02_T1w.nii.gz'},
                'functional_scan': {
                    'bold_task-STUDY1': 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY1_bold.nii.gz',
                    'bold_task-STUDY2': 'data/RawDataBIDS/sub-A0200/ses-A02/func/sub-A0200_ses-A02_task-STUDY2_bold.nii.gz'}}},
        'sub-A0300': {
            'ses-A01': {
                'anatomical_scan': {
                    'T1w': 'data/RawDataBIDS/sub-A0300/ses-A01/anat/sub-A0300_ses-A01_T1w.nii.gz'},
                'functional_scan': {
                    'bold_task-STUDY1': 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY1_bold.nii.gz',
                    'bold_task-STUDY2': 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY2_bold.nii.gz',
                    'bold_task-STUDY3': 'data/RawDataBIDS/sub-A0300/ses-A01/func/sub-A0300_ses-A01_task-STUDY3_bold.nii.gz'}}}}

    subdict = extract_bids_data(BIDS_s3_list)

    corr_ref_subdict = {
        'sub-0003001': {
            'ses-1': {
                'anatomical_scan': {
                    'T1w': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/anat/sub-0003001_ses-1_run-1_T1w.nii.gz'},
                'functional_scan': {
                    'bold_task-rest_run-1': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-1_bold.nii.gz',
                    'bold_task-rest_run-2': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003001/ses-1/func/sub-0003001_ses-1_task-rest_run-2_bold.nii.gz'}}},
        'sub-0003002': {
            'ses-1': {
                'anatomical_scan': {
                    'T1w': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/anat/sub-0003002_ses-1_run-1_T1w.nii.gz'},
                'functional_scan': {
                    'bold_task-rest_run-1': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-1_bold.nii.gz',
                    'bold_task-rest_run-2': 's3://fcp-indi/data/Projects/CORR/RawDataBIDS/BMB_1/sub-0003002/ses-1/func/sub-0003002_ses-1_task-rest_run-2_bold.nii.gz'}}}}

    corr_subdict = extract_bids_data(BIDS_CoRR)

    assert ref_subdict == subdict
    assert corr_ref_subdict == corr_subdict

