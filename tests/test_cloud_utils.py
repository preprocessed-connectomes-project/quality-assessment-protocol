
import os
import unittest
import pytest
import random
import shutil
import string
import tempfile

import boto
import boto3
from moto import mock_s3

import qap.cloud_utils as cloud_utils


@pytest.mark.long()
class TestCloudUtils(unittest.TestCase):

    mock_s3 = mock_s3()

    def setUp(self):
        self.mock_s3.start()
        
        self.local_dir = os.path.dirname(os.path.abspath(__file__))
        self.working_dir = os.path.join(self.local_dir, 'workflow_output')

        self.bucket_name = 'open-dataset'

        conn = boto3.resource('s3')
        self.bucket = conn.create_bucket(Bucket=self.bucket_name)

    def tearDown(self):
        self.mock_s3.stop()

    def test_download_single_s3_path(self):

        file_name = 'anat.nii.gz'
        file_content = 'test_data'
        s3 = boto3.client('s3')
        s3.put_object(Bucket=self.bucket_name, Key=file_name, Body=file_content)

        cfg_dict = {"working_directory": self.working_dir}
        cloud_utils.download_single_s3_path("s3://{}/{}".format(self.bucket_name, file_name), cfg_dict)

        with open(os.path.join(self.working_dir, file_name), 'r') as f:
            content = f.read()

            assert content == file_content

    def test_copy_directory_to_s3(self):

        test_dir = os.path.join(self.working_dir, 'test_copy_directory_to_s3')

        shutil.rmtree(test_dir, ignore_errors=True)
        os.mkdir(test_dir)

        def file_name_gen(N):
            return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

        files = [
            os.path.join(*[file_name_gen(6) for _ in range(random.choice(range(1, 5)))])
            for _ in range(10)
        ]

        for f in files:
            fpath = os.path.join(test_dir, f) 
            if not os.path.exists(os.path.dirname(fpath)):
                os.makedirs(os.path.dirname(fpath))
            with open(fpath, 'w') as fd:
                fd.write(f[::-1])  # reversed filename string

        credential_fd, credential = tempfile.mkstemp()
        with open(credential, 'w') as f:
            f.write("AWSAccessKeyId=ASDASDASD\n")
            f.write("AWSSecretKey=ASDASDASD\n")
        os.close(credential_fd)

        cloud_utils.copy_directory_to_s3(
            test_dir,
            's3://{}/test_copy_directory_to_s3'.format(self.bucket_name),
            {'s3_write_credentials': credential})

        for f in files:
            bobj = self.bucket.Object(key='test_copy_directory_to_s3' + '/' + f)
            bobj.get()
