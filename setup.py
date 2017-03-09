
def main():

    from glob import glob
    from setuptools import setup

    exec(open('qap/version.py').read())
    setup(name='qap',
          version=__version__,
          description='A collection of three quality assessment pipelines '
                      'for anatomical MRI and functional MRI scans.',
          author='Cameron Craddock, Zarrar Shehzad, Steven Giavasis,'
                 'Daniel Clark',
          author_email='cameron.craddock@childmind.org',
          url='https://github.com/preprocessed-connectomes-project/'
              'quality-assessment-protocol',
          download_url='https://github.com/preprocessed-connectomes-project/'
                       'quality-assessment-protocol/archive/v1.0.8b.zip',
          license='',
          packages=['qap', 'qap.viz'],
          package_data={'qap': ['inpoint*.txt',
                                'viz/html/*.html',
                                'test_data/input_data/site_1/sub_001/session_01/anatomical_scan/anat_1/*.nii.gz',
                                'test_data/input_data/site_1/sub_001/session_02/anatomical_scan/anat_1/*.nii.gz',
                                'test_data/*.nii.gz',
                                'test_data/*.dot',
                                'test_data/*.1D',
                                'test_data/*.p',
                                'test_data/*.csv']},
          scripts=glob("scripts/*"),
          requires=['INDI_Tools (>=0.0.6)', 'Pillow (>=3.2.0)', 
                    'PyPDF2 (>=1.26.0)', 'PyYAML (>=3.11)', 
                    'argparse (>=1.2.1)', 'boto3 (>=1.3.1)', 
                    'botocore (>=1.4.22)', 'cycler (>=0.10.0)', 
                    'decorator (>=4.0.9)', 'docutils (>=0.12)', 
                    'future (>=0.15.2)', 'futures (>=3.0.5)', 
                    'html5lib (>=1.0b8)', 'httplib2 (>=0.9.2)', 
                    'jmespath (>=0.9.0)', 'matplotlib (>=1.5.1)', 
                    'networkx (>=1.11)', 'nibabel (>=2.0.2)', 
                    'nitime (>=0.6)', 'nipype (>=0.12.1)', 
                    'nose (>=1.3.7)', 'numpy (>=1.11.0)', 'pandas (>=0.18.1)',
                    'prov (>=1.4.0)', 'pyparsing (>=2.1.4)', 
                    'python_dateutil (>=2.5.3)', 'pytz (>=2016.4)', 
                    'reportlab (>=3.3.0)', 'scipy (>=0.17.1)', 
                    'seaborn (>=0.7.0)', 'simplejson (>=3.8.2)', 
                    'six (>=1.10.0)', 'traits (>=4.5.0)', 'wsgiref (>=0.1.2)',
                    'xhtml2pdf (>=0.1a4)', 'configparser (>=3.5.0)',
                    'lockfile (>=0.12)'],
          install_requires=['numpy >=1.11.0', 'INDI-Tools >=0.0.6', 
                    'PyPDF2 >=1.26.0', 'PyYAML >=3.11', 'Pillow >=3.2.0', 
                    'argparse >=1.2.1', 'boto3 >=1.3.1', 
                    'botocore >=1.4.22', 'cycler >=0.10.0', 
                    'decorator >=4.0.9', 'docutils >=0.12', 
                    'future >=0.15.2', 'futures >=3.0.5', 
                    'html5lib >=1.0b8', 'httplib2 >=0.9.2', 
                    'jmespath >=0.9.0', 'matplotlib >=1.5.1', 
                    'networkx >=1.11', 'nibabel >=2.0.2', 
                    'nitime >=0.6', 'nipype >=0.12.1', 
                    'nose >=1.3.7', 'pandas >=0.18.1',
                    'prov >=1.4.0', 'pyparsing >=2.1.4', 
                    'python-dateutil >=2.5.3', 'pytz >=2016.4', 
                    'reportlab >=3.3.0', 'scipy >=0.17.1', 
                    'seaborn >=0.7.0', 'simplejson >=3.8.2', 
                    'six >=1.10.0', 'traits >=4.5.0', 'wsgiref >=0.1.2',
                    'xhtml2pdf >=0.1a4', 'configparser >=3.5.0', 
                    'numpy >=1.11.0', 'lockfile >=0.12'],
          zip_safe=False)

if __name__ == "__main__":

    import os
    import sys

    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)

    main()
