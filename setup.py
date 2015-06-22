
def main():

    from glob import glob
    from setuptools import setup

    setup(name='qap',
          version='0.1.0',
          description='A collection of three quality assessment pipelines ' \
                      'for anatomical MRI and functional MRI scans.',
          url='https://github.com/preprocessed-connectomes-project/' \
              'quality-assessment-protocol',
          author='Cameron Craddock, Zarrar Shehzad, Steven Giavassis,' \
                 'Daniel Clark',
          author_email='cameron.craddock@childmind.org',
          license='',
          packages=['qap'],
          package_data={'qap': ['inpoint*.txt',
                                'test_data/*.nii.gz',
                                'test_data/*/*/*/*/*']},
          scripts=glob("scripts/*"),
          zip_safe=False)
          
          
          
if __name__ == "__main__":

    import os
    import sys
    
    local_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    os.chdir(local_path)
    sys.path.insert(0, local_path)
    
    main()
