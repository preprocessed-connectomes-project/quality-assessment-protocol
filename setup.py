import os
from setuptools import setup

exec(open('qap/version.py').read())

def list_files(d):
    files = []
    for root, _, filenames in os.walk('qap/' + d):
        root = root[len('qap/'):]
        files += map(lambda f: os.path.join(root, f), filenames)
    return files

package_data = { 'qap': list_files("viz/html") + list_files("configs") }

setup(name='qap',
      version=__version__,
      description='A collection of three quality assessment pipelines '
                  'for anatomical MRI and functional MRI scans.',
      author='Cameron Craddock, Zarrar Shehzad, Steven Giavasis, Daniel Clark',
      author_email='cameron.craddock@childmind.org',
      url='https://github.com/preprocessed-connectomes-project/'
          'quality-assessment-protocol',
      download_url='https://github.com/preprocessed-connectomes-project/'
                  'quality-assessment-protocol/archive/v1.0.8b.zip',
      license='',
      packages=['qap', 'qap.viz', 'qap.workflows'],
      package_data=package_data,
      scripts=list_files("scripts"),
      entry_points={
          'console_scripts': [
              'qap = qap.__main__:main'
          ]
      },
      install_requires=[
          'INDI_Tools>=0.0.6',
          'PyPDF2>=1.26.0',
          'PyYAML>=3.11',
          'argparse>=1.2.1',
          'boto3>=1.3.1',
          'docutils>=0.12',
          'matplotlib>=1.5.1',
          'nibabel>=2.0.2',
          'nitime>=0.6',
          'nipype==0.13.1',
          'nose>=1.3.7',
          'numpy>=1.11.0',
          'pandas>=0.18.1',
          'scipy>=0.17.1',
          'seaborn>=0.7.0',
          'six>=1.10.0',
          'traits>=4.5.0',
          'xhtml2pdf>=0.1a4'])