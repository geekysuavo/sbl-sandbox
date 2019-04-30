
# required imports.
from doit.tools import run_once

# task_eigen3: task generator used to download the eigen3 sources.
def task_eigen3():
  eigen_url = 'http://bitbucket.org/eigen/eigen/get/3.3.7.tar.gz'
  return {
    'actions': [f'wget -q {eigen_url} -O eigen3.tgz',
                'tar xf eigen3.tgz',
                'mv eigen-* eigen3'],
    'uptodate': [run_once],
    'targets': ['eigen3']
  }

