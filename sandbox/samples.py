
# required imports.
import os

# task_nmr1d: task generator used to run the nmr1d sample experiment.
def task_nmr1d():
  sample_dir = os.path.join('samples', 'nmr1d')
  bin_ista = os.path.join(sample_dir, 'ista')
  bin_fmf = os.path.join(sample_dir, 'fmf')

  def buildstr(B):
    return f'g++ -std=c++14 -O3 {B}.cc -o {B} -lfftw3 -lm'

  def runstr(B):
    return f'{B} > {B}.out'

  yield {
    'name': 'ista',
    'actions': [buildstr(bin_ista),
                runstr(bin_ista)],
    'file_dep': [f'{bin_ista}.cc'],
    'targets': [bin_ista, f'{bin_ista}.out']
  }

  yield {
    'name': 'fmf',
    'actions': [buildstr(bin_fmf),
                runstr(bin_fmf)],
    'file_dep': [f'{bin_fmf}.cc'],
    'targets': [bin_fmf, f'{bin_fmf}.out']
  }

# task_nmr2d: task generator used to run the nmr2d sample experiment.
def task_nmr2d():
  sample_dir = os.path.join('samples', 'nmr2d')
  bin_ista = os.path.join(sample_dir, 'ista')
  bin_fmf = os.path.join(sample_dir, 'fmf')

  def buildstr(B):
    return f'g++ -std=c++14 -O3 {B}.cc -o {B} -lfftw3_threads -lfftw3 -lm'

  def runstr(B):
    return f'{B} > {B}.out'

  yield {
    'name': 'ista',
    'actions': [buildstr(bin_ista),
                runstr(bin_ista)],
    'file_dep': [f'{bin_ista}.cc'],
    'targets': [bin_ista, f'{bin_ista}.out']
  }

  yield {
    'name': 'fmf',
    'actions': [buildstr(bin_fmf),
                runstr(bin_fmf)],
    'file_dep': [f'{bin_fmf}.cc'],
    'targets': [bin_fmf, f'{bin_fmf}.out']
  }

