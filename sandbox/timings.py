
# required modules.
from doit.tools import create_folder
from sandbox.util import solvers
import subprocess
import timeit
import pickle
import gzip
import os

# task_timing: task generator for measuring runtimes.
def task_timing():
  # timing: task action function for measuring a single runtime.
  def timing(meth, bins, targets):
    # run all required timings.
    results = [None] * len(bins)
    for idx in range(len(bins)):
      # get the number of replicates and size info.
      nreps = 5
      T = [0] * nreps
      (m, n, binary) = bins[idx]

      # build the execution statement.
      prep = 'import subprocess'
      stmt = (f'proc = subprocess.run("{binary}",' +
              f' stdout=subprocess.DEVNULL,' +
              f' stderr=subprocess.DEVNULL)')

      # measure timings.
      for rep in range(nreps):
        T[rep] = timeit.timeit(stmt=stmt, setup=prep, number=100) / 100

      # store the timings.
      results[idx] = {'method': meth, 'm': m, 'n': n, 'times': tuple(T)}

    # write the results to the target file.
    with gzip.open(targets[0], 'wb') as f:
      pickle.dump(results, f)

  # ---

  # loop over each method.
  methods = tuple({s[0] for s in solvers()})
  for meth in methods:
    bins = [(m, n, os.path.join('bin', f'{meth}-{m}-{n}'))
            for M, m, n in solvers() if M == meth]
    target_dir = os.path.join('expts', 'timing')
    target = os.path.join(target_dir, f'{meth}.gz')

    yield {
      'name': meth,
      'actions': [(create_folder, [target_dir]),
                  (timing, [meth, bins])],
      'file_dep': [b[-1] for b in bins],
      'targets': [target]
    }

