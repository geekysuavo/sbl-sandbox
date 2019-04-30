
# required modules.
from doit.tools import create_folder
from sandbox.util import expts, solve
import subprocess
import pickle
import gzip
import os

# task_experiment: task generator for executing all experiments.
def task_experiment():
  # experiment: task action function for building per-task pickles.
  def experiment(ident, parms, targets):
    # prepare a dictionary to hold the results.
    methods = parms['methods']
    results = {meth: None for meth in methods}

    # execute the solver for each method.
    for meth in methods:
      (out, err) = solve(parms, meth)
      results[meth] = (out, err)

    # write the results to the target file.
    parms['results'] = results
    with gzip.open(targets[0], 'wb') as f:
      pickle.dump(parms, f)

  # ---

  # loop over each task of each experiment.
  for expt, tasks in expts.items():
    expt_dir = os.path.join('expts', expt)
    for idx in range(len(tasks)):
      ident = f'{idx+1:06d}'
      name = f'{expt}-{ident}'
      target = os.path.join('expts', expt, f'{ident}.gz')

      parms = tasks[idx]
      (m, n) = (parms['m'], parms['n'])
      binaries = [os.path.join('bin', f'{meth}-{m}-{n}')
                  for meth in parms['methods']]

      yield {
        'name': name,
        'actions': [(create_folder, [expt_dir]),
                    (experiment, [ident, parms])],
        'file_dep': binaries,
        'targets': [target]
      }

