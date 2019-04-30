
# required modules.
from sandbox.util import expts
import pickle
import gzip
import os

# task_package: task generator for combining all experimental results.
def task_package():
  # package: task action function for repacking per-task data
  # into a larger pickle that is easier to transport.
  def package(expt, targets):
    # build the list of pickle files to load.
    deps = [os.path.join('expts', expt, f'{idx+1:06d}.gz')
            for idx in range(len(tasks))]

    # combine the results.
    results = []
    for dep in deps:
      # load the pickle file.
      with gzip.open(dep) as f:
        data = pickle.load(f)

      # add the loaded data into the results.
      results.append(data)

    # write the results to the target file.
    with gzip.open(targets[0], 'wb') as f:
      pickle.dump(tuple(results), f)

  # ---

  # loop over each experiment.
  for expt, tasks in expts.items():
    target = os.path.join('expts', f'{expt}.gz')
    task_names = [f'experiment:{expt}-{idx+1:06d}'
                  for idx in range(len(tasks))]

    yield {
      'name': expt,
      'actions': [(package, [expt])],
      'task_dep': task_names,
      'targets': [target]
    }

