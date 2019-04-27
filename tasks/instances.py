
def task_instances():
  import os
  from util import expts
  from actions.instances import instance

  # yield tasks from each experiment.
  for expt, gen in expts.items():
    insts = tuple(gen())
    for idx in range(len(insts)):
      name = f'{expt}{idx:06d}'
      hdr = f'inst/{name}.hh'
      parms = insts[idx]
      yield {
        'name': name,
        'actions': [(instance, [parms])],
        'uptodate': [os.path.exists(hdr)],
        'targets': [hdr]
      }

