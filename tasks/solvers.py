
def task_solvers():
  import os
  from util import expts
  from actions.solvers import solver

  # yield tasks from each experiment and solver combination.
  for expt, gen in expts.items():
    insts = tuple(gen())
    for sol in ('gs', 'smf', 'mf', 'if', 'fmf'):
      for idx in range(len(insts)):
        name = f'{expt}{idx:06d}'
        hdr = f'inst/{name}.hh'
        src = f'src/{sol}.cc'
        exe = f'bin/{sol}/{name}'
        parms = insts[idx]
        yield {
          'name': f'{sol}_{name}',
          'actions': [(solver, [parms, hdr, src])],
          'file_dep': [hdr, src],
          'targets': [exe]
        }

