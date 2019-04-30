
# required imports.
import random

# return a consistent dictionary of seed values.
def seeds():
  num_seeds = 100
  random.seed(4735921)
  return {random.randint(99, 9999999) for i in range(num_seeds)}

# phases: yields a grid of instances from which phase diagrams
# can be (roughly) estimated.
#
def phases():
  for seed in seeds():
    for unif in (True, False):
      for sigma in (0, 0.001, 0.01, 0.1):
        for m in range(5, 95+1, 5):
          for k in range(5, m+1, 5):
            yield {'m': m, 'n': 100, 'k': k,
                   'seed': seed, 'unif': unif, 'sigma': sigma,
                   'methods': ('smf', 'mf', 'fmf', 'if')}

# convergence: yields instances with 10x longer iteration counts
# to compare convergence of each method.
#
def convergence():
  for seed in list(seeds())[:5]:
    for unif in (True, False):
      yield {'m': 50, 'n': 100,
             'seed': seed, 'iters': 10000, 'burn': 1000, 'thin': 10,
             'methods': ('gs', 'smf', 'mf', 'fmf', 'if')}

# combine groups of tasks into named experiments.
expts = {
  'phases': tuple(phases()),
  'converg': tuple(convergence())
}

# solvers: yields all unique (method, m, n) combinations required by
# the experiments.
def solvers():
  S = {}
  for tasks in expts.values():
    for task in tasks:
      S = {*S, *{(meth, task['m'], task['n'])
                 for meth in task['methods']}}

  return sorted(tuple(S))

# solve: run a solver for a set of parameters.
def solve(parms, meth):
  # required imports.
  import subprocess
  import os

  # tupify: parse a string into a tuple of tuples.
  def tupify(lines):
    return tuple(tuple(float(field) for field in line.split())
                 for line in lines.strip().split('\n'))

  # create a reduced parameter dictionary.
  P = dict(parms)
  del P['methods']
  (m, n) = (P.pop(key) for key in ('m', 'n'))

  # build the binary filename and arguments list.
  binary = os.path.join('bin', f'{meth}-{m}-{n}')
  args = [f'{key}={str(val).lower()}' for key, val in P.items()]

  # execute the solver.
  proc = subprocess.run([binary, *args], encoding='utf-8',
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)

  # parse standard output and standard error.
  out = tupify(proc.stdout)
  err = tupify(proc.stderr)

  # return the final result.
  return (out, err)

