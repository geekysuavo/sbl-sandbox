#!/usr/bin/env python3

# /* ... */
def typed_comment(L, key, val):
  pass

# constexpr bool
def typed_bool(L, key, val):
  sval = 'true' if val else 'false'
  L += [f'constexpr bool {key} = {sval};']

# constexpr std::size_t
def typed_size_t(L, key, val):
  L += [f'constexpr std::size_t {key} = {val};']

# constexpr double
def typed_double(L, key, val):
  L += [f'constexpr double {key} = {val};']

# set of acceptable instance parameters, their accompanying
# writer functions and default values.
instance_parms = {
  '1': (typed_comment, 'problem sizes'),
  'm': (typed_size_t, 50),
  'n': (typed_size_t, 100),

  '2':     (typed_comment, 'problem data initializers'),
  'unif':  (typed_bool, True),
  'k':     (typed_size_t, 10),
  'sigma': (typed_double, 0.005),
  'seed':  (typed_size_t, 47351),

  '3':      (typed_comment, 'weight prior parameters'),
  'alpha0': (typed_double, 0.001),
  'beta0':  (typed_double, 0.001),

  '4':       (typed_comment, 'noise prior parameters'),
  'nu0':     (typed_double, 0.001),
  'lambda0': (typed_double, 0.001),

  '5':     (typed_comment, 'algorithm parameters'),
  'iters': (typed_size_t, 250),
  'burn':  (typed_size_t, 10),
  'thin':  (typed_size_t, 1),

  'solvers': (typed_comment, ())
}

# defaults: dictionary of default parameter values.
defaults = {key: val[1] for key, val in instance_parms.items()}

# return a consistent dictionary of seed values.
def seeds():
  import random
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
            yield {**defaults,
                   'm': m, 'k': k, 'seed': seed,
                   'unif': unif, 'sigma': sigma,
                   'solvers': ('smf', 'mf', 'fmf', 'if')}

# convergence: yields instances with 10x longer iteration counts
# to compare convergence of each method.
#
def convergence():
  for seed in list(seeds())[:5]:
    for unif in (True, False):
      yield {**defaults, 'seed': seed, 'iters': 10000, 'burn': 1000,
             'solvers': ('gs', 'smf', 'mf', 'fmf', 'if')}

# combine groups of tasks into named experiments.
expts = {
  'phases': tuple(phases()),
  'converg': tuple(convergence())
}

# task_solve: task generator used by doit to manage experimental results.
def task_solve():
  # solve: task python action.
  def solve(ident, parms, targets):
    # import necessary modules.
    import subprocess
    import pickle
    import gzip
    import os

    # ensure the working directory exists.
    target = targets[0]
    wd = os.path.dirname(target)
    if not os.path.isdir(wd):
      os.makedirs(wd)

    # determine the header filename.
    header = f'{wd}/{ident}.hh'

    # build the header string.
    hstr = ['#pragma once', '#include <cstddef>']
    for parm, value in parms.items():
      fn = instance_parms[parm][0]
      fn(hstr, parm, value)

    # move to working with strings, include the core instance header.
    hstr = '\n'.join(hstr) + '\n'
    with open('src/inst.hh') as f:
      hstr += f.read()

    # write the header file.
    with open(header, 'wb') as f:
      f.write(hstr.encode('utf-8'))

    # prepare a dictionary to hold the results.
    solvers = parms['solvers']
    results = {s: None for s in solvers}

    # loop over each solution method.
    for solver in solvers:
      # determine the source and binary filenames.
      source = f'src/{solver}.cc'
      binary = f'{wd}/{ident}.{solver}'

      # execute the compilation command.
      args = ['g++', '-std=c++14',
              '-I.', '-I/usr/local/include/eigen3',
              '-include', header,
              source, '-o', binary]
      proc = subprocess.run(args)

      # execute the solver.
      proc = subprocess.run([binary], capture_output=True, text=True)

      # tupify: parse a string into a tuple of tuples.
      def tupify(lines):
        return tuple(tuple(float(field) for field in line.split())
                     for line in lines.strip().split('\n'))

      # parse standard output and standard error.
      out = tupify(proc.stdout)
      err = tupify(proc.stderr)

      # store the results.
      results[solver] = (out, err)

      # remove the binary file.
      os.remove(binary)

    # remove the header file.
    os.remove(header)

    # write the results to the target file.
    parms['results'] = results
    with gzip.open(target, 'wb') as f:
      pickle.dump(parms, f)

  # ---

  # yield tasks from each experiment.
  for expt, insts in expts.items():
    for idx in range(len(insts)):
      # get the instance parameters, especially (m, k).
      parms = insts[idx]
      m = parms['m']
      k = parms['k']

      # build the name and target filename strings.
      ident = f'{idx:06d}'
      name = f'{expt}{ident}'
      target = f'expts/{expt}/m{m}/k{k}/{ident}.gz'

      # yield the task data.
      yield {
        'name': name,
        'actions': [(solve, [ident, parms])],
        'targets': [target]
      }

# script entry point.
if __name__ == '__main__':
  import doit
  doit.run(globals())

