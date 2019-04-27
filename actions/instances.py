
def instance(parms, targets):
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

  # define the set of acceptable parameters with their default values.
  options = {
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
    'thin':  (typed_size_t, 1)
  }

  # substitute all passed parameters into the options dictionary.
  for key, val in parms.items():
    (fn, default) = options[key]
    options[key] = (fn, val)

  # initialize the output line list with preprocessor statements.
  out = ['#pragma once', '#include <cstddef>']

  # loop over the options dictionary.
  for key in options:
    # get the type function and the value.
    (fn, val) = options[key]

    # write the parameter and its value.
    fn(out, key, val)

  # include the base instance code.
  out += ['#include "src/inst.hh"', '']

  # write the final output string to the target file.
  with open(targets[0], 'w') as f:
    f.write('\n'.join(out))

