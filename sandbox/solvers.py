
# required imports.
from doit.tools import create_folder
from sandbox.util import solvers
import subprocess
import os

# task_solver: task generator for building all needed
# solver binaries from source.
def task_solver():
  # solver: task action function for task_solver().
  def solver(meth, m, n, targets):
    # get the source, binary and (temporary) header filenames.
    source = os.path.join('src', f'{meth}.cc')
    binary = targets[0]
    header = f'{binary}.hh'

    # write the header that defines (m, n).
    with open(header, 'w') as f:
      f.write('\n'.join([
        '#pragma once',
        '#include <cstddef>',
        f'constexpr std::size_t m = {m};',
        f'constexpr std::size_t n = {n};',
        '#include "src/inst.hh"',
        '']))

    # compile the binary and remove the header file.
    args = ['g++', '-std=c++14', '-O3', '-I.', '-Ieigen3',
            '-include', header, source, '-o', binary]
    proc = subprocess.run(args)
    os.remove(header)

  # ---

  # loop over each required solver.
  for meth, m, n in solvers():
    name = f'{meth}-{m}-{n}'
    target = os.path.join('bin', name)
    header = os.path.join('src', 'inst.hh')
    source = os.path.join('src', f'{meth}.cc')

    yield {
      'name': name,
      'actions': [(create_folder, ['bin']),
                  (solver, [meth, m, n])],
      'file_dep': [source, header],
      'task_dep': ['eigen3'],
      'targets': [target]
    }

