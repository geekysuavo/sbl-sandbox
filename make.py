#!/usr/bin/env python3

# task generator imports.
from sandbox.eigen import task_eigen3
from sandbox.solvers import task_solver
from sandbox.timings import task_timing
from sandbox.packages import task_package
from sandbox.experiments import task_experiment
from sandbox.samples import task_nmr1d, task_nmr2d

# script entry point.
if __name__ == '__main__':
  import doit
  doit.run(globals())

