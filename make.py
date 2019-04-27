#!/usr/bin/env python3

# task generator imports.
from tasks import task_instances, task_solvers

# script entry point.
if __name__ == '__main__':
  import doit
  doit.run(globals())
