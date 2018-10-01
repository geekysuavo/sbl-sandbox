#!/usr/bin/env python3

# import the required modules.
import sys
from builders.make_clean import *
from builders.make_solver import *
from builders.make_instance import *

# define the set of scripts that can be executed.
scripts = {
  'clean': make_clean,
  'solver': make_solver,
  'instance': make_instance
}

# print_usage(): print a usage statement and exit.
def print_usage(prog):
  print('Usage: {} script [script-args]'.format(prog))
  sys.exit(1)

# check for -v, --verbose.
verbose = False
for flag in ('-v', '--verbose'):
  if flag in sys.argv:
    verbose = True
    sys.argv.remove(flag)

# check for -h, --help
usage = False
for flag in ('-h', '--help'):
  if flag in sys.argv:
    usage = True
    sys.argv.remove(flag)

# check that a builder script was specified.
if len(sys.argv) < 2:
  print_usage(sys.argv[0])

# get the name of the builder script.
script = sys.argv[1]

# check if the specified script is available.
if script not in scripts:
  print_usage(sys.argv[0])

# execute the builder function.
fn = scripts[script]
fn(sys.argv[2:], usage, verbose)

