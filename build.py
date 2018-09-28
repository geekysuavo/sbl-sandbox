
# import the required modules.
import sys, subprocess

# check for -v, --verbose.
verbose = False
for flag in ['-v', '--verbose']:
  if flag in sys.argv:
    verbose = True
    sys.argv.remove(flag)

# check for -h, --help
usage = False
for flag in ['-h', '--help']:
  if flag in sys.argv:
    usage = True
    sys.argv.remove(flag)

# check if the usage statement was requested.
if usage or len(sys.argv) not in (2, 3):
  print('Usage: {} solver [ident]'.format(sys.argv[0]))
  sys.exit()

# get the solver string.
solver = sys.argv[1]

# get the optional identifier string.
ident = ''
if len(sys.argv) == 3:
  ident = sys.argv[2]

# build some necessary strings.
source = '{}.cc'.format(solver)
binary = '{}{}'.format(solver, ident)
instance = 'instance{}.hh'.format(ident)

# make the list of arguments used for the unity build.
args = ['g++', '-std=c++17',
        '-I/usr/local/include/eigen3',
        '-include', instance,
        source, '-o', binary]

# if verbose output was requested, print the unity build command.
if verbose:
  print(' '.join(args))

# execute the unity build command.
proc = subprocess.run(args)

