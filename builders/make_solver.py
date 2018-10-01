
# import the required modules.
import subprocess

def make_solver(argv = (), usage = False, verbose = False):
  # check if the usage statement should be printed.
  if usage or len(argv) not in (1, 2):
    print('Arguments: solver-name [identifier]')
    return

  # get the solver string.
  solver = argv[0]

  # get the optional identifier string.
  ident = ''
  if len(argv) == 2:
    ident = argv[1]

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

