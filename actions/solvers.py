
def solver(parms, header, source, targets):
  import subprocess
  import os

  # get the binary filename and directory name.
  binary = targets[0]
  bindir = os.path.dirname(binary)

  # ensure the directory exists.
  if not os.path.isdir(bindir):
    os.makedirs(bindir)

  # build the argument list for a compilation command.
  args = ['g++', '-std=c++14',
          '-I.', '-I/usr/local/include/eigen3',
          '-include', header,
          source, '-o', binary]

  # execute the compilation command.
  proc = subprocess.run(args)

