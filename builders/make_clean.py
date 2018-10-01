
# import the required modules.
import os, glob

def make_clean(argv = (), usage = False, verbose = False):
  # check if the usage statement should be printed.
  if usage or len(argv) not in (0,):
    print('No arguments permitted')
    return

  # get a list of source files and their basenames.
  sources = glob.glob('*.cc')
  names = [os.path.splitext(src)[0]
           for src in sources]

  # get a list of instance header files and their identifiers.
  instances = glob.glob('instance*.hh')
  idents = [os.path.splitext(i)[0][len('instance'):]
            for i in instances]

  # build a list of targets to be removed.
  targets = [n + i for n in names for i in idents]

  # remove each target.
  for target in targets:
    if os.path.exists(target):
      os.remove(target)
      if verbose:
        print('removed "{}"'.format(target))

