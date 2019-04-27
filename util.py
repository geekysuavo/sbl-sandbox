
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
            yield {'m': m, 'k': k, 'seed': seed,
                   'unif': unif, 'sigma': sigma}

# convergence: yields instances with 10x longer iteration counts
# to compare convergence of each method.
#
def convergence():
  for seed in list(seeds())[:5]:
    for unif in (True, False):
      yield {'seed': seed, 'iters': 10000}

# combine groups of tasks into named experiments.
expts = {
  'phases': phases,
  'converg': convergence
}

