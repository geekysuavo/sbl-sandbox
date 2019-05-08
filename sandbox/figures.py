
# required imports.
from doit.tools import run_once, create_folder
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import itertools
import pickle
import gzip
import os
import gc

# method_names: dictionary mapping internally used method names
# to the names used when plotting for publication.
method_names = {
  'smf': 'VRVM',
  'mf': 'SAVE',
  'if': 'IF-SBL',
  'gs': 'GS-SBL',
  'fmf': 'FMF-SBL'
}

# setparams: set some base parameters used for every figure.
def setparams():
  mpl.rcParams.update({'font.size': 8,
                       'legend.fontsize': 'small',
                       'lines.linewidth': 1})

# mean: compute the arithmetic mean of the elements of an iterable.
def mean(x):
  return sum(x) / len(x)

# median: compute the median of the elements of an iterable.
def median(x):
  return sorted(x)[len(x)//2]

# nmse: compute the normalized mean squared error between iterables.
def nmse(x, x0):
  ss = sum((a - b)**2 for a, b in zip(x, x0))
  ss0 = sum(b**2 for b in x0)
  return ss / ss0

# kldiv: compute the kl-divergence from uv to uv_hat.
def kldiv(uv, uv_hat):
  def kl(u, v, u_hat, v_hat):
    return (u - u_hat)**2 / v_hat + v / v_hat - np.log(v / v_hat) - 1

  return 0.5 * sum(kl(*args) for args in zip(*zip(*uv), *zip(*uv_hat)))

# get_nmse: get the normalized mean square error of the estimate
# for a single method, at a single instance.
def get_nmse(R, method, seed, m, k):
  unif = True

  L = [(r['results']['oracle'][0], r['results'][method][0])
       for r in R if r['unif'] == unif and r['seed'] == seed
       and r['m'] == m and r['k'] == k][0]

  l0 = tuple(l[0] for l in L[0])
  lm = tuple(l[0] for l in L[1])
  return nmse(lm, l0)

# get_kldiv: get the kl-divergence of the estimate for a single
# method, at a single instance.
def get_kldiv(R, method, seed, m, k):
  unif = True

  L = [(r['results']['gs'][0], r['results'][method][0])
       for r in R if r['unif'] == unif and r['seed'] == seed
       and r['m'] == m and r['k'] == k][0]

  (uv, uv_hat) = L
  return kldiv(uv, uv_hat)

# get_runtimes: get the (x,y) coordinates of the runtime curve
# for a single method.
def get_runtimes(R, method):
  L = sorted([(r['m'], mean(r['times'])) for r in R
              if r['method'] == method])

  (x, y) = zip(*L)
  return (np.array(x), np.array(y))

# get_objectives: get the (x,y) coordinates of the objective function
# curve for a single method.
def get_objectives(R, method):
  seeds = tuple({r['seed'] for r in R})
  seed = seeds[1]
  unif = True

  L = [r['results'][method][1] for r in R
       if r['unif'] == unif and r['seed'] == seed][0]

  (x, y) = zip(*L)
  return (np.array(x), np.array(y))

# get_phasediag: get the (x,y,z) coordinates of the phase diagram
# surface for a single method.
def get_phasediag(R, method):
  kv = sorted(tuple({r['k'] for r in R}))
  mv = sorted(tuple({r['m'] for r in R}))
  seeds = tuple({r['seed'] for r in R})

  (M, K) = np.meshgrid(mv, kv)
  P = np.zeros(M.shape)

  for k in kv:
    idx_k = kv.index(k)

    ms = sorted(tuple({r['m'] for r in R
                       if r['k'] == k
                       and r['unif'] == True
                       and r['seed'] == seeds[0]}))
    recs = tuple([mean([get_nmse(R, method, seed, m, k) < 0.01
                        for seed in seeds]) for m in ms])

    for i in range(len(ms)):
      idx_m = mv.index(ms[i])
      P[idx_m, idx_k] = recs[i]

  return (M, K, P)

# get_errors: get the (x,y) coordinates of the error curve
# for a single method, at a single k.
def get_errors(R, method, k):
  seeds = tuple({r['seed'] for r in R})
  ms = sorted(tuple({r['m'] for r in R
                     if r['k'] == k
                     and r['unif'] == True
                     and r['seed'] == seeds[0]}))

  L = [median([get_nmse(R, method, seed, m, k)
               for seed in seeds]) for m in ms]

  return (np.array(ms), np.array(L))

# get_divergences: get the (x,y) coordinates of the kl-divergence
# curve for a single method, at a single k.
def get_divergences(R, method, k):
  seeds = tuple({r['seed'] for r in R})
  ms = sorted(tuple({r['m'] for r in R
                     if r['k'] == k
                     and r['unif'] == True
                     and r['seed'] == seeds[0]}))

  L = [median([get_kldiv(R, method, seed, m, k)
               for seed in seeds]) for m in ms]

  return (np.array(ms), np.array(L))

# loadnmr: load nmr data from a text file.
def loadnmr(filename, n):
  def fields(F, n):
    ints = tuple(int(f) for f in F[:n])
    flts = tuple(float(f) for f in F[n:])
    return ints + flts

  with open(os.path.join('samples', f'nmr{n}d', filename)) as f:
    return [fields(line.strip().split(), n) for line in f.readlines()]

# ---

# figure1: task action for rendering the timing/convergence figure.
def figure1(targets):
  # load the timing result pickles.
  methods = ('smf', 'mf', 'if', 'fmf', 'gs', 'oracle')
  results = []
  for method in methods:
    filename = os.path.join('expts', 'timing', f'{method}.gz')
    with gzip.open(filename) as f:
      results += pickle.load(f)

  # collect the relevant timing information.
  runtimes = {method: get_runtimes(results, method)
              for method in methods}

  # free the pickle.
  results = None
  gc.collect()

  # load the convergence results pickle.
  with gzip.open(os.path.join('expts', 'convergence.gz')) as f:
    results = pickle.load(f)

  # collect the relevant convergence information.
  objectives = {method: get_objectives(results, method)
                for method in methods if method is not 'oracle'}

  # free the pickle.
  results = None
  gc.collect()

  # prepare the figure for plotting.
  setparams()
  fig, (left, right) = plt.subplots(ncols=2, figsize=(6, 2.75))

  # plot the timing panel data.
  x = runtimes['oracle'][0]
  Y = {method: runtimes[method][1] for method in methods}
  for method in methods:
    if method is not 'oracle':
      _ = left.semilogy(x, Y[method] - Y['oracle'],
                        label=method_names[method])

  # set the timing panel properties.
  _ = left.set_xlabel('Measurements')
  _ = left.set_ylabel('Runtime / s')
  _ = left.legend()
  _ = left.grid()

  # plot the convergence panel data.
  x = objectives['gs'][0]
  Y = {method: objectives[method][1] for method in methods
       if method is not 'oracle'}
  for method in methods:
    if method is not 'oracle':
      _ = right.semilogx(x, Y[method][:1000], label=method_names[method])

  # set the convergence panel properties.
  _ = right.set_xlabel('Iteration')
  _ = right.set_ylabel('Objective')
  _ = right.grid()

  # render the figure.
  fig.tight_layout()
  fig.savefig(targets[0], dpi=600, format='eps')


# figure2: task action for rendering the phase diagram figure.
def figure2(targets):
  # load the noiseless error results pickle.
  with gzip.open(os.path.join('expts', 'errors0.gz')) as f:
    results = pickle.load(f)

  # collect the relevant recovery information.
  methods = ('smf', 'fmf', 'if')
  diagrams = {method: get_phasediag(results, method)
              for method in methods}

  # free the pickle.
  results = None
  gc.collect()

  # prepare the figure for plotting.
  setparams()
  fig, (left, middle, right) = plt.subplots(ncols=3, figsize=(6, 2.75),
                                            sharex=True, sharey=True)

  # plot the left panel.
  method = 'smf'
  (Y, X, Z) = diagrams[method]
  _ = left.contourf(X, Y, Z, levels=10, cmap='jet')
  _ = left.set_title(method_names[method])
  _ = left.set_ylabel('Nonzeros')

  # plot the middle panel.
  method = 'fmf'
  (Y, X, Z) = diagrams[method]
  _ = middle.contourf(X, Y, Z, levels=10, cmap='jet')
  _ = middle.set_title(method_names[method])
  _ = middle.set_xlabel('Measurements')

  # plot the right panel.
  method = 'if'
  (Y, X, Z) = diagrams[method]
  _ = right.contourf(X, Y, Z, levels=10, cmap='jet')
  _ = right.set_title(method_names[method])

  # render the figure.
  fig.tight_layout()
  fig.savefig(targets[0], dpi=600, format='eps')


# figure3: task action for rendering the error figure.
def figure3(targets):
  # loop over two massive pickles. :(
  sigmas = (0.01, 0.1)
  klist = (20, 30, 40)
  methods = ('smf', 'mf', 'if', 'fmf', 'gs')
  errors = {sigma: None for sigma in sigmas}
  for sigma in sigmas:
    # load the error results pickle.
    with gzip.open(os.path.join('expts', f'errors{sigma}.gz')) as f:
      results = pickle.load(f)

    # collect the relevant error information.
    errors[sigma] = {method: {k: get_errors(results, method, k)
                              for k in klist} for method in methods}

    # free the pickle.
    results = None
    gc.collect()

  # prepare the figure for plotting.
  setparams()
  fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(6, 6),
                         sharex=True, sharey=True)

  # define a helper function to plot for a given (sigma, k).
  def plot_errs(ax, sigma, k):
    for method in methods:
      (x, y) = errors[sigma][method][k]
      _ = ax.plot(x, y, label=method_names[method])

    _ = ax.set_title(f'$\sigma = {sigma}, k = {k}$')
    _ = ax.grid()

    if k == min(klist):
      _ = ax.legend()

    if sigma == min(sigmas) and k == median(klist):
      _ = ax.set_ylabel('Median NMSE')

    if k == max(klist):
      _ = ax.set_xlabel('Measurements')

  # plot the results.
  for rowidx in range(len(klist)):
    for colidx in range(len(sigmas)):
      plot_errs(ax[rowidx, colidx], sigmas[colidx], klist[rowidx])

  # render the figure.
  fig.tight_layout()
  fig.savefig(targets[0], dpi=600, format='eps')


# figure4: task action for rendering the divergence figure.
def figure4(targets):
  # load an error results pickle.
  klist = (20, 30, 40)
  methods = ('smf', 'mf', 'if', 'fmf')
  with gzip.open(os.path.join('expts', 'errors0.001.gz')) as f:
    results = pickle.load(f)

  # collect the relevant recovery information.
  divs = {method: {k: get_divergences(results, method, k)
                   for k in klist} for method in methods}

  # free the pickle.
  results = None
  gc.collect()

  # prepare the figure for plotting.
  setparams()
  fig, (top, middle, bottom) = plt.subplots(nrows=3, figsize=(2.5, 6),
                                            sharex=True, sharey=True)

  # define a helper function to plot for a given k.
  def plot_divs(ax, k):
    for method in methods:
      (x, y) = divs[method][k]
      _ = ax.semilogy(x, y, label=method_names[method])

    if k == max(klist):
      _ = ax.legend()

    _ = ax.set_title(f'k = {k}')
    _ = ax.grid()

  # plot the top panel (k=20).
  plot_divs(top, 20)

  # plot the middle panel (k=30).
  plot_divs(middle, 30)
  _ = middle.set_ylabel('Median KL(p||q)')

  # plot the bottom panel (k=40).
  plot_divs(bottom, 40)
  _ = bottom.set_xlabel('Measurements')

  # render the figure.
  fig.tight_layout()
  fig.savefig(targets[0], dpi=600, format='eps')


# figure5: task action for rendering the nmr2d sample figure.
def figure5(targets):
  # load the relevant data.
  truth = loadnmr('truth.dat', 2)
  ista = loadnmr('ista.out', 2)
  fmf = loadnmr('fmf.out', 2)

  # prepare the figure for plotting.
  setparams()
  fig, (left, middle, right) = plt.subplots(ncols=3, figsize=(6, 2.25),
                                            sharex=True, sharey=True)

  # define the shape and level generation function.
  shape = (1024, 2048)
  def levs(a, b, n=10):
    return np.logspace(np.log10(a), np.log10(b), n)

  # plot the ground truth panel.
  (x, y, z) = zip(*truth)
  X = np.array(x).reshape(shape)
  Y = np.array(y).reshape(shape)
  Z = np.array(z).reshape(shape)
  _ = left.contour(X, Y, Z, levels=levs(0.02, 0.55))
  _ = left.set_title('Ground truth')

  # plot the ista panel.
  (x, y, z) = zip(*ista)
  X = np.array(x).reshape(shape)
  Y = np.array(y).reshape(shape)
  Z = np.array(z).reshape(shape)
  _ = middle.contour(X, Y, Z, levels=levs(0.07, 1.5))
  _ = middle.set_title('ISTA')

  # plot the fmf-sbl panel.
  (x, y, z, dz) = zip(*fmf)
  X = np.array(x).reshape(shape)
  Y = np.array(y).reshape(shape)
  Z = np.array(z).reshape(shape)
  _ = right.contour(X, Y, Z, levels=levs(0.01, 1))
  _ = right.set_title('FMF-SBL')

  # set the panel limits.
  _ = plt.xlim(200, 840)
  _ = plt.ylim(500, 2000)

  # drop tick labels.
  fmt = mpl.ticker.NullFormatter()
  _ = left.xaxis.set_major_formatter(fmt)
  _ = left.yaxis.set_major_formatter(fmt)
  _ = middle.xaxis.set_major_formatter(fmt)
  _ = right.xaxis.set_major_formatter(fmt)

  # render the figure.
  fig.tight_layout()
  fig.savefig(targets[0], dpi=600, format='eps')


# figure6: task action for rendering the nmr1d sample figure.
def figure6(targets):
  # load the relevant data.
  truth = loadnmr('truth.dat', 1)
  ista = loadnmr('ista.out', 1)
  fmf = loadnmr('fmf.out', 1)

  # prepare the figure for plotting.
  setparams()
  fig, (top, middle, bottom) = plt.subplots(nrows=3, figsize=(2.5, 6),
                                            sharex=True)
  
  # plot the ground truth panel.
  (x, y, _) = zip(*truth)
  _ = top.plot(x, np.fft.fftshift(y))
  _ = top.set_title('Ground truth')
  _ = top.grid()

  # plot the ista panel.
  (x, y, _) = zip(*ista)
  _ = middle.plot(x, np.fft.fftshift(y))
  _ = middle.set_title('ISTA')
  _ = middle.grid()

  # plot the fmf-sbl panel.
  (x, y, _, dy) = zip(*fmf)
  lb = np.fft.fftshift(y - np.sqrt(dy))
  ub = np.fft.fftshift(y + np.sqrt(dy))
  _ = bottom.fill_between(x, lb, ub, color=(0.50, 0.80, 1.00))
  _ = bottom.plot(x, np.fft.fftshift(y))
  _ = bottom.set_title('FMF-SBL')
  _ = bottom.grid()

  # set figure properties.
  _ = plt.xlim(0, 2047)
  _ = bottom.set_xlabel('Frequency')
  _ = middle.set_ylabel('Amplitude')
  _ = plt.xticks([0, 1024, 2047], ['$-\pi$', '0', '$\pi$'])

  # render the figure.
  fig.tight_layout()
  fig.savefig(targets[0], dpi=600, format='eps')

# ---

# task_figure: task generator for rendering final figures.
def task_figure():
  GST = globals()
  fig_dir = 'figures'
  for n in range(1, 7):
    yield {
      'name': str(n),
      'actions': [(create_folder, [fig_dir]),
                  (GST[f'figure{n}'])],
      'uptodate': [run_once],
      'targets': [os.path.join(fig_dir, f'fig{n}.eps')]
    }

