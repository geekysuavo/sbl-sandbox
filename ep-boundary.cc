
/* Copyright (c) 2018 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

#include "gamma.hh"

int main (int argc, char **argv) {
  /* declare the gamma distribution parameters. */
  double alpha = 0.001;
  double beta = 0.001;

  /* if available, set the distribution parameters
   * from command-line arguments.
   */
  if (argc >= 2) alpha = atof(argv[1]);
  if (argc >= 3) beta  = atof(argv[2]);

  /* declare a gamma utility class for transforming distributions. */
  auto f = [alpha, beta] (double mean, double var) -> double {
    static gamma_util<1> fxi(alpha, beta);
    auto [m, v] = fxi.map(mean, var);
    return (1 / v) - (1 / var);
  };

  /* find the boundary on a uniformly-spaced grid of mean values. */
  for (double m = -10; m <= 10; m += 0.01) {
    /* initialize the lower variance bound. */
    double va = 1e-6;
    double fa = f(m, va);
    while (!isfinite(fa)) {
      va *= 1.01;
      fa = f(m, va);
    }

    /* initialize the upper variance bound. */
    double vb = 1e6;
    double fb = f(m, vb);

    /* declare the midpoint and its function value. */
    double vc, fc;

    /* ensure the lower bound is less than the upper. */
    if (fa > fb) {
      std::swap(va, vb);
      std::swap(fa, fb);
    }

    /* skip if no crossing exists. */
    if (fa > 0 && fb > 0) {
      std::cout << m << " " << std::min(va, vb) << "\n";
      continue;
    }

    /* find the intercept using the bisection method. */
    do {
      vc = 0.5 * (va + vb);
      fc = f(m, vc);

      if (fc < 0) {
        va = vc;
        fa = fc;
      }
      else {
        vb = vc;
        fb = fc;
      }
    }
    while (std::abs(va - vb) > 1e-6);

    /* write the intercept. */
    std::cout << m << " " << vc << "\n";
  }
}

