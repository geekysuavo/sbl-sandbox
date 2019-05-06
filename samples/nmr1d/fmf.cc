
/* Copyright (c) 2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 *
 * Compilation:
 *   g++ -std=c++14 -O3 fmf.cc -o fmf -lfftw3 -lm
 */

#include "inst.hh"

int main () {
  /* read the instance data from disk. */
  auto data = load("inst.dat");

  /* problem sizes. */
  const std::size_t m = data.size();
  const std::size_t n = 2048;

  /* prepare fftw. */
  fft<n> F;
  auto dy = F.data();

  /* number of iterations. */
  const std::size_t iters = 1000;

  /* weight prior parameters. */
  const double alpha0 = 0.001;
  const double beta0 = 0.001;

  /* noise prior parameters. */
  const double nu0 = 50.0;
  const double lambda0 = 0.005;

  /* posterior shape parameter values. */
  const double alpha = alpha0 + 0.5;
  const double nu = nu0 + 0.5 * m;
  const double a = double(m) / double(n);

  /* construct the schedule and measured
   * vectors from the data table.
   */
  auto S = schedule_vector(data, n);
  auto y = measured_vector(data, n);

  /* allocate u, z, v. */
  complex_vector u{new double[n][K]};
  complex_vector z{new double[n][K]};
  real_vector v{new double[n]};

  /* allocate xi, initialize tau. */
  real_vector xi{new double[n]};
  double tau = nu0 / lambda0;

  /* initialize u, v, xi. */
  for (std::size_t i = 0; i < n; i++) {
    /* initialize u, v. */
    u[i][0] = u[i][1] = 0;
    v[i] = 1;

    /* initialize xi. */
    xi[i] = alpha0 / beta0;
  }

  /* update the error vector: dy == A(u) - y
   *  dy <- u
   *  ifft(dy)
   *  dy <- S .* (dy * sqrt(n) - y)
   */
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < K; k++)
      dy[i][k] = u[i][k];
  F.inv();
  for (std::size_t i = 0; i < n; i++)
    for (std::size_t k = 0; k < K; k++)
      dy[i][k] = S[i] * (dy[i][k] / std::sqrt(n) - y[i][k]);

  /* iterate. */
  for (std::size_t it = 0; it < iters; it++) {
    /* update z. */
    F.fwd();
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t k = 0; k < K; k++)
        z[i][k] = tau * (u[i][k] - dy[i][k] / std::sqrt(n));

    /* update v, u, xi. */
    double trAAV = 0;
    for (std::size_t i = 0; i < n; i++) {
      /* update v. */
      v[i] = 1 / (xi[i] + tau * a);
      trAAV += a * v[i];

      /* update u. */
      const double w = 1 / (tau + xi[i]);
      for (std::size_t k = 0; k < K; k++)
        u[i][k] = w * z[i][k];

      /* update xi. */
      const double u2 = std::pow(u[i][0], 2) + std::pow(u[i][1], 2);
      xi[i] = alpha / (beta0 + 0.5 * (u2 + v[i]));
    }

    /* update the error vector. */
    for (std::size_t i = 0; i < n; i++)
      for (std::size_t k = 0; k < K; k++)
        dy[i][k] = u[i][k];
    F.inv();
    double ess = 0;
    for (std::size_t i = 0; i < n; i++) {
      for (std::size_t k = 0; k < K; k++) {
        dy[i][k] = S[i] * (dy[i][k] / std::sqrt(n) - y[i][k]);
        ess += std::pow(dy[i][k], 2);
      }
    }

    /* update tau. */
    tau = nu / (lambda0 + 0.5 * (ess + trAAV));
  }

  /* output the results. */
  std::cout.precision(9);
  std::cout << std::scientific;
  for (std::size_t i = 0; i < n; i++)
    std::cout << i << " "
              << u[i][0] << " "
              << u[i][1] << " "
              << v[i] << "\n";
}

