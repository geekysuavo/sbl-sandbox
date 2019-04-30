
/* Copyright (c) 2018-2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

#pragma once
#include <random>
#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>

/* problem data initializers:
 *  @unif: whether or not to use uniform-amplitude impulses.
 *  @k: number of impulses in the weight vector.
 *  @sigma: measurement noise standard deviation.
 *  @seed: pseudorandom number generator seed.
 */
bool unif = true;
std::size_t k = 10;
double sigma = 0.001;
std::size_t seed = 47351;

/* weight prior parameters:
 *  @alpha0: shape.
 *  @beta0: rate.
 */
double alpha0 = 0.001;
double beta0 = 0.001;

/* noise prior parameters:
 *  @nu0: shape.
 *  @lambda0: rate.
 */
double nu0 = 0.001;
double lambda0 = 0.001;

/* algorithm parameters:
 *  @iters: iteration count.
 *  @burn: (monte carlo) burn-in iteration count.
 *  @thin: (monte carlo) thinning iteration count.
 */
std::size_t iters = 1000;
std::size_t burn = 10;
std::size_t thin = 1;

/* instance-constant expressions:
 *  @alpha: weight posterior shape parameter.
 *  @nu: noise posterior shape parameter.
 */
constexpr double pi = 3.14159265358979323846264338327950288;
double alpha;
double nu;

/* global pseudorandom number generator:
 */
std::default_random_engine gen;

/* instance data:
 *  @A: measurement matrix.
 *  @x0: true weight vector.
 *  @y: data vector.
 */
Eigen::Matrix<double, m, n> A;
Eigen::Matrix<double, n, 1> x0;
Eigen::Matrix<double, m, 1> y;

/* precomputed data:
 *  @a: vector of diagonal elements of the measurement gramian.
 *  @L: twice the maximal eigenvalue of the measurement gramian.
 *  @phi0: constant offset term of the universal sbl objective.
 */
Eigen::Matrix<double, n, 1> a;
double L, phi0;

/* instance_init(): initialize the current problem instance.
 */
static void instance_init (int argc, char **argv) {
  /* parse runtime arguments. */
  for (std::size_t i = 1; i < argc; i++) {
    /* get the current argument. */
    std::string arg(argv[i]);
    auto idx = arg.find_first_of('=');
    if (idx == std::string::npos)
      continue;

    /* split the argument into key=val. */
    auto key = arg.substr(0, idx);
    auto val = arg.substr(idx + 1);

    /* run some dirty argument parsing. */
    if (key.compare("unif") == 0) {
      if (val.compare("true") == 0)       unif = true;
      else if (val.compare("false") == 0) unif = false;
    }
    else if (key.compare("k")     == 0) { k = std::stoi(val); }
    else if (key.compare("sigma") == 0) { sigma = std::stod(val); }
    else if (key.compare("seed")  == 0) { seed = std::stoull(val); }
    else if (key.compare("alpha0") == 0) { alpha0 = std::stod(val); }
    else if (key.compare("beta0")  == 0) { beta0 = std::stod(val); }
    else if (key.compare("nu0")     == 0) { nu0 = std::stod(val); }
    else if (key.compare("lambda0") == 0) { lambda0 = std::stod(val); }
    else if (key.compare("iters") == 0) { iters = std::stoul(val); }
    else if (key.compare("burn")  == 0) { burn = std::stoul(val); }
    else if (key.compare("thin")  == 0) { thin = std::stoul(val); }
  }

  /* update the posterior shape parameters. */
  alpha = alpha0 + 0.5;
  nu = nu0 + 0.5 * m;

  /* prepare to sample from two distributions:
   *  @idx: {0, 1, 2, ..., n-1}.
   *  @nrm: N(0, 1).
   */
  std::uniform_int_distribution<std::size_t> idx{0, n - 1};
  std::uniform_int_distribution<std::size_t> bin{0, 1};
  std::normal_distribution<double> nrm{0, 1};
  gen.seed(seed);

  /* compute each row of the measurement matrix. */
  for (std::size_t i = 0; i < m; i++) {
    /* sample the row elements from a standard normal. */
    for (std::size_t j = 0; j < n; j++)
      A(i,j) = nrm(gen);

    /* normalize the row to unit length. */
    A.row(i).normalize();
  }

  /* fill the feature vector with spikes. */
  std::size_t spikes = 0;
  x0.setZero();
  do {
    /* sample a new element index. */
    std::size_t j = idx(gen);
    if (x0(j) != 0)
      continue;

    /* sample a random spike intensity. */
    double xj = 0;
    if (unif)
      xj = (bin(gen) ? 1 : -1);
    else
      xj = nrm(gen);

    /* store the spike intensity. */
    x0(j) = xj;
    spikes++;
  }
  while (spikes < k);

  /* compute the noise-free data vector. */
  y = A * x0;

  /* check if the noise is nonzero. */
  if (sigma > 0) {
    /* add noise to the data vector. */
    for (std::size_t i = 0; i < m; i++)
      y(i) += sigma * nrm(gen);
  }

  /* construct an eigenvalue solver for the measurement matrix gramian. */
  Eigen::Matrix<double, n, n> AtA = A.transpose() * A;
  Eigen::EigenSolver<decltype(AtA)> es(AtA, false);

  /* get the diagonal elements of the gramian. */
  a = AtA.diagonal();

  /* identify the maximal eigenvalue of the gramian. */
  L = es.eigenvalues()(0).real();
  for (std::size_t j = 1; j < n; j++)
    L = std::max(L, es.eigenvalues()(j).real());

  /* double the result (should be twice the maximal eigenvalue). */
  L *= 2;

  /* compute the constant offset to the sbl objective. */
  phi0 = 0.5 * (m + n) * std::log(2 * pi)
       - nu0 * std::log(lambda0)
       - n * alpha0 * std::log(beta0)
       - (std::lgamma(nu) - std::lgamma(nu0))
       - n * (std::lgamma(alpha) - std::lgamma(alpha0));
}

