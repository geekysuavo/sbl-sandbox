
/* Copyright (c) 2018 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

#pragma once
#include <random>
#include <iostream>
#include <Eigen/Dense>

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

/* instance_init(): initialize the current problem instance.
 */
static void instance_init () {
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
}

