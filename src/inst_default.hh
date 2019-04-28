
/* Copyright (c) 2018-2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

#pragma once
#include <cstddef>

/* problem sizes:
 *  @m: measurement count.
 *  @n: feature count.
 */
constexpr std::size_t m = 50;
constexpr std::size_t n = 100;

/* problem data initializers:
 *  @unif: whether or not to use uniform-amplitude impulses.
 *  @k: number of impulses in the weight vector.
 *  @sigma: measurement noise standard deviation.
 *  @seed: pseudorandom number generator seed.
 */
constexpr bool unif = true;
constexpr std::size_t k = 10;
constexpr double sigma = 0.001;
constexpr std::size_t seed = 47351;

/* weight prior parameters:
 *  @alpha0: shape.
 *  @beta0: rate.
 */
constexpr double alpha0 = 0.001;
constexpr double beta0 = 0.001;

/* noise prior parameters:
 *  @nu0: shape.
 *  @lambda0: rate.
 */
constexpr double nu0 = 0.001;
constexpr double lambda0 = 0.001;

/* algorithm parameters:
 *  @iters: iteration count.
 *  @burn: (monte carlo) burn-in iteration count.
 *  @thin: (monte carlo) thinning iteration count.
 */
constexpr std::size_t iters = 250;
constexpr std::size_t burn = 10;
constexpr std::size_t thin = 1;

/* include the base instance functions. */
#include "src/inst.hh"

