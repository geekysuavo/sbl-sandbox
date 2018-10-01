
/* Copyright (c) 2018 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

#pragma once

#include <cmath>
#include <random>
#include <cstddef>
#include <stdexcept>
#include <boost/math/special_functions/gamma.hpp>

/* gamma_util: class holding extra functionalities of gamma distributions
 * not available from std::gamma_distribution<double>.
 */
template<std::size_t n>
class gamma_util {
public:
  gamma_util () : gamma_util(1, 1) {}

  gamma_util (double a, double b) : _alpha(0), _beta(0), gen(), rdist() {
    /* set the distribution parameters. */
    set(a, b);

    /* seed the distribution with a random number. */
    std::random_device rdev;
    gen.seed(rdev());
  }

  /* alpha(), beta(): parameter getter functions.
   */
  double alpha () const { return _alpha; }
  double beta () const { return _beta; }

  /* alpha(double), beta(double): parameter setter functions.
   */
  void alpha (double a) { set(a, _beta); }
  void beta (double b) { set(_alpha, b); }

  /* set(): main parameter setter function.
   */
  void set (double a, double b) {
    /* bounds-check the parameters. */
    if (a <= 0) throw std::domain_error("alpha must be positive");
    if (b <= 0) throw std::domain_error("beta must be positive");

    /* update the parameters of the random distribution. */
    typename decltype(rdist)::param_type p{a, 1 / b};
    rdist.param(p);

    /* store the new parameter values. */
    _alpha = a;
    _beta = b;

    /* update the quadrature rule. */
    compute_rule();
  }

  /* draw(): sample a value from the distribution.
   */
  double draw () {
    return rdist(gen);
  }

  /* mean(): return the first moment, E[x], of the distribution.
   */
  double mean () const {
    return _alpha / _beta;
  }

  /* var(): return the second moment, E[x^2], of the distribution.
   */
  double var () const {
    return _alpha / (_beta * _beta);
  }

  /* F(): compute the F-function for the distribution,
   *
   *  F(x) := E[(1 + theta / x)**(-1)]
   *
   * where theta ~ Gamma(a,b) and x >= 0.
   */
  double F (double x) const {
    if (x < 0)
      throw std::domain_error("x must be non-negative");

    const double s = 1 - _alpha;
    const double z = _beta * x;

    return std::exp(z) * std::pow(z, _alpha) * boost::math::tgamma(s, z);
  }

  /* G(): compute the G-function for the distribution,
   *
   *  G(x) := E[(1 + theta / x)**(-2)]
   *
   * where theta ~ Gamma(a,b) and x >= 0.
   */
  double G (double x) const {
    const double s = 1 - _alpha;
    const double z = _beta * x;

    return (s - z) * F(x) + z;
  }

  /* weight(): function to get quadrature weights.
   */
  double weight (std::size_t k) const {
    if (k >= n) throw std::out_of_range("k");
    return w[k];
  }

  /* node(): function to get quadrature nodes.
   */
  double node (std::size_t k) const {
    if (k >= n) throw std::out_of_range("k");
    return x[k];
  }

  /* expect(): function to compute quadrature estimates of expected
   * values of random functions.
   */
  template<typename F>
  auto expect (const F& f) const {
    auto I = w[0] * f(x[0]);
    for (std::size_t k = 1; k < n; k++)
      I += w[k] * f(x[k]);

    return I;
  }

private:
  /* distribution parameters:
   *  @_alpha: shape.
   *  @_beta: rate.
   */
  double _alpha, _beta;

  /* random sampling members:
   *  @gen: pseudorandom number generator.
   *  @rdist: random distribution.
   */
  std::default_random_engine gen;
  std::gamma_distribution<double> rdist;

  /* quadrature members:
   *  @w: weights.
   *  @x: nodes.
   */
  double w[n], x[n];

  /* fsgn(): compute the sign function of a double-precision float.
   */
  static inline double fsgn (double x) {
    return (x >= 0 ? 1 : -1);
  }

  /* imtqlx(): implicit ql algorithm for symmetric tridiagonal matrices.
   *
   * arguments:
   *  @u: matrix main diagonal.
   *  @v: matrix sub-diagonal.
   *  @w: pointer to the output weight (eigenvalue) array.
   *  @x: pointer to the output node (eigenvector first element) array.
   */
  static void imtqlx (double u[n], double v[n], double w[n], double x[n]) {
    /* define the maximum number of iterations. */
    constexpr std::size_t itmax = 30;

    /* get the machine epsilon. */
    constexpr auto eps = std::numeric_limits<double>::epsilon();

    /* return if a one-element matrix was given. */
    if (n == 1)
      return;

    v[n-1] = 0;

    for (std::size_t l = 1; l <= n; l++) {
      std::size_t it = 0;

      while (1) {
        std::size_t m;
        for (m = l; m <= n; m++) {
          if (m == n)
            break;

          if (std::abs(v[m-1]) <= eps * (std::abs(u[m-1]) + std::abs(u[m])))
            break;
        }

        double p = u[l-1];

        if (m == l)
          break;

        /* fail if the iteration limit is exceeded. */
        it++;
        if (it > itmax)
          throw std::runtime_error("jacobi diagonalization failed");

        double g = (u[l] - p) / (2 * v[l-1]);
        double r = std::sqrt(g * g + 1);
        g = u[m-1] - p + v[l-1] / (g + std::abs(r) * fsgn(g));
        double s = 1;
        double c = 1;
        p = 0;

        for (std::size_t ii = 1; ii <= m - l; ii++) {
          std::size_t i = m - ii;
          double f = s * v[i-1];
          double b = c * v[i-1];

          if (std::abs(g) <= std::abs(f)) {
            c = g / f;
            r = std::sqrt(c * c + 1);
            v[i] = f * r;
            s = 1 / r;
            c = c * s;
          }
          else {
            s = f / g;
            r = std::sqrt(s * s + 1);
            v[i] = g * r;
            c = 1 / r;
            s = s * c;
          }

          g = u[i] - p;
          r = (u[i-1] - g) * s + 2 * c * b;
          p = s * r;
          u[i] = g + p;
          g = c * r - b;

          f = w[i];
          w[i] = s * w[i-1] + c * f;
          w[i-1] = c * w[i-1] - s * f;
        }

        u[l-1] = u[l-1] - p;
        v[l-1] = g;
        v[m-1] = 0;
      }
    }

    /* sort by the eigenvalues, in ascending order. */
    for (std::size_t ii = 2; ii <= n; ii++) {
      std::size_t i = ii - 1;
      std::size_t k = i;
      double p = u[i-1];

      for (std::size_t j = ii; j <= n; j++) {
        if (u[j-1] < p) {
          p = u[j-1];
          k = j;
        }
      }

      if (k != i) {
        u[k-1] = u[i-1];
        u[i-1] = p;

        p = w[i-1];
        w[i-1] = w[k-1];
        w[k-1] = p;
      }
    }

    /* copy the eigenvalues into their output vector. */
    for (std::size_t i = 0; i < n; i++)
      x[i] = u[i];
  }

  /* compute_rule(): compute the quadrature rule for a distribution.
   */
  void compute_rule () {
    /* declare temporary working variables. */
    double W[n] = {0}, X[n] = {0}, u[n] = {0}, v[n] = {0};

    /* compute the zero-order moment. */
    const double Z = std::exp(std::lgamma(_alpha));
    W[0] = std::sqrt(Z);

    /* compute the jacobi matrix elements. */
    for (std::size_t i = 0; i < n; i++) {
      u[i] = 2 * i + _alpha;
      v[i] = std::sqrt((i + 1) * (i + _alpha));
    }

    /* diagonalize the jacobi matrix. */
    imtqlx(u, v, W, X);

    /* finalize and store the quadrature rule. */
    for (std::size_t i = 0; i < n; i++) {
      w[i] = (W[i] * W[i]) / Z;
      x[i] = X[i] / _beta;
    }
  }
};

