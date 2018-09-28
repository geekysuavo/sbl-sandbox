
/* Copyright (c) 2018 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

int main () {
  /* initialize the problem instance. */
  instance_init();

  /* compute the diagonal of the measurement matrix gramian. */
  Eigen::Matrix<double, n, n> AtA = A.transpose() * A;
  Eigen::Matrix<double, n, 1> a = AtA.diagonal();

  /* compute the projected data vector. */
  Eigen::Matrix<double, n, 1> Aty = A.transpose() * y;

  /* initialize the weight means and variances. */
  Eigen::Matrix<double, n, 1> u, v;
  u.setZero();
  v.setOnes();

  /* initialize the precision means. */
  Eigen::Matrix<double, n, 1> xi;
  xi.setConstant(alpha / beta);

  /* initialize the noise mean. */
  double tau = nu / lambda;

  /* compute the updated shape parameters. */
  const double t_alpha = alpha + 0.5;
  const double t_nu = nu + 0.5 * m;

  /* iterate. */
  for (std::size_t it = 0; it < iters; it++) {
    /* update each weight factor. */
    for (std::size_t i = 0; i < n; i++) {
      /* update the weight variance. */
      v(i) = 1 / (xi(i) + tau * a(i));

      /* update the weight mean. */
      const double AtAx = AtA.row(i).dot(u) - a(i) * u(i);
      u(i) = tau * v(i) * (Aty(i) - AtAx);
    }

    /* update the precision means. */
    for (std::size_t i = 0; i < n; i++) {
      const double ex2 = std::pow(u(i), 2) + v(i);
      const double t_beta = beta + 0.5 * ex2;
      xi(i) = t_alpha / t_beta;
    }

    /* update the noise. */
    const double mess = (y - A * u).squaredNorm();
    const double trAAV = (AtA * v.asDiagonal()).trace();
    const double t_lambda = lambda + 0.5 * mess + 0.5 * trAAV;
    tau = t_nu / t_lambda;
  }

  /* output the final means. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << u(i) << (i + 1 == n ? "\n" : " ");

  /* output the final variances. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << v(i) << (i + 1 == n ? "\n" : " ");
}

