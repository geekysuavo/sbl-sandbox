
/* Copyright (c) 2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

int main (int argc, char **argv) {
  /* initialize the problem instance. */
  instance_init(argc, argv);

  /* initialize the weight means, variances, and bounding variables. */
  Eigen::Matrix<double, n, 1> u, v, z;
  u.setZero();
  v.setOnes();

  /* initialize the precision means. */
  Eigen::Matrix<double, n, 1> xi;
  xi.setConstant(alpha0 / beta0);

  /* initialize the noise mean. */
  double tau = nu0 / lambda0;

  /* compute the updated shape parameters. */
  const double alpha = alpha0 + 0.5;
  const double nu = nu0 + 0.5 * m;

  /* iterate. */
  for (std::size_t it = 0; it < iters; it++) {
    /* initialize the objective. */
    double phi = phi0;

    /* compute the new "right-hand-side" vector. */
    z = 0.5 * L * tau * u - tau * A.transpose() * (A * u - y);

    /* update each weight factor. */
    for (std::size_t i = 0; i < n; i++) {
      /* update the weight variance. */
      v(i) = 1 / (xi(i) + tau * a(i));

      /* update the weight mean. */
      u(i) = z(i) / (0.5 * L * tau + xi(i));

      /* update the objective. */
      phi -= 0.5 * std::log(v(i));
    }

    /* update the precision means. */
    for (std::size_t i = 0; i < n; i++) {
      /* update xi(i). */
      const double ex2 = std::pow(u(i), 2) + v(i);
      const double beta = beta0 + 0.5 * ex2;
      xi(i) = alpha / beta;

      /* update the objective. */
      phi += alpha * std::log(beta);
    }

    /* update the noise. */
    const double trAAV = a.dot(v);
    const double mess = (y - A * u).squaredNorm();
    const double lambda = lambda0 + 0.5 * mess + 0.5 * trAAV;
    tau = nu / lambda;

    /* output the objective. */
    phi += nu * std::log(lambda);
    std::cerr << it << " " << phi << "\n";
  }

  /* output the final means and variances. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << u(i) << " " << v(i) << "\n";
}

