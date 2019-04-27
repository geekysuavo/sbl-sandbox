
/* Copyright (c) 2018-2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

int main () {
  /* initialize the problem instance. */
  instance_init();

  /* initialize intermediate quantities. */
  Eigen::Matrix<double, m, m> Q, K;
  Eigen::Matrix<double, m, 1> t;
  Eigen::Matrix<double, n, 1> S;

  /* initialize the weight means and variances. */
  Eigen::Matrix<double, n, 1> mu, sigma;
  mu.setZero();
  sigma.setOnes();

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

    /* K = A * diag(1 ./ xi) * A'
     * Q = inv(eye(m) ./ tau + K)
     */
    K = A * xi.cwiseInverse().asDiagonal() * A.transpose();
    auto L = (Eigen::Matrix<double, m, m>::Identity() / tau + K).llt();
  //Q = (Eigen::Matrix<double, m, m>::Identity() / tau + K).inverse();
    Q = L.solve(Eigen::Matrix<double, m, m>::Identity());

    /* update the weight means. */
    t = Q * K * y;
    mu = tau * (A.transpose() * (y - t)).cwiseQuotient(xi);

    /* update the weight variances. */
    S = (A.transpose() * Q * A).diagonal();
    sigma = xi.cwiseInverse() - S.cwiseQuotient(xi.cwiseAbs2());

    /* compute the log-determinant of the kernel matrix. */
    double lndetQ = 0;
    for (std::size_t j = 0; j < m; j++)
      lndetQ += 2 * std::log(L.matrixL()(j,j));

    /* compute the log-determinant of the full covariance matrix. */
    double lndetGamma = lndetQ + m * std::log(tau);
    for (std::size_t i = 0; i < n; i++)
      lndetGamma += std::log(xi(i));

    /* update the objective. */
    phi -= 0.5 * lndetGamma;

    /* update the precision means. */
    for (std::size_t i = 0; i < n; i++) {
      const double mu2 = std::pow(mu(i), 2) + sigma(i);
      const double beta = beta0 + 0.5 * mu2;
      xi(i) = alpha / beta;

      /* update the objective. */
      phi += alpha * std::log(beta);
    }

    /* update the noise. */
    const double mess = (y - A * mu).squaredNorm();
    const double trAAS = K.trace() - (K * Q * K).trace();
    const double lambda = lambda0 + 0.5 * mess + 0.5 * trAAS;
    tau = nu / lambda;

    /* output the objective. */
    phi += nu * std::log(lambda);
    std::cerr << it << " " << phi << "\n";
  }

  /* output the final means and variances. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << mu(i) << " " << sigma(i) << "\n";
}

