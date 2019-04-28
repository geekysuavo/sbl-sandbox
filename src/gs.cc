
/* Copyright (c) 2018-2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

int main () {
  /* initialize the problem instance. */
  instance_init();

  /* initialize the weight means, variances, and running variables. */
  Eigen::Matrix<double, n, 1> M1, M2, mu, s2;
  M1.setZero();
  M2.setZero();
  mu.setZero();
  s2.setZero();

  /* declare variables for sampling tau, xi. */
  std::gamma_distribution<double> gam;
  using gam_param_t = typename decltype(gam)::param_type;
  auto gamrnd = [&] (double a, double b) -> double {
    gam_param_t gpar{a, 1 / b};
    gam.param(gpar);
    return gam(gen);
  };

  /* initialize the noise sample. */
  double tau = 0;

  /* initialize the precision sample. */
  Eigen::Matrix<double, n, 1> xi;
  xi.setZero();

  /* initialize the weight sample. */
  Eigen::Matrix<double, n, 1> x;
  x = A.transpose() * y;

  /* declare variables for sampling new weights. */
  std::normal_distribution<double> nrm{0, 1};
  Eigen::Matrix<double, n, 1> z2, z3, u;
  Eigen::Matrix<double, m, 1> z1, t;

  /* compute the updated shape parameters. */
  const double alpha = alpha0 + 0.5;
  const double nu = nu0 + 0.5 * m;

  /* randomly re-seed the pseudorandom number generator. */
  std::random_device rdev;
  gen.seed(rdev());

  /* iterate. */
  for (std::size_t it = 0; it < burn + iters; it++) {
    /* compute the error sum of squares. */
    const double ess = (y - A * x).squaredNorm();

    /* sample tau. */
    const double lambda = lambda0 + 0.5 * ess;
    tau = gamrnd(nu, lambda);

    /* sample xi. */
    for (std::size_t j = 0; j < n; j++) {
      const double beta = beta0 + 0.5 * std::pow(x(j), 2);
      xi(j) = gamrnd(alpha, beta);
    }

    /* draw an m-vector of standard normal variates. */
    for (std::size_t i = 0; i < m; i++)
      z1(i) = nrm(gen);

    /* draw an n-vector of standard normal variates. */
    for (std::size_t j = 0; j < n; j++)
      z2(j) = nrm(gen);

    /* z3 = sqrt(tau) .* A' * z1 + sqrt(xi) .* z2 */
    z3 = std::sqrt(tau) * A.transpose() * z1 +
         xi.cwiseSqrt().cwiseProduct(z2);

    /* u = tau .* A' * y + z3
     * t = (eye(m) ./ tau + A * diag(1 ./ xi) * A') \ (A * (u ./ xi))
     * x = (u - A' * t) ./ xi
     */
    u = tau * A.transpose() * y + z3;
    t = (Eigen::Matrix<double, m, m>::Identity() / tau +
         A * xi.cwiseInverse().asDiagonal() * A.transpose())
        .llt().solve(A * u.cwiseQuotient(xi));
    x = (u - A.transpose() * t).cwiseQuotient(xi);

    /* check if the sample should be stored. */
    if (it >= burn && (it - burn) % thin == 0) {
      /* compute the sample count. */
      const double itt = (it - burn) / thin;

      /* update the first moment. */
      M1 = x - mu;
      mu += M1 / itt;

      /* update the second moment. */
      M2 += M1.cwiseProduct(x - mu);
      s2 = M2 / itt;
    }

    /* FIXME: output the current likelihood value. */
    std::cerr << it << " " << 0.0 << "\n";
  }

  /* output the final means and variances. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << mu(i) << " " << s2(i) << "\n";
}

