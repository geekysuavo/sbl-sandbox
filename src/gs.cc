
/* Copyright (c) 2018-2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

int main (int argc, char **argv) {
  /* initialize the problem instance. */
  instance_init(argc, argv);

  /* declare variables for sampling tau, xi. */
  std::gamma_distribution<double> gam;
  using gam_param_t = typename decltype(gam)::param_type;
  auto gamrnd = [&] (double a, double b) -> double {
    gam_param_t gpar{a, 1 / b};
    gam.param(gpar);
    return gam(gen);
  };

  /* complete log-conditionals of tau, xi. */
  auto lgampdf = [] (double z, double a, double b) -> double {
    return a * std::log(b) - std::lgamma(a) + (a - 1) * std::log(z) - b * z;
  };

  /* initialize the weight means, variances, and running variables. */
  Eigen::Matrix<double, n, 1> M1, M2, xbar, s2;
  xbar.setZero();
  s2.setZero();
  M1.setZero();
  M2.setZero();

  /* initialize the precision (posterior) means. */
  Eigen::Matrix<double, n, 1> xibar;
  double taubar = 0;
  xibar.setZero();

  /* initialize the posterior marginal estimates. */
  double ptau = 0;
  double pxi = 0;
  double px = 0;

  /* initialize the noise sample. */
  double tau = 0;

  /* initialize the precision sample. */
  Eigen::Matrix<double, n, 1> beta, xi;
  beta.setZero();
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
      const double beta_j = beta0 + 0.5 * std::pow(x(j), 2);
      xi(j) = gamrnd(alpha, beta_j);
      beta(j) = beta_j;
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
    auto lltQ = (Eigen::Matrix<double, m, m>::Identity() / tau +
                 A * xi.cwiseInverse().asDiagonal() * A.transpose()).llt();
    t = lltQ.solve(A * u.cwiseQuotient(xi));
    x = (u - A.transpose() * t).cwiseQuotient(xi);

    /* check if the sample should be stored. */
    if (it >= burn && (it - burn) % thin == 0) {
      /* compute the sample count. */
      const double itt = (it - burn) / thin + 1;

      /* update the first moment of x. */
      M1 = x - xbar;
      xbar += M1 / itt;

      /* update the second moment of x. */
      M2 += M1.cwiseProduct(x - xbar);
      s2 = M2 / itt;

      /* update the posterior mean estimate of tau. */
      taubar += (tau - taubar) / itt;

      /* update the posterior mean estimate of xi. */
      xibar += (xi - xibar) / itt;

      /* compute the log-determinant of the covariance matrix
       * of the complete conditional of x.
       */
      double lndetS = m * std::log(tau);
      for (std::size_t i = 0; i < m; i++)
        lndetS += 2 * std::log(lltQ.matrixL()(i,i));
      for (std::size_t j = 0; j < n; j++)
        lndetS += std::log(xi(j));

      /* compute the complete log-conditional of x. */
      const double Q1 = tau * (A * M1).squaredNorm();
      const double Q2 = M1.transpose() * xi.cwiseProduct(M1);
      const double Q = -0.5 * (Q1 + Q2) - 0.5 * lndetS
                       - 0.5 * n * std::log(2 * pi);

      /* compute the reduced complete conditional of x. (no log!) */
      px += (std::exp(Q) - px) / itt;

      /* compute the reduced complete log-conditional of tau. */
      const double essbar = (y - A * xbar).squaredNorm();
      const double l = lambda0 + 0.5 * essbar;
      ptau = lgampdf(taubar, nu, l);

      /* compute the reduced complete log-conditional of xi. */
      pxi = 0;
      for (std::size_t j = 0; j < n; j++) {
        const double b = beta0 + 0.5 * std::pow(xbar(j), 2);
        pxi += lgampdf(xibar(j), alpha, b);
      }

      /* log-likelihood: ln p(y|xbar,taubar). */
      double lml = 0;
      lml += -0.5 * taubar * essbar
             + 0.5 * m * std::log(taubar / (2 * pi));

      /* conditional weight prior: ln p(xbar|xibar). */
      lml += -0.5 * xbar.transpose() * xibar.cwiseProduct(xbar)
             - 0.5 * n * std::log(2 * pi);
      for (std::size_t j = 0; j < n; j++)
        lml += 0.5 * std::log(xibar(j));

      /* weight precision prior: ln p(xibar). */
      for (std::size_t j = 0; j < n; j++)
        lml += lgampdf(xibar(j), alpha0, beta0);

      /* noise precision prior: ln p(taubar). */
      lml += lgampdf(taubar, nu0, lambda0);

      /* finally, subtract the log-posterior ordinal estimate. */
      const double lpo = std::log(px) + pxi + ptau;
      lml -= lpo;

      /* output the current negated log-marginal likelihood estimate. */
      std::cerr << itt << " " << -lml << "\n";
    }
  }

  /* output the final means and variances. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << xbar(i) << " " << s2(i) << "\n";
}

