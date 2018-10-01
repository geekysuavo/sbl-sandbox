
/* Copyright (c) 2018 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

#include "gamma.hh"

int main (int argc, char **argv) {
  /* declare the gamma distribution parameters. */
  double alpha = 0.001;
  double beta = 0.001;

  /* if available, set the distribution parameters
   * from command-line arguments.
   */
  if (argc >= 2) alpha = atof(argv[1]);
  if (argc >= 3) beta  = atof(argv[2]);

  /* define a function that maps a univariate normal distribution
   * into another using the F- and G-functions.
   */
  auto map = [alpha, beta] (double mean, double var) 
           -> std::pair<double, double> {
    static Gamma<0> fxi(alpha, beta);

    const double F = fxi.F(1 / var);
    const double G = fxi.G(1 / var);

    const double new_mean = mean * F;
    const double new_var = var * F + std::pow(mean, 2) * (G - F * F);

    return {new_mean, new_var};
  };

  /* loop over all initial mean and variance values. */
  for (double i_var = 1e2; i_var <= 1e6; i_var *= 2) {
    for (double i_mean = -10; i_mean <= 10; i_mean += 20) {
      double mean = i_mean;
      double var = i_var;

      /* trace the path of distributions from the initial values. */
      bool feasible = false;
      do {
        auto [new_mean, new_var] = map(mean, var);
        feasible = (new_var > 0.001 && 1 / new_var > 1 / var);

        std::cout << mean << " " << var << "\n";

        mean = new_mean;
        var = new_var;
      }
      while (feasible);

      std::cout << std::endl;
    }
  }
}

