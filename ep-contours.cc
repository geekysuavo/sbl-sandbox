
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

  /* declare a gamma utility class for transforming distributions. */
  gamma_util<1> fxi(alpha, beta);

  /* loop over all initial mean and variance values. */
  for (double i_var = 1e-6; i_var <= 1e6; i_var *= 2) {
    for (double i_mean = -10; i_mean <= 10; i_mean += 20) {
      double mean = i_mean;
      double var = i_var;

      /* trace the path of distributions from the initial values. */
      while (true) {
        const auto [new_mean, new_var] = fxi.map(mean, var);
        const bool feasible = (new_var > 0.001 && 1 / new_var > 1 / var);

        if (!feasible)
          break;

        std::cout << mean << " " << var << "\n";

        mean = new_mean;
        var = new_var;
      }

      std::cout << std::endl;
    }
  }
}

