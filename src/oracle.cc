
/* Copyright (c) 2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

int main (int argc, char **argv) {
  /* initialize the problem instance. */
  instance_init(argc, argv);

  /* output the ground-truth signal with zero variance. */
  for (std::size_t i = 0; i < n; i++)
    std::cout << x0(i) << " " << 0 << "\n";
}

