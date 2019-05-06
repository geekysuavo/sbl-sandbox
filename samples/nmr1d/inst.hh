
/* Copyright (c) 2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

/* standard library headers. */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tuple>
#include <cmath>

/* fast fourier transform header. */
#include <fftw3.h>

/* define the data type of the input data table. */
using data_item = std::tuple<std::size_t, double, double>;
using data_table = std::vector<data_item>;

/* define the data types of real and complex vectors. */
constexpr std::size_t K = 2;
using real_vector = std::unique_ptr<double[]>;
using complex_vector = std::unique_ptr<double[][K]>;

/* fft: fast fourier transform class.
 */
template<std::size_t N>
struct fft {
public:
  /* fft(): allocate the data vector and plans. */
  fft () {
    X = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    Pfwd = fftw_plan_dft_1d(N, X, X, FFTW_FORWARD, FFTW_ESTIMATE);
    Pinv = fftw_plan_dft_1d(N, X, X, FFTW_BACKWARD, FFTW_ESTIMATE);
  }

  /* ~fft(): free all allocated fftw memory. */
  ~fft () {
    fftw_free(X);
    fftw_destroy_plan(Pfwd);
    fftw_destroy_plan(Pinv);
  }

  /* data(): return the data vector. */
  auto data () { return X; }

  /* fwd(), inv(): apply the forward, inverse transforms. */
  void fwd () { fftw_execute(Pfwd); }
  void inv () { fftw_execute(Pinv); }

private:
  fftw_plan Pfwd, Pinv;
  fftw_complex *X;
};

/* load(): construct a data table from a file.
 */
auto load (const std::string& filename) {
  std::ifstream ifs;
  ifs.open(filename, std::ifstream::in);

  std::size_t len = 0;
  ifs >> len;

  data_table v;
  v.reserve(len);

  for (std::size_t i = 0; i < len; i++) {
    std::size_t s = 0;
    double re = 0;
    double im = 0;

    ifs >> s;
    ifs >> re;
    ifs >> im;

    v.push_back(data_item{s, re, im});
  }

  ifs.close();
  return v;
}

/* schedule_vector(): construct a schedule vector from a data table.
 */
auto schedule_vector (const data_table& data, std::size_t n) {
  real_vector sched{new double[n]};

  for (std::size_t i = 0; i < n; i++)
    sched[i] = 0;

  for (auto &item : data)
    sched[std::get<0>(item)] = 1;

  return std::move(sched);
}

/* measured_vector(): construct a measured vector from a data table.
 */
auto measured_vector (const data_table& data, std::size_t n) {
  complex_vector y{new double[n][K]};

  for (std::size_t i = 0; i < n; i++)
    y[i][0] = y[i][1] = 0;

  for (auto& item : data) {
    y[std::get<0>(item)][0] = std::get<1>(item);
    y[std::get<0>(item)][1] = std::get<2>(item);
  }

  return std::move(y);
}

