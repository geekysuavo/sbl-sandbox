
/* Copyright (c) 2019 Bradley Worley <geekysuavo@gmail.com>
 * Released under the MIT License.
 */

/* standard library headers. */
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <cmath>

/* fast fourier transform header. */
#include <fftw3.h>

/* define the data type of the input data table. */
using size_array = std::array<std::size_t, 2>;
using data_item = std::tuple<std::size_t, std::size_t,
                             double, double, double, double>;
using data_table = std::vector<data_item>;

/* define the data types of real and complex vectors. */
constexpr std::size_t K = 4;
using real_vector = std::unique_ptr<double[]>;
using complex_vector = std::unique_ptr<double[][K]>;

/* fft2_complex: scalar type used by the fft2 class. */
typedef double fft2_complex[K];

/* fft2: two-dimensional fast fourier transform class.
 */
template<std::size_t N0, std::size_t N1>
struct fft2 {
public:
  /* fft2(): allocate the data vector and plans. */
  fft2 () {
    fftw_init_threads();
    fftw_plan_with_nthreads(4);

    X = (fft2_complex*) fftw_malloc(sizeof(fft2_complex) * N0 * N1);
    fftw_iodim dims[2], vdims[2];
    double *xp = (double*) X;

    dims[0].n = N0;
    dims[0].is = dims[0].os = K;

    dims[1].n = N1;
    dims[1].is = dims[1].os = K * N0;

    vdims[0].n = N1;
    vdims[0].is = vdims[0].os = K * N0;

    vdims[1].n = N0;
    vdims[1].is = vdims[1].os = K;

    /* construct the forward plans. */
    Pfwd[0] = fftw_plan_guru_split_dft(1, dims, 1, vdims,
                                       xp + 0, xp + 1,
                                       xp + 0, xp + 1,
                                       FFTW_ESTIMATE);

    Pfwd[1] = fftw_plan_guru_split_dft(1, dims, 1, vdims,
                                       xp + 2, xp + 3,
                                       xp + 2, xp + 3,
                                       FFTW_ESTIMATE);

    Pfwd[2] = fftw_plan_guru_split_dft(1, dims + 1, 1, vdims + 1,
                                       xp + 0, xp + 2,
                                       xp + 0, xp + 2,
                                       FFTW_ESTIMATE);

    Pfwd[3] = fftw_plan_guru_split_dft(1, dims + 1, 1, vdims + 1,
                                       xp + 1, xp + 3,
                                       xp + 1, xp + 3,
                                       FFTW_ESTIMATE);

    /* construct the inverse plans. */
    Pinv[0] = fftw_plan_guru_split_dft(1, dims, 1, vdims,
                                       xp + 1, xp + 0,
                                       xp + 1, xp + 0,
                                       FFTW_ESTIMATE);

    Pinv[1] = fftw_plan_guru_split_dft(1, dims, 1, vdims,
                                       xp + 3, xp + 2,
                                       xp + 3, xp + 2,
                                       FFTW_ESTIMATE);

    Pinv[2] = fftw_plan_guru_split_dft(1, dims + 1, 1, vdims + 1,
                                       xp + 2, xp + 0,
                                       xp + 2, xp + 0,
                                       FFTW_ESTIMATE);

    Pinv[3] = fftw_plan_guru_split_dft(1, dims + 1, 1, vdims + 1,
                                       xp + 3, xp + 1,
                                       xp + 3, xp + 1,
                                       FFTW_ESTIMATE);
  }

  /* ~fft2(): free all allocated fftw memory. */
  ~fft2 () {
    fftw_free(X);

    for (std::size_t k = 0; k < K; k++) {
      fftw_destroy_plan(Pfwd[k]);
      fftw_destroy_plan(Pinv[k]);
    }

    fftw_cleanup_threads();
  }

  /* data(): return the data vector. */
  auto data () { return X; }

  /* fwd(): apply the forward transform. */
  void fwd () {
    double *xp = (double*) X;
    fftw_execute_split_dft(Pfwd[0], xp + 0, xp + 1, xp + 0, xp + 1);
    fftw_execute_split_dft(Pfwd[1], xp + 2, xp + 3, xp + 2, xp + 3);
    fftw_execute_split_dft(Pfwd[2], xp + 0, xp + 2, xp + 0, xp + 2);
    fftw_execute_split_dft(Pfwd[3], xp + 1, xp + 3, xp + 1, xp + 3);
  }

  /* inv(): apply the inverse transform. */
  void inv () {
    double *xp = (double*) X;
    fftw_execute_split_dft(Pinv[0], xp + 1, xp + 0, xp + 1, xp + 0);
    fftw_execute_split_dft(Pinv[1], xp + 3, xp + 2, xp + 3, xp + 2);
    fftw_execute_split_dft(Pinv[2], xp + 2, xp + 0, xp + 2, xp + 0);
    fftw_execute_split_dft(Pinv[3], xp + 3, xp + 1, xp + 3, xp + 1);
  }

private:
  fftw_plan Pfwd[K], Pinv[K];
  fft2_complex *X;
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
    std::size_t si = 0;
    std::size_t sj = 0;
    double rr = 0;
    double ri = 0;
    double ir = 0;
    double ii = 0;

    ifs >> si;
    ifs >> sj;
    ifs >> rr;
    ifs >> ri;
    ifs >> ir;
    ifs >> ii;

    v.push_back(data_item{si, sj, rr, ri, ir, ii});
  }

  ifs.close();
  return v;
}

/* schedule_vector(): construct a schedule vector from a data table.
 */
auto schedule_vector (const data_table& data, const size_array& N) {
  const std::size_t n = N[0] * N[1];
  real_vector sched{new double[n]};

  for (std::size_t i = 0; i < n; i++)
    sched[i] = 0;

  for (auto& item : data) {
    const std::size_t si = std::get<0>(item);
    const std::size_t sj = std::get<1>(item);
    const std::size_t idx = si + N[0] * sj;

    sched[idx] = 1;
  }

  return std::move(sched);
}

/* measured_vector(): construct a measured vector from a data table.
 */
auto measured_vector (const data_table& data, const size_array& N) {
  const std::size_t n = N[0] * N[1];
  complex_vector y{new double[n][K]};

  for (std::size_t i = 0; i < n; i++)
    y[i][0] = y[i][1] = 0;

  for (auto& item : data) {
    const std::size_t si = std::get<0>(item);
    const std::size_t sj = std::get<1>(item);
    const std::size_t idx = si + N[0] * sj;

    y[idx][0] = std::get<2>(item);
    y[idx][1] = std::get<3>(item);
    y[idx][2] = std::get<4>(item);
    y[idx][3] = std::get<5>(item);
  }

  return std::move(y);
}

