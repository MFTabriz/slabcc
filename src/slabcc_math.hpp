// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include "arma_io.hpp"
#include "general_io.hpp"
#include "slabcc_consts.hpp"
#include "spline.hpp"
#include <armadillo>
#include <fftw3.h>

// These functions extend the functionality of the included Armadillo library by
// providing:
// (arma,iostream) << and >> operators for reading from iostreams to armadillo
// types and vice versa fft/ifft functions for cube files: Wrappers for FFTW
// with MATLAB/Octave scaling convention in FFT functions the forward FFT scales
// by 1, and the reverse scales by 1/N. 3D ndgrid and 3D meshgrid 3D spline
// interpolation 3D shift: by number of the elements along one axis or with a
// relative 3D shift vector irowvec to matrix/cube size conversion scalar triple
// product

// a (very inefficient) 3D spline interpolation!
arma::cube interp3(const arma::rowvec &x, const arma::rowvec &y,
                   const arma::rowvec &z, const arma::cube &v,
                   const arma::rowvec &xi, const arma::rowvec &yi,
                   const arma::rowvec &zi);
arma::cube interp3(const arma::cube &v, const arma::rowvec &xi,
                   const arma::rowvec &yi, const arma::rowvec &zi);

// Rectangular grid in 3D space
std::tuple<arma::cube, arma::cube, arma::cube>
ndgrid(const arma::rowvec &v1, const arma::rowvec &v2, const arma::rowvec &v3);

// 3D meshgrid
std::tuple<arma::cube, arma::cube, arma::cube> meshgrid(const arma::rowvec &v1,
                                                        const arma::rowvec &v2,
                                                        const arma::rowvec &v3);

// shifts a cube by a relative 3D vector [0 1]
arma::cube shift(arma::cube cube_in, arma::rowvec3 shifts);

// Planar average of a cube in the defined direction
// direction: 0,1,2 > x,y,z
arma::vec planar_average(const arma::uword &direction,
                         const arma::cube &cube_in);

// 1D FFT of complex data.
// no normalization for forward FFT
arma::cx_vec fft(arma::cx_vec X);

// 1D FFT of real data.
// no normalization for forward FFT
arma::cx_vec fft(arma::vec X);

// 3D FFT of complex data.
// no normalization for forward FFT
arma::cx_cube fft(arma::cx_cube X);

// 3D FFT of real data.
// no normalization for forward FFT
arma::cx_cube fft(arma::cube X);

// 1D inverse FFT of complex data.
// normalized by N = X.n_elem
arma::cx_vec ifft(arma::cx_vec X);

// 3D inverse FFT of complex data.
// normalized by N = X.n_elem
arma::cx_cube ifft(arma::cx_cube X);

// returns a cube size object from the values inside a vector
arma::SizeCube as_size(const arma::urowvec3 &vec);

// returns a matrix size object from the values inside a vector
arma::SizeMat as_size(const arma::urowvec2 &vec);

// element-wise fmod
arma::mat fmod(arma::mat mat_in, const double &denom) noexcept;

// element-wise positive fmod
arma::mat fmod_p(arma::mat mat_in, const double &denom) noexcept;

// positive fmod
double fmod_p(double num, const double &denom) noexcept;

// just a simple square! May cause overflows!!
inline double square(const double &input) noexcept { return input * input; }

// Poisson solver in 3D with anisotropic dielectric profiles
// diel is the N*3 matrix of variations in dielectric tensor elements in
// direction normal to the surface
arma::cx_cube poisson_solver_3D(const arma::cx_cube &rho, arma::mat diel,
                                arma::rowvec3 lengths,
                                arma::uword normal_direction);

// generate a copy of the cube with the elements shifted by N positions along:
// dim=0: each row
// dim=1: each column
// dim=2: each slice
template <typename T>
arma::Cube<T> shift(const arma::Cube<T> &A, const arma::sword &N,
                    const arma::uword &dim) {
  arma::Cube<T> B(arma::size(A));
  const auto index_init = arma::regspace<arma::uvec>(0, A.n_elem - 1);
  arma::imat sub_shift =
      arma::conv_to<arma::imat>::from(arma::ind2sub(arma::size(A), index_init));
  const arma::uword size = arma::size(A)(dim);
  sub_shift.row(dim).for_each([&N, &size ](arma::sword & i) noexcept {
    i += N;
    while (i < 0)
      i += size;
    i = i % size;
  });

  B(arma::sub2ind(arma::size(A), arma::conv_to<arma::umat>::from(sub_shift))) =
      A(index_init);

  return B;
}

// Undo a fftshift
template <typename T> arma::Row<T> ifftshift(const arma::Row<T> &A) {
  return arma::shift(A, -1 * (A.n_elem / 2));
}

// Undo a fftshift
template <typename T> arma::Cube<T> ifftshift(arma::Cube<T> A) {
  for (arma::uword i = 0; i < 3; ++i) {
    A = arma::shift(A, -1 * (arma::size(A)(i) / 2), i);
  }

  return A;
}

// returns the size of a cube as a rowvec
template <typename T> arma::urowvec3 SizeVec(const arma::Cube<T> &c) {
  const arma::SizeCube size = arma::size(c);
  return arma::urowvec({size(0), size(1), size(2)});
}

// returns the size of a matrix as a rowvec
template <typename T> arma::urowvec2 SizeVec(const arma::Mat<T> &c) {
  const arma::SizeMat size = arma::size(c);
  return arma::urowvec({size(0), size(1)});
}

// sign of the val as -1/0/+1
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
