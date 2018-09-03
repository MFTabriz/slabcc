// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include <armadillo>
#include <fftw3.h>
#include "general_io.hpp"
#include "spline.hpp"

#define PI datum::pi


//These functions extend the functionality of the included Armadillo library (version 8.5) by providing:
// (arma,iostream) << and >> operators for reading from iostreams to armadillo types and vice versa
// fft/ifft functions for cube files: Wrappers for FFTW with MATLAB/Octave scaling convention
//				in FFT functions the forward FFT scales by 1, and the reverse scales by 1/N. 
// 3D ndgrid and 3D meshgrid
// 3D spline interpolation
// 3D shift: by number of the elements along one axis or with a relative 3D shift vector
// 3D swap_axes
// irowvec to matrix/cube size conversion
// scalar triple product


using namespace arma;

// a (very inefficient) 3D spline interpolation!
cube interp3(const rowvec& x, const rowvec& y, const rowvec& z, const cube& v, const rowvec& xi, const rowvec& yi, const rowvec& zi);
cube interp3(const cube& v, const rowvec& xi, const rowvec& yi, const rowvec& zi);

//Rectangular grid in 3D space
tuple<cube, cube, cube> ndgrid(const rowvec& v1, const rowvec& v2, const rowvec& v3);

//3D meshgrid
tuple<cube, cube, cube> meshgrid(const rowvec& v1, const rowvec& v2, const rowvec& v3);

//shifts a cube by a relative 3D vector [0 1]
cube shift(cube cube_in, rowvec3 shifts);

//Planar average of a cube in defined direction
//direction: 0,1,2 > x,y,z
vector<double> planar_average(const uword& direction, const cube& cube_in);


// 1D FFT of complex data.
// no normalization for forward FFT
cx_vec fft(cx_vec X);

//1D FFT of real data.
// no normalization for forward FFT
cx_vec fft(vec X);


//3D FFT of complex data.
//no normalization for forward FFT
cx_cube fft(cx_cube X);

//3D FFT of real data.
//no normalization for forward FFT
cx_cube fft(cube X);

//1D inverse FFT of complex data.
//normalized by N = X.n_elem
cx_vec ifft(cx_vec X);

//3D inverse FFT of complex data.
//normalized by N = X.n_elem
cx_cube ifft(cx_cube X);

//returns a cube size object from the values inside a vector
SizeCube as_size(const urowvec3& vec);

//returns a matrix size object from the values inside a vector
SizeMat as_size(const urowvec2& vec);

//element-wise fmod
mat fmod(mat mat_in, const double& denom) noexcept;

//element-wise positive fmod
mat fmod_p(mat mat_in, const double& denom) noexcept;

//positive fmod
double fmod_p(double num, const double& denom) noexcept;

void write_mat2file(const mat& input, const string& output_file);

//generate a copy of the cube with the elements shifted by N positions along:
//dim=0: each row
//dim=1: each column
//dim=2: each slice
template <typename T>
Cube<T> shift(const Cube<T>& A, const sword& N, const uword& dim) {
	Cube<T> B(arma::size(A));
	auto index_init = regspace<uvec>(0, A.n_elem - 1);
	imat sub_shift = conv_to<imat>::from(ind2sub(arma::size(A), index_init));
	uword size = arma::size(A)(dim);
	sub_shift.row(dim).for_each([&N, &size](sword& i) noexcept {
		i += N;
		while (i < 0) i += size;
		i = i % size;
	});

	B(sub2ind(arma::size(A), conv_to<umat>::from(sub_shift))) = A(index_init);

	return B;
}

//Undo a fftshift
template <typename T>
Row<T> ifftshift(const Row<T> &A) {
	return shift(A, -1 * (A.n_elem / 2));
}

//Undo a fftshift
template <typename T>
Cube<T> ifftshift(Cube<T> A) {
	for (uword i = 0; i < 3; ++i) {
		A = shift(A, -1 * (arma::size(A)(i) / 2), i);
	}

	return A;
}

//returns the size of a cube as a rowvec
template<typename T>
urowvec3 SizeVec(const Cube<T>& c) {
	const SizeCube size = arma::size(c);
	return urowvec({ size(0), size(1), size(2) });
}

//returns the size of a matrix as a rowvec
template<typename T>
urowvec2 SizeVec(const Mat<T>& c) {
	SizeMat size = arma::size(c);
	return urowvec({ size(0), size(1) });
}

template <typename T>
Cube<T> swap_axes(const Cube<T> &cube_in, const uword& axis1, const uword& axis2) {

	if (cube_in.is_empty()) {
		return {};
	}
	if (axis1 == axis2) {
		return cube_in;
	}

	auto index_init = regspace<uvec>(0, cube_in.n_elem - 1);
	auto sub_swap = ind2sub(arma::size(cube_in), index_init);
	sub_swap.swap_rows(axis1, axis2);

	urowvec3 output_size = SizeVec(cube_in);
	output_size.swap_cols(axis1, axis2);
	Cube<T> cube_out(as_size(output_size), fill::zeros);

	cube_out(sub2ind(arma::size(cube_out), sub_swap)) = cube_in(index_init);

	return cube_out;
}


//single-line output for vec
template<typename T>
ostream &operator << (ostream &o, const Row<T> &vec) {
	vec.for_each([&o](const T &elem) { o << " " << elem; });
	return o;
}

//single-line output for mat
template<typename T>
ostream &operator << (ostream &o, const Mat<T> &mat) {
	mat.each_row([&o](const Row<T> &row) { o << row << ";"; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, Row<T> &vec) {
	vec.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, subview_row<T> vec) {
	vec.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, Cube<T> &c) {
	c.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, Mat<T> &m) {
	m.for_each([&o](T &elem) { o >> elem; });
	return o;
}

//sign of the val as -1/0/+1
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}
