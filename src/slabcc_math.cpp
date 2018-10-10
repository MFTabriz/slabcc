// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_math.hpp"

cube interp3(const rowvec& x, const rowvec& y, const rowvec& z, const cube& v, const rowvec& xi, const rowvec& yi, const rowvec& zi) {
	cube v_x = zeros(xi.n_elem, y.n_elem, z.n_elem);
	cube v_xy = zeros(xi.n_elem, yi.n_elem, z.n_elem);
	cube v_xyz = zeros(xi.n_elem, yi.n_elem, zi.n_elem);

	vec vec_i = zeros<vec>(x.n_elem);

	for (uword i = 0; i < y.n_elem; ++i) {
		for (uword j = 0; j < z.n_elem; ++j) {
			vec_i = vectorise(v(span(), span(i), span(j)));
			Spline<double> sp(x, vec_i);
			v_x(span(), span(i), span(j)) = sp.interpolate(xi);
		}
	}

	vec_i.set_size(y.n_elem);
	for (uword i = 0; i < xi.n_elem; ++i) {
		for (uword j = 0; j < z.n_elem; ++j) {
			vec_i = vectorise(v_x(span(i), span(), span(j)));
			Spline<double> sp(y, vec_i);
			v_xy(span(i), span(), span(j)) = sp.interpolate(yi);
		}
	}

	vec_i.set_size(z.n_elem);
	for (uword i = 0; i < xi.n_elem; ++i) {
		for (uword j = 0; j < yi.n_elem; ++j) {
			vec_i = vectorise(v_xy(span(i), span(j), span()));
			Spline<double> sp(z, vec_i);
			v_xyz(span(i), span(j), span()) = sp.interpolate(zi);
		}
	}
	return v_xyz;
}

cube interp3(const cube& v, const rowvec& xi, const rowvec& yi, const rowvec& zi) {
	return interp3(regspace<rowvec>(1, v.n_rows), regspace<rowvec>(1, v.n_cols), regspace<rowvec>(1, v.n_slices),
		v, xi, yi, zi);
}

tuple<cube, cube, cube> ndgrid(const rowvec& v1, const rowvec& v2, const rowvec& v3) {

	cube x(v1.n_elem, v2.n_elem, v3.n_elem, fill::ones), y = x, z = x;
	mat s_mat(v1.n_elem, v2.n_elem, fill::zeros);
	s_mat.each_col() = v1.t();
	x.each_slice() = s_mat;

	s_mat.each_row() = v2;
	y.each_slice() = s_mat;

	for (uword i = 0; i < v3.n_elem; ++i) {
		z.slice(i) *= v3(i);
	}

	return make_tuple(x, y, z);
}

tuple<cube, cube, cube> meshgrid(const rowvec& v1, const rowvec& v2, const rowvec& v3) {
	cube x2, y2, z2;
	tie(y2, x2, z2) = ndgrid(v2, v1, v3);
	return make_tuple(x2, y2, z2);
}

cube shift(cube cube_in, rowvec3 shifts) {
	if (cube_in.is_empty()) {
		return {};
	}

	shifts = round(rowvec3(SizeVec(cube_in) % shifts));
	for (uword i = 0; i < 3; ++i) {
		cube_in = shift(cube_in, shifts(i), i);
	}

	return cube_in;
}

vector<double> planar_average(const uword& direction, const cube& cube_in) {
	const uword dim_size = arma::size(cube_in)(direction);
	vector<double> average(dim_size, 0);
	for (uword i = 0; i < dim_size; ++i) {
		switch (direction) {
		case 0:
			average.at(i) = accu(cube_in(span(i), span(), span()));
			break;
		case 1:
			average.at(i) = accu(cube_in(span(), span(i), span()));
			break;
		case 2:
			average.at(i) = accu(cube_in(span(), span(), span(i)));
			break;
		}
	}
	return average;
}

cx_vec fft(vec X)
{
	//TODO: should come up with a better solution than reinterpret_cast
	cx_vec out(X.n_elem);
	fftw_plan plan = fftw_plan_dft_r2c_1d(X.n_elem, X.memptr(), reinterpret_cast<fftw_complex*>(out.memptr()), FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	for (uword i = out.n_elem / 2 + 1; i < out.n_elem; ++i)
		out(i) = conj(out(X.n_rows - i));

	return out;
}

cx_vec fft(cx_vec X)
{
	cx_vec out(X.n_elem);
	fftw_plan plan = fftw_plan_dft_1d(X.n_elem, reinterpret_cast<fftw_complex*>(X.memptr()), reinterpret_cast<fftw_complex*>(out.memptr()), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	return out;
}

cx_cube fft(cube X)
{
	cx_cube out(X.n_rows / 2 + 1, X.n_cols, X.n_slices);
	fftw_plan plan = fftw_plan_dft_r2c_3d(X.n_slices, X.n_cols, X.n_rows, X.memptr(), reinterpret_cast<fftw_complex*>(out.memptr()), FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);
	out.resize(X.n_rows, X.n_cols, X.n_slices);

	for (uword i = X.n_rows / 2 + 1; i < X.n_rows; ++i) {
		out(i, 0, 0) = conj(out(X.n_rows - i, 0, 0));
		for (uword j = 1; j < X.n_cols; ++j) {
			out(i, j, 0) = conj(out(X.n_rows - i, X.n_cols - j, 0));
		}
		for (uword k = 1; k < X.n_slices; ++k) {
			out(i, 0, k) = conj(out(X.n_rows - i, 0, X.n_slices - k));
		}
		for (uword j = 1; j < X.n_cols; ++j) {
			for (uword k = 1; k < X.n_slices; ++k) {
				out(i, j, k) = conj(out(X.n_rows - i, X.n_cols - j, X.n_slices - k));
			}
		}
	}

	return out;
}

cx_cube fft(cx_cube X)
{
	cx_cube fft(X.n_rows, X.n_cols, X.n_slices);
	fftw_plan plan = fftw_plan_dft_3d(X.n_slices, X.n_cols, X.n_rows, reinterpret_cast<fftw_complex*>(X.memptr()), reinterpret_cast<fftw_complex*>(fft.memptr()), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	return fft;
}

cx_vec ifft(cx_vec X)
{
	cx_vec out(X.n_elem);
	fftw_plan plan = fftw_plan_dft_1d(X.n_elem, reinterpret_cast<fftw_complex*>(X.memptr()), reinterpret_cast<fftw_complex*>(out.memptr()), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	return out / out.n_elem;
}

cx_cube ifft(cx_cube X)
{
	cx_cube ifft(X.n_rows, X.n_cols, X.n_slices);
	fftw_plan plan = fftw_plan_dft_3d(X.n_slices, X.n_cols, X.n_rows, reinterpret_cast<fftw_complex*>(X.memptr()), reinterpret_cast<fftw_complex*>(ifft.memptr()), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan);
	fftw_destroy_plan(plan);

	return ifft / X.n_elem;
}

SizeCube as_size(const urowvec3& vec) {
	return SizeCube(vec(0), vec(1), vec(2));
}
SizeMat as_size(const urowvec2& vec) {
	return SizeMat(vec(0), vec(1));
}

mat fmod(mat mat_in, const double& denom) noexcept {
	mat_in.for_each([&denom](double& val) noexcept {
		val = fmod(val, denom);
	});
	return mat_in;
}
mat fmod_p(mat mat_in, const double& denom) noexcept {
	mat_in.for_each([&denom](double& val) noexcept {
		val = fmod_p(val, denom);
	});
	return mat_in;
}

double fmod_p(double num, const double& denom) noexcept {
	num = fmod(num, denom);
	if (num < 0) num += abs(denom);
	return num;
}

void write_mat2file(const mat& input, const string& output_file) {
	ofstream out_file;
	out_file.open(output_file);
	out_file << fixed << showpos << setprecision(15);
	input.each_row([&out_file](const rowvec &row) { 
		row.for_each([&out_file](const double& val) { out_file << val << " "; });
		out_file << endl; 
	});
	out_file.close();
}

