// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_math.hpp"

arma::cube interp3(const arma::rowvec &x, const arma::rowvec &y,
                   const arma::rowvec &z, const arma::cube &v,
                   const arma::rowvec &xi, const arma::rowvec &yi,
                   const arma::rowvec &zi) {
  arma::cube v_x = arma::zeros(xi.n_elem, y.n_elem, z.n_elem);
  arma::cube v_xy = arma::zeros(xi.n_elem, yi.n_elem, z.n_elem);
  arma::cube v_xyz = arma::zeros(xi.n_elem, yi.n_elem, zi.n_elem);

  arma::vec vec_i = arma::zeros<arma::vec>(x.n_elem);

  for (arma::uword i = 0; i < y.n_elem; ++i) {
    for (arma::uword j = 0; j < z.n_elem; ++j) {
      vec_i = arma::vectorise(v(arma::span(), arma::span(i), arma::span(j)));
      Spline<double> sp(x, vec_i);
      v_x(arma::span(), arma::span(i), arma::span(j)) = sp.interpolate(xi);
    }
  }

  vec_i.set_size(y.n_elem);
  for (arma::uword i = 0; i < xi.n_elem; ++i) {
    for (arma::uword j = 0; j < z.n_elem; ++j) {
      vec_i = vectorise(v_x(arma::span(i), arma::span(), arma::span(j)));
      Spline<double> sp(y, vec_i);
      v_xy(arma::span(i), arma::span(), arma::span(j)) = sp.interpolate(yi);
    }
  }

  vec_i.set_size(z.n_elem);
  for (arma::uword i = 0; i < xi.n_elem; ++i) {
    for (arma::uword j = 0; j < yi.n_elem; ++j) {
      vec_i = vectorise(v_xy(arma::span(i), arma::span(j), arma::span()));
      Spline<double> sp(z, vec_i);
      v_xyz(arma::span(i), arma::span(j), arma::span()) = sp.interpolate(zi);
    }
  }
  return v_xyz;
}

arma::cube interp3(const arma::cube &v, const arma::rowvec &xi,
                   const arma::rowvec &yi, const arma::rowvec &zi) {
  return interp3(arma::regspace<arma::rowvec>(1, v.n_rows),
                 arma::regspace<arma::rowvec>(1, v.n_cols),
                 arma::regspace<arma::rowvec>(1, v.n_slices), v, xi, yi, zi);
}

std::tuple<arma::cube, arma::cube, arma::cube>
ndgrid(const arma::rowvec &v1, const arma::rowvec &v2, const arma::rowvec &v3) {

  arma::cube x(v1.n_elem, v2.n_elem, v3.n_elem, arma::fill::ones), y = x, z = x;
  arma::mat s_mat(v1.n_elem, v2.n_elem, arma::fill::zeros);
  s_mat.each_col() = v1.t();
  x.each_slice() = s_mat;

  s_mat.each_row() = v2;
  y.each_slice() = s_mat;

  for (arma::uword i = 0; i < v3.n_elem; ++i) {
    z.slice(i) *= v3(i);
  }

  return std::make_tuple(x, y, z);
}

std::tuple<arma::cube, arma::cube, arma::cube>
meshgrid(const arma::rowvec &v1, const arma::rowvec &v2,
         const arma::rowvec &v3) {
  arma::cube x2, y2, z2;
  std::tie(y2, x2, z2) = ndgrid(v2, v1, v3);
  return std::make_tuple(x2, y2, z2);
}

arma::cube shift(arma::cube cube_in, arma::rowvec3 shifts) {
  if (cube_in.is_empty()) {
    return {};
  }

  shifts = round(arma::rowvec3(SizeVec(cube_in) % shifts));
  for (arma::uword i = 0; i < 3; ++i) {
    cube_in = shift(cube_in, shifts(i), i);
  }

  return cube_in;
}

arma::vec planar_average(const arma::uword &direction,
                         const arma::cube &cube_in) {
  const arma::uword dim_size = arma::size(cube_in)(direction);
  arma::vec average = arma::vec(dim_size);
  for (arma::uword i = 0; i < dim_size; ++i) {
    switch (direction) {
    case 0:
      average(i) =
          arma::accu(cube_in(arma::span(i), arma::span(), arma::span()));
      break;
    case 1:
      average(i) =
          arma::accu(cube_in(arma::span(), arma::span(i), arma::span()));
      break;
    case 2:
      average(i) =
          arma::accu(cube_in(arma::span(), arma::span(), arma::span(i)));
      break;
    }
  }
  return average;
}

arma::cx_vec fft(arma::vec X) {
  arma::cx_vec out(X.n_elem);
  fftw_plan plan = fftw_plan_dft_r2c_1d(
      X.n_elem, X.memptr(), reinterpret_cast<fftw_complex *>(out.memptr()),
      FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  for (arma::uword i = out.n_elem / 2 + 1; i < out.n_elem; ++i)
    out(i) = conj(out(X.n_rows - i));

  return out;
}

arma::cx_vec fft(arma::cx_vec X) {
  arma::cx_vec out(X.n_elem);
  fftw_plan plan =
      fftw_plan_dft_1d(X.n_elem, reinterpret_cast<fftw_complex *>(X.memptr()),
                       reinterpret_cast<fftw_complex *>(out.memptr()),
                       FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  return out;
}

arma::cx_cube fft(arma::cube X) {
  arma::cx_cube out(X.n_rows / 2 + 1, X.n_cols, X.n_slices);
  fftw_plan plan = fftw_plan_dft_r2c_3d(
      X.n_slices, X.n_cols, X.n_rows, X.memptr(),
      reinterpret_cast<fftw_complex *>(out.memptr()), FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  out.resize(X.n_rows, X.n_cols, X.n_slices);

  for (arma::uword i = X.n_rows / 2 + 1; i < X.n_rows; ++i) {
    out(i, 0, 0) = conj(out(X.n_rows - i, 0, 0));
    for (arma::uword j = 1; j < X.n_cols; ++j) {
      out(i, j, 0) = conj(out(X.n_rows - i, X.n_cols - j, 0));
    }
    for (arma::uword k = 1; k < X.n_slices; ++k) {
      out(i, 0, k) = conj(out(X.n_rows - i, 0, X.n_slices - k));
    }
    for (arma::uword j = 1; j < X.n_cols; ++j) {
      for (arma::uword k = 1; k < X.n_slices; ++k) {
        out(i, j, k) = conj(out(X.n_rows - i, X.n_cols - j, X.n_slices - k));
      }
    }
  }

  return out;
}

arma::cx_cube fft(arma::cx_cube X) {
  arma::cx_cube fft(X.n_rows, X.n_cols, X.n_slices);
  fftw_plan plan =
      fftw_plan_dft_3d(X.n_slices, X.n_cols, X.n_rows,
                       reinterpret_cast<fftw_complex *>(X.memptr()),
                       reinterpret_cast<fftw_complex *>(fft.memptr()),
                       FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  return fft;
}

arma::cx_vec ifft(arma::cx_vec X) {
  arma::cx_vec out(X.n_elem);
  fftw_plan plan =
      fftw_plan_dft_1d(X.n_elem, reinterpret_cast<fftw_complex *>(X.memptr()),
                       reinterpret_cast<fftw_complex *>(out.memptr()),
                       FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  return out / out.n_elem;
}

arma::cx_cube ifft(arma::cx_cube X) {
  arma::cx_cube ifft(X.n_rows, X.n_cols, X.n_slices);
  fftw_plan plan =
      fftw_plan_dft_3d(X.n_slices, X.n_cols, X.n_rows,
                       reinterpret_cast<fftw_complex *>(X.memptr()),
                       reinterpret_cast<fftw_complex *>(ifft.memptr()),
                       FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  return ifft / X.n_elem;
}

arma::SizeCube as_size(const arma::urowvec3 &vec) {
  return arma::SizeCube(vec(0), vec(1), vec(2));
}
arma::SizeMat as_size(const arma::urowvec2 &vec) {
  return arma::SizeMat(vec(0), vec(1));
}

arma::mat fmod(arma::mat mat_in, const double &denom) noexcept {
  mat_in.for_each([&denom](double &val) noexcept { val = fmod(val, denom); });
  return mat_in;
}
arma::mat fmod_p(arma::mat mat_in, const double &denom) noexcept {
  mat_in.for_each([&denom](double &val) noexcept { val = fmod_p(val, denom); });
  return mat_in;
}

double fmod_p(double num, const double &denom) noexcept {
  num = fmod(num, denom);
  if (num < 0)
    num += abs(denom);
  return num;
}

arma::cx_cube poisson_solver_3D(const arma::cx_cube &rho, arma::mat diel,
                                arma::rowvec3 lengths,
                                arma::uword normal_direction) {
  auto n_points = SizeVec(rho);

  if (normal_direction != 2) {
    n_points.swap_cols(normal_direction, 2);
    lengths.swap_cols(normal_direction, 2);
    diel.swap_cols(normal_direction, 2);
  }

  const arma::rowvec Gs = 2.0 * PI / lengths;

  arma::rowvec Gx0 = arma::ceil(arma::regspace<arma::rowvec>(
                         -0.5 * n_points(0), 0.5 * n_points(0) - 1)) *
                     Gs(0);
  arma::rowvec Gy0 = arma::ceil(arma::regspace<arma::rowvec>(
                         -0.5 * n_points(1), 0.5 * n_points(1) - 1)) *
                     Gs(1);
  arma::rowvec Gz0 = arma::ceil(arma::regspace<arma::rowvec>(
                         -0.5 * n_points(2), 0.5 * n_points(2) - 1)) *
                     Gs(2);

  Gx0 = ifftshift(Gx0);
  Gy0 = ifftshift(Gy0);
  Gz0 = ifftshift(Gz0);

  // 4PI is for the atomic units
  const auto rhok = fft(arma::cx_cube(4.0 * PI * rho));
  const arma::cx_mat dielsG = fft(diel);
  const arma::cx_mat eps11 = arma::circ_toeplitz(dielsG.col(0)) / Gz0.n_elem;
  const arma::cx_mat eps22 = arma::circ_toeplitz(dielsG.col(1)) / Gz0.n_elem;
  const arma::cx_mat eps33 = arma::circ_toeplitz(dielsG.col(2)) / Gz0.n_elem;
  const arma::mat GzGzp = Gz0.t() * Gz0;
  const arma::cx_mat Az = eps33 % GzGzp;
  arma::cx_cube Vk(arma::size(rhok));

#pragma omp parallel for firstprivate(Az, eps11, eps22, rhok)
  for (arma::uword k = 0; k < Gx0.n_elem; ++k) {
    const arma::cx_mat eps11_Gx0k2 = eps11 * square(Gx0(k));
    for (arma::uword m = 0; m < Gy0.n_elem; ++m) {
      std::vector<arma::span> spans = {arma::span(k), arma::span(m),
                                       arma::span()};
      std::swap(spans[normal_direction], spans[2]);
      arma::cx_mat AG = Az + eps11_Gx0k2 + eps22 * square(Gy0(m));
      if ((k == 0) && (m == 0)) {
        AG(0, 0) = 1;
      }
      Vk(spans[0], spans[1], spans[2]) =
          arma::solve(AG, arma::vectorise(rhok(spans[0], spans[1], spans[2])));
    }
  }
  // 0,0,0 in k-space corresponds to a constant in the real space: average
  // potential over the supercell.
  Vk(0, 0, 0) = 0;
  const arma::cx_cube V = ifft(Vk);

  return V;
}
