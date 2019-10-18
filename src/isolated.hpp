#pragma once
#include "arma_io.hpp"
#include "madelung.hpp"
#include "slabcc_math.hpp"

struct nonlinear_fit_data {
  arma::rowvec &energies, &sizes;
  double &madelung_term;
};

// evaluates the MSE for fitting to the 2nd-order + exponential function as
// described in the Erratum of the paper
double fit_eval(const std::vector<double> &x, std::vector<double> &grad,
                void *data);

// fit the extrapolated energies to customized (2nd-order + exponential)
// function form (needed for non-linear energies of the extrapolate_2D)
std::vector<double> nonlinear_fit(const double &opt_tol,
                                  nonlinear_fit_data &fit_data);
