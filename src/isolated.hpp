#pragma once
#include "arma_io.hpp"
#include "slabcc_math.hpp"
#include "madelung.hpp"

using namespace arma;
using namespace std;

struct nonlinear_fit_data {
	rowvec &energies, &sizes;
	double &madelung_term;
};


// evaluates the MSE for fitting to the 2nd-order + exponential function as described in the Erratum of the paper
double fit_eval(const vector<double> &x, vector<double> &grad, void *data);

// fit the extrapolated energies to customized (2nd-order + exponential) function form (needed for non-linear energies of the extrapolate_2D)
vector<double> nonlinear_fit(const double& opt_tol, nonlinear_fit_data& fit_data);


