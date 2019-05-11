#pragma once
#include "arma_io.hpp"
#include "slabcc_math.hpp"
#include "slabcc_core.hpp"
#include "madelung.hpp"

using namespace arma;
using namespace std;
extern slabcc_cell model_cell;

struct nonlinear_fit_data {
	rowvec &energies, &sizes;
	double &madelung_term;
};

// returns the extrapolated sizes and the energies
tuple <rowvec, rowvec> extrapolate_3D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &charge_q, const mat &charge_sigma, const mat &charge_rotations, const double &grid_multiplier, const bool &trivariate);
tuple <rowvec, rowvec> extrapolate_2D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &charge_q, const mat &charge_sigma, const mat &charge_rotations, const double &grid_multiplier, const bool &trivariate);

// evaluates the MSE for fitting to the 2nd-order + exponential function as described in the Erratum of the paper
double fit_eval(const vector<double> &x, vector<double> &grad, void *data);

// fit the extrapolated energies to customized (2nd-order + exponential) function form (needed for non-linear energies of the extrapolate_2D)
vector<double> nonlinear_fit(const double& opt_tol, nonlinear_fit_data& fit_data);

// calculate the isolated energy from the Bessel expansion of the Poisson equation
double Eiso_bessel(double Q, double z0, double sigma, mat diel);
rowvec Uk(rowvec k, double z0, double sigma, mat diel);

