#pragma once
#include "slabcc_cell.hpp"
#include "slabcc_input.hpp"
#include "vasp.hpp"
#include "general_io.hpp"
#include "nlopt.hpp"

extern const double Hartree_to_eV;

// These functions are providing the basic functions needed in various stages of calculation:
// generating the 1D dielectric profiles
// generating a 3D Gaussian charge distribution
// solving the Poisson equation
// estimating the average error of the model potential
// optimizing the model parameters



//input data for the optimizer function
struct opt_data {
	const double &total_vasp_charge, &diel_erf_beta;
	const rowvec &diel_in, &diel_out;
	const cube &defect_potential;
	const bool &trivariate;
	const rowvec3 &rounded_relative_shift;
	mat &dielectric_profiles;
	cx_cube& rhoM;
	cx_cube& V;
	cube& V_diff;
	double &initial_potential_RMSE;
};

struct opt_switch {
	const bool &charge_position, &charge_sigma, &charge_rotation, &charge_fraction, &interfaces;
};

struct opt_variable {
	rowvec2 &interfaces;
	mat &charge_sigma, &charge_rotations;
	rowvec &charge_fraction;
	mat &charge_position;
};

// generates dielectric profile matrix with each column representing the 
// dielectric tensor elements' variation in the normal direction.
mat dielectric_profiles_gen(const rowvec2 &interfaces, const rowvec3 &diel_in, const rowvec3 &diel_out, const double &diel_erf_beta);


//Produces Gaussian charge distribution in real space
// Q is total charge, rel_pos is the relative position of the center of Gaussian charge,
// sigma is the Gaussian width in x/y/z direction (for simple Gaussians only the 1st element is used)
// the generated charge distribution data is in (e/bohr^3)
cx_cube gaussian_charge(const double& Q, const rowvec3& rel_pos, const rowvec3& sigma, const rowvec3& rotation_angle, const bool& trivariate);


//Poisson solver in 3D with anisotropic dielectric profiles
//diel is the N*3 matrix of variations in dielectric tensor elements in direction normal to the surface
cx_cube poisson_solver_3D(const cx_cube &rho, mat diel);


//calculates local: V, V_diff, rhoM (without jellium), diels, Q
//returns: mean squared error (MSE) of the model charge potential 
double potential_eval(const vector<double> &x, vector<double> &grad, void *slabcc_data);

// runs the NLOPT with:
// algorithm: "opt_algo"
// tolerance for error: "opt_tol"
// max number of evaluations: "max_eval"
// reference to the data: "opt_data"
// reference to the variables to be optimized: "opt_vars"
double do_optimize(const string& opt_algo, const double& opt_tol, const int &max_eval, const int &max_time, opt_data& opt_data, opt_variable& opt_vars, opt_switch& optimize);

//pack the optimization variable structure and their lower and upper boundaries into std::vector<double> for NLOPT
//returned vectors are "optimization parameters", "lower boundaries", "upper boundaries"
tuple<vector<double>, vector<double>, vector<double>, vector<double>> optimizer_packer(const opt_variable& opt_vars, opt_switch optimize = { false,false,false,false,false }, const bool trivariate = false);

//unpack the optimization variable structure into opt_vars struct variables
void optimizer_unpacker(const vector<double> &optimizer_vars_vec, opt_variable &opt_vars);

//check conditions and consistency of the supercell grid sizes and the shape
void verify_cells(const supercell& Neutral_supercell, const supercell& Charged_supercell);

//calculate the changes in the interface positions and warn the user if the changes are too big
void verify_interface_optimization(const rowvec2& initial_interfaces, const rowvec2& optimized_interfaces);

//check for large charge_sigma and report the delocalized charges
void verify_charge_sigma_optimization(const rowvec& charge_q, const mat& charge_sigma);

//optimization constraint to ensure all the Gaussian charges have the same sign
double opt_charge_constraint(const vector<double> &x, vector<double> &grad, void *data);
