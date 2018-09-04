#pragma once
#include "vasp.hpp"
#include "general_io.hpp"
#include "nlopt.hpp"
#include "INIReader.h"
#include "madelung.hpp"

extern const double Hartree_to_eV;

// These functions are providing the basic functions needed in various stages of calculation:
// generating a 1D dielectric profile
// generating a 3D Gaussian charge distribution
// solving the Poisson equation
// estimating the average error of the model potential
// optimizing the model parameters

struct slabcc_cell_type {

	// lengths of the basis vectors in Bohr
	rowvec3 vec_lengths = { 0, 0, 0 };

	// size of the supercell in Bohr
	mat33 size = zeros(3, 3);

	//calculation grid size
	urowvec3 grid = { 1, 1, 1 };

	uword normal_direction = 0;

	// voxel volume of the grid in Bohr^3
	double voxel_vol = 0;
};

//input data for checker
struct input_data {
	string &CHGCAR_NEU, &LOCPOT_CHG, &LOCPOT_NEU, &CHGCAR_CHG, &opt_algo;
	mat &charge_position;
	rowvec &Qd, &sigma;
	rowvec3 &slabcenter;
	rowvec &diel_in, &diel_out;
	uword &normal_direction;
	rowvec2 &interfaces;
	double &diel_erf_beta, &opt_tol;
	bool &optimize_charge, &optimize_interface, &extrapol_slab;
	double &opt_grid_x, &extrapol_grid_x;
	int &max_eval, &max_time, &extrapol_steps_num;
	double &extrapol_steps_size;
};

//input data for the optimizer function
struct opt_data {
	const double &Q0, &diel_erf_beta;
	const rowvec &diel_in, &diel_out;
	const cube &defect_potential;
	mat &diels;
	cx_cube& rhoM;
	cx_cube& V;
	cube& V_diff;
	double &initial_pot_MSE;
};

struct opt_vars {
	rowvec2 &interfaces;
	rowvec &sigma, &Qd;
	mat &charge_position;
};

struct nonlinear_fit_data {
	rowvec &energies, &sizes;
	double &madelung_term;
};
// generates dielectric profile matrix with each column representing the 
// dielectric tensor elements' variation in the normal direction.
mat dielectric_profiles(const rowvec2 &interfaces, const rowvec3 &diel_in, const rowvec3 &diel_out, const double &diel_erf_beta);

//Sets the global struct slabcc_cell parameters from cell vectors "size" (in Bohr) and grid density "grid"
void UpdateCell(const mat33& size, const urowvec3& grid);

//Produces Gaussian charge distribution in real space
// Q is total charge, rel_pos is the relative position of the center of Gaussian charge,
// sigma is the Gaussian width
// the generated charge distribution data is in (e/bohr^3)
cx_cube gaussian_charge(const double& Q, const vec3& rel_pos, const double& sigma);


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
double do_optimize(const string& opt_algo, const double& opt_tol, const int &max_eval, const int &max_time, opt_data& opt_data, opt_vars& opt_vars, const bool &optimize_charge, const bool &optimize_interfaces);

//pack the optimization variable structure and their lower and upper boundaries into std::vector<double> for NLOPT
//returned vectors are "optimization parameters", "lower boundaries", "upper boundaries"
tuple<vector<double>, vector<double>, vector<double>> optimizer_packer(const opt_vars& opt_vars, const bool optimize_charge = false, const bool optimize_interface = false);

//unpack the optimization variable structure into opt_vars struct variables
void optimizer_unpacker(const vector<double> &optimizer_vars_vec, opt_vars &opt_vars);

//sanity check on the input parameters
void check_inputs(input_data input_set);

//check or enforce some conditions on the supercell grid sizes and the shape
void check_cells(supercell& Neutral_supercell, supercell& Charged_supercell, input_data input_set);

//parse the parameters from the input file
void parse_input_params(const string& input_file, ofstream& output_fstream, const input_data& input_set);

//optimization constraint to ensure all the Gaussian charges have the same sign
double opt_charge_constraint(const vector<double> &x, vector<double> &grad, void *data);

// returns the extrapolated sizes and the energies
tuple <rowvec, rowvec> extrapolate_3D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &Qd, const rowvec &sigma, const double &grid_multiplier);
tuple <rowvec, rowvec> extrapolate_2D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &Qd, const rowvec &sigma, const double &grid_multiplier);

// evaluates the MSE for fitting to the 2nd-order + exponential function as described in the Erratum of the paper
double fit_eval(const vector<double> &x, vector<double> &grad, void *data);

// fit the extrapolated energies to customized (2nd-order + exponential) function form (needed for non-linear energies of the extrapolate_2D)
vector<double> nonlinear_fit(const double& opt_tol, nonlinear_fit_data& fit_data);