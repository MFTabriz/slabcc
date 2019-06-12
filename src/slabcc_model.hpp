#pragma once
#include "slabcc_math.hpp"
#include "slabcc_input.hpp"

extern const double Hartree_to_eV;
extern const double ang_to_bohr;

struct opt_variable {
	rowvec2& interfaces;
	mat& charge_sigma, & charge_rotations;
	rowvec& charge_fraction;
	mat& charge_position;
};


struct opt_switches {
	const bool& charge_position, & charge_sigma, & charge_rotation, & charge_fraction, & interfaces;
};

struct slabcc_model {
	// supercell info
	rowvec3 cell_vectors_lengths = { 0, 0, 0 };	// length of the basis vectors in Bohr
	mat33 cell_vectors = zeros(3, 3);		// supercell vectors in Bohr
	urowvec3 cell_grid = { 1, 1, 1 };		// grid size
	double voxel_vol = 0;				// voxel volume of the grid in Bohr^3
	double cell_volume = 0;			//cell volume in bohr ^ 3
	
	// slab info
	uword normal_direction = 0;
	bool model_2d = false;
	rowvec2 interfaces = { 0, 0 };
	rowvec diel_in;
	rowvec diel_out;
	double diel_erf_beta = 1;
	rowvec3 rounded_relative_shift = {0,0,0};

	// extra charge info
	mat charge_position;
	mat charge_sigma;
	mat charge_rotations;
	rowvec charge_fraction;
	double total_charge = 0;
	double defect_charge = 0;
	bool trivariate_charge = false;

	//calculated data
	double potential_RMSE = 0;
	double initial_potential_RMSE = -1;
	cx_cube CHG; // model charge distribution (e/bohr^3), negative for presence of the electron 

	//potential resulted from the model charge (Hartree)
	cx_cube V;

	//difference of the potential resulted from the model charge and the QM calculation (VASP) results (eV)
	cube V_diff;

	//reference potential of the extra charge in the input files (eV) from the input files
	cube V_ref0;
	//reference potential of the extra charge with the adjusted grid size (eV)
	cube V_ref;
	
	mat dielectric_profiles;

	//Sets the global struct slabcc_cell parameters from the new cell vectors (in Bohr) and the grid density
	void init_supercell(const mat33& new_vectors, const urowvec3& new_grid);

	//updates the V_ref from the V_ref0 to the model grid size
	void update_Vref();

	//Sets the global struct slabcc_cell parameters from the new cell vectors (in Bohr) and the grid density
	void update_supercell(const mat33& new_vectors, const urowvec3& new_grid);

	// must be checked before!
	void set_input_variables(const input_data& inputfile_variables);

	// returns the extrapolated sizes and the energies
	tuple <rowvec, rowvec> extrapolate(int extrapol_steps_num, double extrapol_steps_size);

	// generates dielectric profile matrix with each column representing the 
	// dielectric tensor elements' variation in the normal direction.
	void dielectric_profiles_gen();

	// produces Gaussian charge distribution in real space
	// the generated charge distribution data is in (e/bohr^3)
	void gaussian_charges_gen();


	//pack the optimization variable structure and their lower and upper boundaries into std::vector<double> for NLOPT
	//returned vectors are "optimization parameters", "lower boundaries", "upper boundaries"
	//TODO: to keep the compatibility with MSVC, the default values must be replaced with std::optional in the C++17 standard version of slabcc
	tuple<vector<double>, vector<double>, vector<double>, vector<double>> data_packer(opt_switches optimize = { false,false,false,false,false }) const;


	//writes the parameters from the optimization variable structure into the model
	void data_unpacker(const vector<double>& optimizer_vars_vec);

	//calculate the changes in the interface positions and warn the user if the changes are too big
	void verify_interface_optimization(const rowvec2& initial_interfaces) const;

	// check for large charge_sigma and report the delocalized charges
	void verify_charge_optimization() const;


	// calculate the isolated energy from the Bessel expansion of the Poisson equation
	double Eiso_bessel() const;

	// increases the grid size if there is huge discretization error in the model charge
	bool had_discretization_error();

	//check for the discretization error and adjust the extrapol_grid_x
	//for the 2d models extrapol_steps_num and extrapol_steps_size may also be adjusted
	void adjust_extrapolation_params(int& extrapol_steps_num, double& extrapol_steps_size, double& extrapol_grid_x);

private:
	rowvec Uk(rowvec k) const;

};



//input data for the optimizer function
struct opt_data {
	slabcc_model& model;
};

// runs the NLOPT with:
// algorithm: "opt_algo"
// tolerance for error: "opt_tol"
// max number of evaluations: "max_eval"
// reference to the data: "opt_data"
// reference to the variables to be optimized: "opt_vars"
void optimize(const string& opt_algo, const double& opt_tol, const int& max_eval, const int& max_time, opt_data& opt_data, const opt_switches& optimize);


//calculates local: V, V_diff, rhoM (without jellium), diels, Q
//returns: root mean squared error (RMSE) of the model charge potential 
double potential_error(const vector<double>& x, vector<double>& grad, void* slabcc_data);
