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

enum class model_type :int {
	slab, bulk, monolayer
};

struct slabcc_model {

	bool in_optimization = false;
	// supercell info
	rowvec3 cell_vectors_lengths = { 0, 0, 0 };	// length of the basis vectors (Bohr)
	mat33 cell_vectors = zeros(3, 3);			// supercell vectors (Bohr)
	urowvec3 cell_grid = { 1, 1, 1 };			// grid size
	double voxel_vol = 0;						// voxel volume of the grid in (Bohr^3)
	double cell_volume = 0;						// supercell volume in (Bohr^3)
	
	// slab info
	uword normal_direction = 0;					// index of the normal direction (0/1/2)
	model_type type;
	rowvec2 interfaces = { 0, 0 };
	rowvec diel_in;								// diagonal elements of slab dielectric tensor
	rowvec diel_out;							// diagonal elements of enviroment dielectric tensor
	double diel_erf_beta = 1;
	rowvec3 rounded_relative_shift = {0,0,0};

	// extra charge info
	mat charge_position;			// center of each Gaussian model charge
	mat charge_sigma;				// width of each Gaussian model charges
	mat charge_rotations;			// rotation angles along each axis for the trivariate Gaussians
	rowvec charge_fraction;			// charge fraction in each Gaussian
	double total_charge = 0;		// total charge in CHG
	double defect_charge = 0;		// difference in the charge of the input files
	bool trivariate_charge = false;
	double last_charge_error = 0;		// error in the total charge of the model in the last check

	//calculated data
	double potential_RMSE = 0;
	double initial_potential_RMSE = -1;
	cx_cube CHG; // model charge distribution (e/bohr^3), negative for presence of the electron 

	//potential resulted from the model charge (Hartree)
	cx_cube V;

	//difference of the potential resulted from the model charge (V) and the target potential from QM calculations (V_target)  (eV)
	cube V_diff;

	//original target potential of the extra charge in the input files (eV)
	cube V_target0;
	//reference target of the extra charge with the adjusted grid size (eV)
	cube V_target;
	
	mat dielectric_profiles;

	//sets the cell_vectors, cell_grid, and updates the voxel_vol
	void init_supercell(const mat33& new_vectors, const urowvec3& new_grid);

	// change cell_grid and update the voxel_vol
	void change_grid(const urowvec3& new_cell_grid);

	// change cell_vectors
	// update "interfaces", "charge_position", "cell_vectors_lengths", and "voxel_vol"
	void change_size(const mat33& new_cell_vectors);

	//updates the V_target from the V_target0 to the model grid size
	void update_V_target();

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

	// check for the discretization error and adjust the grid_size
	void adjust_extrapolation_grid(const int& extrapol_steps_num, const double& extrapol_steps_size);
	
	// runs the NLOPT with:
	// algorithm: "opt_algo"
	// tolerance for error: "opt_tol"
	// max number of evaluations: "max_eval"
	// reference to the data: "opt_data"
	// reference to the variables to be optimized: "opt_vars"
	void optimize(const string& opt_algo, const double& opt_tol, const int& max_eval, const int& max_time, const opt_switches& optimize);

	//calculates local: V, V_diff, rhoM (without jellium), diels, Q
	//returns: root mean squared error (RMSE) of the model charge potential 
	double potential_error(const vector<double>& x, vector<double>& grad);

private:
	rowvec Uk(rowvec k) const;
	//updates the voxel_vol from the "cell_vectors_lengths" and "cell_grid"
	void update_voxel_vol();
	//updates the cell_vectors_lengths from the cell_vectors
	void update_cell_vectors_lengths();

	void set_model_type(const bool& model_2d, const rowvec& diel_in, const rowvec& diel_out);
};

//NLOPT wrapper for the slabcc_model::potential_error()
double potential_error(const vector<double>& x, vector<double>& grad, void* slabcc_data);
