#pragma once
#include "slabcc_consts.hpp"
#include "slabcc_input.hpp"
#include "slabcc_math.hpp"
#include "vasp.hpp"

struct opt_variable {
  arma::rowvec2 &interfaces;
  arma::mat &charge_sigma, &charge_rotations;
  arma::rowvec &charge_fraction;
  arma::mat &charge_position;
};

struct opt_switches {
  const bool &charge_position, &charge_sigma, &charge_rotation,
      &charge_fraction, &interfaces;
};

enum class model_type : int { slab, bulk, monolayer };

struct slabcc_model {

  bool in_optimization = false;
  // supercell info
  arma::rowvec3 cell_vectors_lengths = {
      0, 0, 0}; // length of the basis vectors (Bohr)
  arma::mat33 cell_vectors = arma::zeros(3, 3); // supercell vectors (Bohr)
  arma::urowvec3 cell_grid = {1, 1, 1};         // grid size
  double voxel_vol = 0;   // voxel volume of the grid (Bohr^3)
  double cell_volume = 0; // supercell volume (Bohr^3)

  // slab info
  arma::uword normal_direction = 0; // index of the normal direction (0/1/2)
  model_type type;
  arma::rowvec2 interfaces = {0, 0};
  arma::rowvec diel_in;  // diagonal elements of slab dielectric tensor
  arma::rowvec diel_out; // diagonal elements of enviroment dielectric tensor
  double diel_erf_beta = 1;
  arma::rowvec3 rounded_relative_shift = {0, 0, 0};

  // extra charge info
  arma::mat charge_position;    // center of each Gaussian model charge
  arma::mat charge_sigma;       // width of each Gaussian model charges
  arma::mat charge_rotations;   // rotation angles along each axis for the
                                // trivariate Gaussians
  arma::rowvec charge_fraction; // charge fraction in each Gaussian
  double total_charge = 0;      // total charge in CHG
  double defect_charge = 0;     // difference in the charge of the input files
  bool trivariate_charge = false;
  double last_charge_error =
      0; // error in the total charge of the model in the last check

  // calculated data
  double potential_RMSE = 0;
  double initial_potential_RMSE = -1;
  // model charge distribution (e/bohr^3), negative for presence of the electron
  arma::cx_cube CHG;

  // potential resulted from the model charge (Hartree)
  arma::cx_cube POT;

  // difference of the potential resulted from the model charge (POT) and the
  // target potential from QM calculations (POT_target) (eV)
  arma::cube POT_diff;

  // reference target of the extra charge with the adjusted grid size (eV)
  arma::cube POT_target;

  // original target potential of the extra charge in the input files (eV)
  arma::cube POT_target_on_input_grid;

  arma::mat dielectric_profiles;

  // sets the cell_vectors, cell_grid, and updates the voxel_vol
  void init_supercell(const arma::mat33 &new_vectors,
                      const arma::urowvec3 &new_grid);

  // change cell_grid and update the voxel_vol
  void change_grid(const arma::urowvec3 &new_cell_grid);

  // change cell_vectors, update "interfaces", "charge_position",
  // "cell_vectors_lengths", and "voxel_vol"
  void change_size(const arma::mat33 &new_cell_vectors);

  // updates the POT_target from the POT_target_on_input_grid to the model grid
  // size
  void update_V_target();

  // must be checked before!
  void set_input_variables(const input_data &inputfile_variables);

  // returns the extrapolated sizes and the energies
  std::tuple<arma::rowvec, arma::rowvec>
  extrapolate(int extrapol_steps_num, double extrapol_steps_size);

  // generates dielectric profile matrix with each column representing the
  // dielectric tensor elements' variation in the normal direction.
  void dielectric_profiles_gen();

  // produces Gaussian charge distribution in real space
  // the generated charge distribution data is in (e/bohr^3)
  void gaussian_charges_gen();

  // pack the optimization variable structure and their lower and upper
  // boundaries into std::vector<double> for NLOPT returned vectors are
  // "optimization parameters", "lower boundaries", "upper boundaries"
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
             std::vector<double>>
  data_packer(opt_switches optimize = opt_switches{false, false, false, false,
                                                   false}) const;

  // writes the parameters from the optimization variable structure into the
  // model
  void data_unpacker(const std::vector<double> &optimizer_vars_vec);

  // checks the amount of the charge between the interfaces for the model and
  // compares it to the input files
  void verify_CHG(const arma::cube &defect_charge);

  // calculate the changes in the interface positions and warn the user if the
  // changes are too big
  void
  verify_interface_optimization(const arma::rowvec2 &initial_interfaces) const;

  // check for large charge_sigma and report the delocalized charges
  void verify_charge_optimization() const;

  // calculate the isolated energy from the Bessel expansion of the Poisson
  // equation
  double Eiso_bessel() const;

  // increases the grid size if there is huge discretization error in the model
  // charge
  bool had_discretization_error();

  // check for the discretization error and adjust the grid_size
  void adjust_extrapolation_grid(const int &extrapol_steps_num,
                                 const double &extrapol_steps_size);

  // runs the NLOPT with:
  // algorithm: "opt_algo"
  // tolerance for error: "opt_tol"
  // max number of evaluations: "max_eval"
  // reference to the data: "opt_data"
  // reference to the variables to be optimized: "opt_vars"
  void optimize(const std::string &opt_algo, const double &opt_tol,
                const int &max_eval, const int &max_time,
                const opt_switches &optimize);

  // calculates local: POT, POT_diff, rhoM (without jellium), diels, Q
  // returns: root mean squared error (RMSE) of the model charge potential
  double potential_error(const std::vector<double> &x,
                         std::vector<double> &grad);

  // checks the potential_RMSE and its directional values
  void check_V_error();

private:
  arma::rowvec Uk(arma::rowvec k) const;
  // updates the voxel_vol from the "cell_vectors_lengths" and "cell_grid"
  void update_voxel_vol();
  // updates the cell_vectors_lengths from the cell_vectors
  void update_cell_vectors_lengths();

  void set_model_type(const bool &model_2d, const arma::rowvec &diel_in,
                      const arma::rowvec &diel_out);
};

// NLOPT wrapper for the slabcc_model::potential_error()
double potential_error(const std::vector<double> &x, std::vector<double> &grad,
                       void *slabcc_data);
