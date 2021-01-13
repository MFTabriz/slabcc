// Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "isolated.hpp"
#include "slabcc_consts.hpp"
#include "slabcc_model.hpp"
#include "stdafx.h"
#include "vasp.hpp"

int verbosity_level = 0;

int main(int argc, char *argv[]) {
  slabcc_model model;
  std::string input_file = "slabcc.in";
  std::string output_file = "slabcc.out";
  std::string log_file = "slabcc.log";
  bool output_diffs_only = false;
  cli_params parameters_list = {input_file, output_file, log_file,
                                output_diffs_only};
  parameters_list.parse(argc, argv);
  prepare_output_file(output_file);
  initialize_loggers(log_file, output_file);
  auto log = spdlog::get("loggers");
  auto output_log = spdlog::get("output");

  // default values should be defined in the parser function

  std::string CHGCAR_neutral = "";
  std::string LOCPOT_neutral = "";
  std::string LOCPOT_charged = "";
  std::string CHGCAR_charged = "";
  std::string opt_algo = "";      // optimization algorithm
  arma::mat charge_position;      // center of each Gaussian model charge
  arma::rowvec charge_fraction;   // charge fraction in each Gaussian
  arma::mat charge_sigma;         // width of each Gaussian model charges
  arma::mat charge_rotations;     // rotation angles along each axis for the
                                  // trivariate Gaussians
  bool charge_trivariate = false; // use trivariate Gaussians
  arma::rowvec diel_in;           // diagonal elements of slab dielectric tensor
  arma::rowvec diel_out; // diagonal elements of enviroment dielectric tensor
  arma::rowvec3 slabcenter;
  arma::uword normal_direction = 0; // index of the normal direction (0/1/2)
  arma::rowvec2 interfaces;   // interfaces in relative coordinates, ordered as
                              // the user inputs
  double diel_erf_beta = 0;   // beta value of the erf for dielectric profile
  double opt_tol = 0;         // relative optimization tolerance
  double extrapol_grid_x = 0; // extrapolation grid size multiplier
  double opt_grid_x = 0;      // optimization grid size multiplier
  int max_eval = 0;           // maximum number of steps for the optimization
  int max_time = 0;           // maximum time for the optimization (minutes)
  int extrapol_steps_num = 0; // number of extrapolation steps for E_isolated
  double extrapol_steps_size = 0; // size of each extrapolation step with
                                  // respect to the initial supercell size
  bool optimize = false; // optimizer master switch. Overrides the others if
                         // this one is disabled!
  bool optimize_charge_position = false; // optimize the charge_position
  bool optimize_charge_sigma = false;    // optimize the charge_sigma
  bool optimize_charge_rotation = false; // optimize the charge_rotation
  bool optimize_charge_fraction = false; // optimize the charge_fraction
  bool optimize_interfaces = false;      // optimize the position of interfaces
  bool extrapolate = false; // use the extrapolation for E-isolated calculations
  bool model_2D = false;    // the model is 2D

  // parameters read from the input file
  const input_data inputfile_variables = {CHGCAR_neutral,
                                          LOCPOT_charged,
                                          LOCPOT_neutral,
                                          CHGCAR_charged,
                                          opt_algo,
                                          charge_position,
                                          charge_fraction,
                                          charge_sigma,
                                          charge_rotations,
                                          slabcenter,
                                          diel_in,
                                          diel_out,
                                          normal_direction,
                                          interfaces,
                                          diel_erf_beta,
                                          opt_tol,
                                          optimize,
                                          optimize_charge_position,
                                          optimize_charge_sigma,
                                          optimize_charge_rotation,
                                          optimize_charge_fraction,
                                          optimize_interfaces,
                                          extrapolate,
                                          model_2D,
                                          charge_trivariate,
                                          opt_grid_x,
                                          extrapol_grid_x,
                                          max_eval,
                                          max_time,
                                          extrapol_steps_num,
                                          extrapol_steps_size};

  inputfile_variables.parse(input_file);
  if (!output_diffs_only) {
    inputfile_variables.verify();
  }
  model.set_input_variables(inputfile_variables);

  log->debug("SLABCC: version {}.{}.{}", SLABCC_VERSION_MAJOR,
             SLABCC_VERSION_MINOR, SLABCC_VERSION_PATCH);
  log->debug("Armadillo library: version {}.{}.{}", ARMA_VERSION_MAJOR,
             ARMA_VERSION_MINOR, ARMA_VERSION_PATCH);
  log->debug("NLOPT library: version {}.{}.{}", nlopt::version_major(),
             nlopt::version_minor(), nlopt::version_bugfix());
  log->debug("SPDLOG library: version {}.{}.{}", SPDLOG_VER_MAJOR,
             SPDLOG_VER_MINOR, SPDLOG_VER_PATCH);
  log->debug("SLABCC input file: {}", input_file);
  log->debug("SLABCC output file: {}", output_file);
  log->debug("SLABCC log file: {}", log_file);

  std::vector<std::pair<std::string, std::string>> calculation_results;

  // promises for async read of CHGCAR and POTCAR files
  std::vector<std::future<arma::cube>> future_cells;

  // promises for async file write
  std::vector<std::future<void>> future_files;

  future_cells.push_back(
      std::async(std::launch::async, read_VASP_grid_data, CHGCAR_neutral));
  future_cells.push_back(
      std::async(std::launch::async, read_VASP_grid_data, CHGCAR_charged));
  future_cells.push_back(
      std::async(std::launch::async, read_VASP_grid_data, LOCPOT_neutral));
  future_cells.push_back(
      std::async(std::launch::async, read_VASP_grid_data, LOCPOT_charged));

  supercell Neutral_supercell(CHGCAR_neutral);
  supercell Charged_supercell(CHGCAR_charged);

  Neutral_supercell.charge = future_cells.at(0).get();
  Charged_supercell.charge = future_cells.at(1).get();
  Neutral_supercell.potential = future_cells.at(2).get();
  Charged_supercell.potential = future_cells.at(3).get();

  check_slabcc_compatiblity(Neutral_supercell, Charged_supercell);

  // cell vectors of the CHGCAR and LOCPOT files (bohr)
  const arma::mat33 input_cell_vectors =
      arma::abs(Neutral_supercell.cell_vectors) * Neutral_supercell.scaling *
      ang_to_bohr;
  const arma::urowvec3 input_grid_size = SizeVec(Neutral_supercell.charge);
  model.init_supercell(input_cell_vectors, input_grid_size);

  const arma::rowvec3 relative_shift = 0.5 - slabcenter;
  model.rounded_relative_shift =
      arma::round(model.cell_grid % relative_shift) / model.cell_grid;

  model.interfaces = fmod(
      model.interfaces + model.rounded_relative_shift(normal_direction), 1);

  if (!output_diffs_only) {
    Neutral_supercell.shift(model.rounded_relative_shift);
    Charged_supercell.shift(model.rounded_relative_shift);
    model.charge_position += arma::repmat(model.rounded_relative_shift,
                                          model.charge_position.n_rows, 1);
    model.charge_position = fmod_p(model.charge_position, 1);
  }
  log->debug("Slab normal direction index (0-2): {}", model.normal_direction);
  log->trace("Shift to center done!");

  supercell Defect_supercell = Neutral_supercell;
  Defect_supercell.potential =
      Charged_supercell.potential - Neutral_supercell.potential;
  Defect_supercell.charge = Charged_supercell.charge - Neutral_supercell.charge;

  if (is_active(verbosity::write_defect_file) || output_diffs_only) {
    future_files.push_back(std::async(std::launch::async,
                                      &supercell::write_LOCPOT,
                                      Defect_supercell, "slabcc_D.LOCPOT"));
    future_files.push_back(std::async(std::launch::async,
                                      &supercell::write_CHGCAR,
                                      Defect_supercell, "slabcc_D.CHGCAR"));
  }

  // normalize the charges and potentials
  Neutral_supercell.charge *= -1.0 / model.cell_volume;
  Charged_supercell.charge *= -1.0 / model.cell_volume;
  Defect_supercell.charge *= -1.0 / model.cell_volume;
  Defect_supercell.potential *= -1.0;
  model.POT_target_on_input_grid = Defect_supercell.potential;

  if (output_diffs_only) {
    log->debug("Only the extra charge and the potential difference calculation "
               "have been requested!");
    write_planar_avg(Defect_supercell.potential,
                     Defect_supercell.charge * model.voxel_vol, "D",
                     model.cell_vectors_lengths);
    for (auto &promise : future_files) {
      promise.get();
    }
    finalize_loggers();
    exit(0);
  }

  if (is_active(verbosity::write_planarAvg_file)) {
    write_planar_avg(Neutral_supercell.potential,
                     Neutral_supercell.charge * model.voxel_vol, "N",
                     model.cell_vectors_lengths);
    write_planar_avg(Charged_supercell.potential,
                     Charged_supercell.charge * model.voxel_vol, "C",
                     model.cell_vectors_lengths);
    write_planar_avg(Defect_supercell.potential,
                     Defect_supercell.charge * model.voxel_vol, "D",
                     model.cell_vectors_lengths);
  }

  // total extra charge of the VASP calculation
  model.defect_charge = accu(Defect_supercell.charge) * model.voxel_vol;

  if (std::abs(model.defect_charge) < 0.001) {
    log->debug("Total extra charge: {}", model.defect_charge);
    log->warn("Total extra charge seems to be very small. Please make sure the "
              "path to the input CHGCAR files are "
              "set properly!");
  }
  const opt_switches optimizer_activation_switches{
      optimize_charge_position, optimize_charge_sigma, optimize_charge_rotation,
      optimize_charge_fraction, optimize_interfaces};
  const bool optimize_any = optimize_charge_position || optimize_charge_sigma ||
                            optimize_charge_rotation ||
                            optimize_charge_fraction || optimize_interfaces;

  if (optimize_any) {
    const arma::rowvec2 shifted_interfaces0 = model.interfaces;
    const arma::mat charge_position0 = model.charge_position;
    const arma::urowvec3 cell_grid0 = model.cell_grid;
    const arma::rowvec3 optimization_grid_size =
        opt_grid_x * arma::conv_to<arma::rowvec>::from(model.cell_grid);
    const arma::urowvec3 optimization_grid = {
        static_cast<arma::uword>(optimization_grid_size(0)),
        static_cast<arma::uword>(optimization_grid_size(1)),
        static_cast<arma::uword>(optimization_grid_size(2))};
    model.change_grid(optimization_grid);
    model.update_V_target();
    model.optimize(opt_algo, opt_tol, max_eval, max_time,
                   optimizer_activation_switches);

    // write the unshifted optimized values to the file
    output_log->info("\n[Optimized_model_parameters]");
    if (optimize_interfaces) {
      const arma::rowvec2 optimized_interfaces =
          fmod_p(model.interfaces -
                     model.rounded_relative_shift(model.normal_direction),
                 1);
      output_log->info("interfaces_optimized = {}",
                       to_string(optimized_interfaces));
    }

    if (optimize_charge_fraction) {
      output_log->info("charge_fraction_optimized = {}",
                       to_string(model.charge_fraction));
    }
    if (optimize_charge_sigma) {
      const arma::mat opt_charge_sigma = model.trivariate_charge
                                             ? model.charge_sigma
                                             : model.charge_sigma.col(0);
      output_log->info("charge_sigma_optimized = {}",
                       to_string(opt_charge_sigma));
      model.verify_charge_optimization();
    }
    if (optimize_charge_rotation) {
      const arma::mat rotations = model.charge_rotations * 180.0 / PI;
      output_log->info("charge_rotation_optimized = {}", to_string(rotations));
    }
    if (optimize_charge_position) {
      const arma::mat optimized_charge_position = fmod_p(
          model.charge_position - repmat(model.rounded_relative_shift,
                                         model.charge_position.n_rows, 1),
          1);
      output_log->info("charge_position_optimized = {}",
                       to_string(optimized_charge_position));

      // TODO: need a better algorithm to handle the swaps and PBCs similar to
      // verify_interface_optimization()
      const arma::mat charge_position_change =
          arma::abs(charge_position0 - model.charge_position);
      if (charge_position_change.max() > 0.1) {
        log->warn("The optimized position for the extra charge is "
                  "significantly different from the initial value. "
                  "Please make sure that the final position of the extra "
                  "charge have been estimated correctly!");
        log->debug("Charge position changes: ",
                   to_string(charge_position_change));
      }
    }

    if (optimize_interfaces) {
      model.verify_interface_optimization(shifted_interfaces0);
    }

    if (model.initial_potential_RMSE * (opt_tol + 1) < model.potential_RMSE) {
      // Don't panic! either NLOPT seems to be malfunctioning
      // or we are not correctly logging/checking the result
      log->critical("Optimization failed!");
      log->critical("Potential error of the initial parameters seems to be "
                    "smaller than the optimized parameters! "
                    "You may want to change the initial guess for "
                    "charge_position, change the optimization algorithm, or "
                    "turn off the optimization.");
      log->debug("Initial model potential RMSE: {}",
                 model.initial_potential_RMSE);
      log->debug("Optimized model potential RMSE: {}", model.potential_RMSE);
      finalize_loggers();
      exit(1);
    }
    if (cell_grid0(0) > model.cell_grid(0)) {
      model.change_grid(cell_grid0);
      model.update_V_target();
    }
  }

  log->debug("Cell dimensions (bohr): " +
             to_string(model.cell_vectors_lengths));
  log->debug("Volume (bohr^3): {}", model.cell_volume);

  if (is_active(verbosity::write_defect_file)) {
    supercell Model_supercell = Neutral_supercell;
    // charge is normalized to the VASP CHGCAR convention (rho * Vol)
    // Also, positive value for the electron charge! (the probability of finding
    // an electron)
    Model_supercell.charge =
        -arma::real(model.CHG) * model.voxel_vol * model.CHG.n_elem;
    Model_supercell.potential = -arma::real(model.POT) * Hartree_to_eV;
    future_files.push_back(std::async(std::launch::async,
                                      &supercell::write_CHGCAR, Model_supercell,
                                      "slabcc_M.CHGCAR"));
    future_files.push_back(std::async(std::launch::async,
                                      &supercell::write_LOCPOT, Model_supercell,
                                      "slabcc_M.LOCPOT"));
  }

  if (is_active(verbosity::write_dielectric_file)) {
    model.dielectric_profiles.save("slabcc_DIEL.dat", arma::raw_ascii);
  }
  if (is_active(verbosity::write_planarAvg_file)) {
    write_planar_avg(arma::real(model.POT) * Hartree_to_eV,
                     arma::real(model.CHG) * model.voxel_vol, "M",
                     model.cell_vectors_lengths);
  } else if (is_active(verbosity::write_normal_planarAvg)) {
    write_planar_avg(arma::real(model.POT) * Hartree_to_eV,
                     arma::real(model.CHG) * model.voxel_vol, "M",
                     model.cell_vectors_lengths, model.normal_direction);
  }

  // making sure all the files are written in case of a failure
  for (auto &promise : future_files) {
    promise.get();
  }

  auto local_param = model.data_packer();
  std::vector<double> gradients = {};
  model.potential_RMSE =
      potential_error(std::get<0>(local_param), gradients, &model);
  model.check_V_error();

  model.verify_CHG(Defect_supercell.charge);

  // add jellium to the charge (Because the V is normalized, it is not needed in
  // solving the Poisson eq. but it is needed in the energy calculations)
  model.CHG -= model.total_charge / model.cell_volume;

  const auto farthest_element_index = model.total_charge < 0
                                          ? arma::real(model.POT).index_max()
                                          : arma::real(model.POT).index_min();

  const auto dV = model.POT_diff(farthest_element_index);
  log->info("Potential alignment (dV=): {}", to_string(dV));
  calculation_results.emplace_back("dV", to_string(dV));
  const bool isotropic_screening =
      arma::accu(arma::abs(arma::diff(diel_in))) < 0.02;
  if (std::abs(dV) > 0.05) {
    if (model.type == model_type::bulk && isotropic_screening) {
      log->debug("The potential alignment term (dV) is relatively large. But "
                 "in the isotropic bulk models "
                 "this should not make much difference in the total energy "
                 "correction value!");
    } else {
      log->warn("The potential alignment term (dV) is relatively large. The "
                "constructed model may not be accurate!");
    }
  }

  log->debug("Calculation grid point for the potential alignment term: {}",
             to_string(arma::ind2sub(as_size(model.cell_grid),
                                     farthest_element_index)));

  const double EperModel0 =
      0.5 * arma::accu(arma::real(model.POT) % arma::real(model.CHG)) *
      model.voxel_vol * Hartree_to_eV;
  log->info("E_periodic of the model charge: {}", to_string(EperModel0));
  calculation_results.emplace_back("E_periodic of the model charge",
                                   to_string(EperModel0));

  log->debug("Difference of the charge in the input files: {}",
             to_string(model.defect_charge));
  log->debug("Total charge of the model: {}", to_string(model.total_charge));

  double E_isolated = 0;
  double E_correction = 0;
  if (extrapolate) {

    const arma::rowvec3 extrapolation_grid_size =
        extrapol_grid_x * arma::conv_to<arma::rowvec>::from(model.cell_grid);
    const arma::urowvec3 extrapolation_grid = {
        static_cast<arma::uword>(extrapolation_grid_size(0)),
        static_cast<arma::uword>(extrapolation_grid_size(1)),
        static_cast<arma::uword>(extrapolation_grid_size(2))};
    model.change_grid(extrapolation_grid);
    model.adjust_extrapolation_grid(extrapol_steps_num, extrapol_steps_size);
    if (as_size(model.cell_grid) !=
        as_size(extrapolation_grid)) { // discretization error has been detected
      if (model.type != model_type::monolayer) {
        std::string adjusted_parameters = "";
        if (extrapol_steps_num > 4) {
          extrapol_steps_num = 4;
          adjusted_parameters +=
              " extrapolate_steps_number=" + to_string(extrapol_steps_num);
        }
        if (extrapol_steps_size > 0.25) {
          extrapol_steps_size = 0.25;
          adjusted_parameters +=
              " extrapolate_steps_size=" + to_string(extrapol_steps_size);
        }
        if (adjusted_parameters != "") {
          log->debug("Adjusted parameters:{}", adjusted_parameters);
          model.change_grid(extrapolation_grid);
          model.adjust_extrapolation_grid(extrapol_steps_num,
                                          extrapol_steps_size);
        }
      }
    }

    log->debug("--------------------------------------------------------");
    log->debug(
        "Scaling\tE_periodic\t\tmodel charge\t\tinterfaces\t\tcharge position");
    const arma::rowvec2 interface_pos =
        model.interfaces * model.cell_vectors_lengths(model.normal_direction);
    std::string extrapolation_info =
        to_string(1.0) + "\t" + to_string(EperModel0) + "\t" +
        to_string(model.total_charge) + "\t" + to_string(interface_pos);
    for (arma::uword i = 0; i < model.charge_position.n_rows; ++i) {
      extrapolation_info +=
          "\t" + to_string(model.charge_position(i, model.normal_direction) *
                           model.cell_vectors_lengths(model.normal_direction));
    }
    log->debug(extrapolation_info);
    arma::rowvec Es = arma::zeros<arma::rowvec>(extrapol_steps_num - 1),
                 sizes = Es;
    std::tie(Es, sizes) =
        model.extrapolate(extrapol_steps_num, extrapol_steps_size);

    if (model.type == model_type::monolayer) {
      const arma::rowvec3 unit_cell =
          model.cell_vectors_lengths / arma::max(model.cell_vectors_lengths);
      const auto radius = 10.0;
      const auto ewald_shells = generate_shells(unit_cell, radius);
      const auto madelung_const =
          jellium_madelung_constant(ewald_shells, unit_cell, 1);
      auto madelung_term =
          -std::pow(model.total_charge, 2) * madelung_const / 2;
      nonlinear_fit_data fit_data = {Es, sizes, madelung_term};
      const auto cs = nonlinear_fit(1e-10, fit_data);

      log->info("Madelung constant = " + to_string(madelung_const));
      const std::string fit_params =
          "c0= " + to_string(cs.at(0)) + ", c1=" + to_string(cs.at(1)) +
          ", c2=" + to_string(cs.at(2)) + ", c3=" + to_string(cs.at(3));
      log->info("Non-linear fit parameters:" + fit_params);
      calculation_results.emplace_back("Non-linear fit parameters", fit_params);
      calculation_results.emplace_back("Madelung constant",
                                       to_string(madelung_const));

      E_isolated = cs.at(0) + (cs.at(1) - madelung_term) / cs.at(3);
      E_correction = E_isolated - EperModel0 - model.total_charge * dV;
    } else { // bulk and slab models
      const arma::colvec pols = arma::polyfit(sizes, Es, 1);
      const arma::colvec evals = arma::polyval(pols, sizes.t());
      const auto linearfit_MSE =
          arma::accu(arma::square(evals.t() - Es)) / Es.n_elem * 100;
      const arma::rowvec slopes = arma::diff(Es) / arma::diff(sizes);
      const auto extrapol_error_periodic =
          std::abs(slopes(0) - slopes(slopes.n_elem - 1));
      log->debug("--------------------------------------------------------");
      log->debug("Linear fit: Eper(Model) = {}/scaling + {}",
                 to_string(pols(0)), to_string(pols(1)));
      log->debug("Linear fit Root Mean Square Error: {}",
                 to_string(sqrt(linearfit_MSE)));
      log->debug("Polyfit evaluated energies: {}", to_string(evals));
      log->debug("Linear fit error for the periodic model: {}",
                 to_string(extrapol_error_periodic));

      if (extrapol_error_periodic > 0.05) {
        log->debug("Extrapolation energy slopes: {}", to_string(slopes));
        log->critical(
            "The extrapolated energies are not scaling linearly as expected!");
        if (model.type != model_type::bulk) {
          log->critical("The slab thickness may be too small for this "
                        "extrapolation algorithm. "
                        "For calculating the charge correction energy for the "
                        "2D models use \"2D_model = "
                        "yes\" in the input file.");
        }
        finalize_loggers();
        exit(1);
      }

      E_isolated = EperModel0 - pols(0);
      E_correction = -pols(0) - model.total_charge * dV;
    }
    log->info("E_isolated from extrapolation with {}x{} steps: {}",
              to_string(extrapol_steps_num), to_string(extrapol_steps_size),
              to_string(E_isolated));
  } else {
    if (model.type == model_type::monolayer) {
      E_isolated = model.Eiso_bessel();
      E_correction = E_isolated - EperModel0 - model.total_charge * dV;
      log->info(
          "E_isolated from the Bessel expansion of the Poisson equation: {}",
          to_string(E_isolated));
    } else {
      // input parameter checking function must prevent this from happening!
      log->critical("There is no algorithm other than the extrapolation for "
                    "E_isolated calculation of the slab "
                    "models in this version of the slabcc!");
      finalize_loggers();
      exit(1);
    }
  }
  calculation_results.emplace_back("E_isolated of the model charge",
                                   to_string(E_isolated));

  log->info("Energy correction for the model charge (E_iso-E_per-q*dV=): {}",
            to_string(E_correction));
  calculation_results.emplace_back(
      "Energy correction for the model charge (E_iso-E_per-q*dV)",
      to_string(E_correction));
  log->flush();

  finalize_loggers();
  output_log->info("\n[Results]");
  for (const auto &i : calculation_results) {
    output_log->info("{} = {}", i.first, i.second);
  }
  output_log->flush();

  // making sure all the files are written
  for (auto &promise : future_files) {
    promise.get();
  }

  log->trace("Calculations successfully ended!");
  return 0;
}