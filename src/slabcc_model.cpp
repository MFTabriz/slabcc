// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_model.hpp"

void slabcc_model::set_input_variables(const input_data &inputfile_variables) {
  normal_direction = inputfile_variables.normal_direction;
  interfaces = inputfile_variables.interfaces;
  diel_in = inputfile_variables.diel_in;
  diel_out = inputfile_variables.diel_out;
  diel_erf_beta = inputfile_variables.diel_erf_beta;
  charge_position = inputfile_variables.charge_position;
  charge_sigma = inputfile_variables.charge_sigma;
  charge_rotations = inputfile_variables.charge_rotations;
  charge_fraction = inputfile_variables.charge_fraction;
  trivariate_charge = inputfile_variables.trivariate;
  set_model_type(inputfile_variables.model_2D, diel_in, diel_out);
};

void slabcc_model::set_model_type(const bool &model_2d,
                                  const arma::rowvec &diel_in,
                                  const arma::rowvec &diel_out) {
  auto log = spdlog::get("loggers");
  if (arma::approx_equal(diel_in, diel_out, "absdiff", 0.02)) {
    this->type = model_type::bulk;
    log->debug("Model type: Bulk");
  } else if (model_2d) {
    this->type = model_type::monolayer;
    log->debug("Model type: 2D");
  } else {
    this->type = model_type::slab;
    log->debug("Model type: Slab");
  }
}

void slabcc_model::init_supercell(const arma::mat33 &new_vectors,
                                  const arma::urowvec3 &new_grid) {
  cell_vectors = new_vectors;
  cell_grid = new_grid;
  update_cell_vectors_lengths();
  update_voxel_vol();
}

void slabcc_model::change_grid(const arma::urowvec3 &new_cell_grid) {
  cell_grid = new_cell_grid;
  update_voxel_vol();
}

void slabcc_model::change_size(const arma::mat33 &new_cell_vectors) {
  auto log = spdlog::get("loggers");
  const arma::rowvec3 cell_vectors_lengths0 = cell_vectors_lengths;
  cell_vectors = new_cell_vectors;
  update_cell_vectors_lengths();
  update_voxel_vol();

  interfaces *= cell_vectors_lengths0(normal_direction) /
                cell_vectors_lengths(normal_direction);
  charge_position *= cell_vectors_lengths0(normal_direction) /
                     cell_vectors_lengths(normal_direction);

  const arma::rowvec3 scaling = cell_vectors_lengths / cell_vectors_lengths0;
  if (abs(max(scaling) - min(scaling)) > 0.00001) {
    log->warn("Model cell scaling is anisotropic!");
  }
}

void slabcc_model::update_voxel_vol() {
  voxel_vol = prod(cell_vectors_lengths / cell_grid);
}

void slabcc_model::update_cell_vectors_lengths() {
  for (arma::uword i = 0; i < 3; ++i) {
    cell_vectors_lengths(i) = norm(cell_vectors.col(i));
  }

  //(only works in the orthogonal case!)
  cell_volume = arma::prod(cell_vectors_lengths);
}

void slabcc_model::dielectric_profiles_gen() {
  const auto length = cell_vectors_lengths(normal_direction);
  const auto n_points = cell_grid(normal_direction);
  arma::rowvec2 interfaces_cartesian = interfaces * length;
  interfaces_cartesian = arma::sort(interfaces_cartesian);
  const auto positions = arma::linspace<arma::rowvec>(0, length, n_points + 1);
  dielectric_profiles = arma::zeros<arma::mat>(n_points, 3);
  const arma::rowvec3 diel_sum = diel_in + diel_out;
  const arma::rowvec3 diel_diff = diel_out - diel_in;

  for (arma::uword k = 0; k < n_points; ++k) {
    double min_distance = 0;

    double diel_side = 1;
    // min distance of the point positions(k) to the interfaces in PBC
    // [-length/2  length/2]
    const arma::rowvec2 distances =
        fmod_p(positions(k) - interfaces_cartesian + length / 2, length) -
        length / 2;
    if (abs(distances(0)) < abs(distances(1))) {
      min_distance = distances(0);
      diel_side = -1;
    } else {
      min_distance = distances(1);
    }

    const auto diel_edge = erf(min_distance / diel_erf_beta);
    dielectric_profiles.row(k) =
        (diel_diff * diel_side * diel_edge + diel_sum) / 2;
  }
}

void slabcc_model::gaussian_charges_gen() {

  do {
    arma::rowvec x0 = arma::linspace<arma::rowvec>(
        0, cell_vectors_lengths(0) - cell_vectors_lengths(0) / cell_grid(0),
        cell_grid(0));
    arma::rowvec y0 = arma::linspace<arma::rowvec>(
        0, cell_vectors_lengths(1) - cell_vectors_lengths(1) / cell_grid(1),
        cell_grid(1));
    arma::rowvec z0 = arma::linspace<arma::rowvec>(
        0, cell_vectors_lengths(2) - cell_vectors_lengths(2) / cell_grid(2),
        cell_grid(2));

    for (arma::uword i = 0; i < charge_fraction.n_elem; ++i) {
      // shift the axis reference to position of the Gaussian charge center
      arma::rowvec x =
          x0 - arma::accu(cell_vectors.col(0) * charge_position(i, 0));
      arma::rowvec y =
          y0 - arma::accu(cell_vectors.col(1) * charge_position(i, 1));
      arma::rowvec z =
          z0 - arma::accu(cell_vectors.col(2) * charge_position(i, 2));
      // handle the minimum distance from the mirror charges
      for (auto &pos : x) {
        if (abs(pos) > cell_vectors_lengths(0) / 2) {
          pos = cell_vectors_lengths(0) - abs(pos);
        }
      }
      for (auto &pos : y) {
        if (abs(pos) > cell_vectors_lengths(1) / 2) {
          pos = cell_vectors_lengths(1) - abs(pos);
        }
      }
      for (auto &pos : z) {
        if (abs(pos) > cell_vectors_lengths(2) / 2) {
          pos = cell_vectors_lengths(2) - abs(pos);
        }
      }

      arma::cube xs, ys, zs;
      std::tie(xs, ys, zs) = ndgrid(x, y, z);

      const arma::rowvec3 rotation_angle = charge_rotations.row(i);
      // rotate around xyz axis
      if (arma::max(arma::abs(rotation_angle)) > 0.002) {
        const arma::mat33 rot_x = {
            {1, 0, 0},
            {0, cos(rotation_angle(0)), -sin(rotation_angle(0))},
            {0, sin(rotation_angle(0)), cos(rotation_angle(0))}};

        const arma::mat33 rot_y = {
            {cos(rotation_angle(1)), 0, sin(rotation_angle(1))},
            {0, 1, 0},
            {-sin(rotation_angle(1)), 0, cos(rotation_angle(1))}};

        const arma::mat33 rot_z = {
            {cos(rotation_angle(2)), -sin(rotation_angle(2)), 0},
            {sin(rotation_angle(2)), cos(rotation_angle(2)), 0},
            {0, 0, 1}};

        const arma::mat33 rotation_mat = rot_x * rot_y * rot_z;
        for (arma::uword i = 0; i < xs.n_elem; ++i) {
          const arma::vec3 old_coordinates = {xs(i), ys(i), zs(i)};
          const arma::vec3 new_coordinates = rotation_mat * old_coordinates;
          xs(i) = new_coordinates(0);
          ys(i) = new_coordinates(1);
          zs(i) = new_coordinates(2);
        }
      }

      const arma::cube r2 =
          arma::square(xs) + arma::square(ys) + arma::square(zs);

      // this charge distribution is due to the 1st nearest gaussian image.
      // In case of the very small supercells or very diffuse charges (large
      // sigma), the higher order of the image charges must also be included.
      // But the validity of the correction method for these cases must be
      // checked!

      const double Q = charge_fraction(i) * defect_charge;
      CHG = arma::zeros<arma::cx_cube>(as_size(cell_grid));

      if (trivariate_charge) {
        CHG += arma::cx_cube(
            Q / (pow(2 * PI, 1.5) * arma::prod(charge_sigma.row(i))) *
                exp(-arma::square(xs) / (2 * square(charge_sigma(i, 0))) -
                    arma::square(ys) / (2 * square(charge_sigma(i, 1))) -
                    arma::square(zs) / (2 * square(charge_sigma(i, 2)))),
            arma::zeros(as_size(cell_grid)));
      } else {
        CHG += arma::cx_cube(Q / pow((charge_sigma(i, 0) * sqrt(2 * PI)), 3) *
                                 exp(-r2 / (2 * square(charge_sigma(i, 0)))),
                             arma::zeros(as_size(cell_grid)));
      }
    }

    total_charge = arma::accu(arma::real(CHG)) * voxel_vol;
  } while (had_discretization_error());
  update_V_target();
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
           std::vector<double>>
slabcc_model::data_packer(opt_switches optimize) const {
  auto log = spdlog::get("loggers");
  // size of the first step for each parameter
  const double move_step = 1; // Ang
  const double sigma_step = 0.2;
  const double fraction_step = 0.2;
  const double rotation_step = 0.2; // rad
  const arma::rowvec3 relative_move_step =
      move_step * ang_to_bohr / cell_vectors_lengths;
  //------interfaces----
  std::vector<double> step_size{relative_move_step(normal_direction),
                                relative_move_step(normal_direction)};
  std::vector<double> opt_param{interfaces(0), interfaces(1)};
  std::vector<double> low_b{0, 0}; // lower bounds
  std::vector<double> upp_b{1, 1}; // upper bounds
  if (!optimize.interfaces) {
    low_b = opt_param;
    upp_b = opt_param;
  }

  //------charge positions----
  for (arma::uword i = 0; i < charge_position.n_rows; ++i) {
    for (arma::uword j = 0; j < 3; ++j) {
      opt_param.push_back(charge_position(i, j));
      if (optimize.charge_position) {
        low_b.push_back(0);
        upp_b.push_back(1);
      } else {
        low_b.push_back(charge_position(i, j));
        upp_b.push_back(charge_position(i, j));
      }
    }
    step_size.insert(
        step_size.end(),
        {relative_move_step(0), relative_move_step(1), relative_move_step(2)});
    //------charge_sigma----
    opt_param.insert(opt_param.end(), {charge_sigma(i, 0), charge_sigma(i, 1),
                                       charge_sigma(i, 2)});
    step_size.insert(step_size.end(), {sigma_step, sigma_step, sigma_step});
    if (optimize.charge_sigma) {
      if (trivariate_charge) {
        low_b.insert(low_b.end(), {0.1, 0.1, 0.1});
        upp_b.insert(upp_b.end(), {7, 7, 7});
      } else {
        low_b.insert(low_b.end(),
                     {0.1, charge_sigma(i, 1), charge_sigma(i, 2)});
        upp_b.insert(upp_b.end(), {7, charge_sigma(i, 1), charge_sigma(i, 2)});
      }
    } else {
      low_b.insert(low_b.end(), {charge_sigma(i, 0), charge_sigma(i, 1),
                                 charge_sigma(i, 2)});
      upp_b.insert(upp_b.end(), {charge_sigma(i, 0), charge_sigma(i, 1),
                                 charge_sigma(i, 2)});
    }
    //-----charge_rotations-----
    opt_param.insert(opt_param.end(),
                     {charge_rotations(i, 0), charge_rotations(i, 1),
                      charge_rotations(i, 2)});
    step_size.insert(step_size.end(),
                     {rotation_step, rotation_step, rotation_step});
    if (optimize.charge_rotation) {
      if (trivariate_charge) {
        low_b.insert(low_b.end(), {-0.5 * PI, -0.5 * PI, -0.5 * PI});
        upp_b.insert(upp_b.end(), {0.5 * PI, 0.5 * PI, 0.5 * PI});
      } else {
        // this must be caught in the checker!
        log->warn("Optimizing the rotation angles for the simple Gaussian "
                  "charges is not possible!");
        low_b.insert(low_b.end(),
                     {charge_rotations(i, 0), charge_rotations(i, 1),
                      charge_rotations(i, 2)});
        upp_b.insert(upp_b.end(),
                     {charge_rotations(i, 0), charge_rotations(i, 1),
                      charge_rotations(i, 2)});
      }
    } else {
      low_b.insert(low_b.end(), {charge_rotations(i, 0), charge_rotations(i, 1),
                                 charge_rotations(i, 2)});
      upp_b.insert(upp_b.end(), {charge_rotations(i, 0), charge_rotations(i, 1),
                                 charge_rotations(i, 2)});
    }
    //------charge fractions----
    opt_param.push_back(charge_fraction(i));
    step_size.push_back(fraction_step);
    if (optimize.charge_fraction) {
      low_b.push_back(0);
      upp_b.push_back(1);
    } else {
      low_b.push_back(charge_fraction(i));
      upp_b.push_back(charge_fraction(i));
    }
  }
  // remove the last charge_fraction from the parameters
  // this gives a better chance to optimizer and also removes the
  // charge_fraction if we have only one charge
  opt_param.pop_back();
  upp_b.pop_back();
  low_b.pop_back();
  step_size.pop_back();
  return make_tuple(opt_param, low_b, upp_b, step_size);
}

void slabcc_model::data_unpacker(
    const std::vector<double> &optimizer_vars_vec) {

  // optimizer_vars_vec is ordered as: 2x interface, :|[x, y, z, 3x sigma, 3x
  // rotation, fraction]|:
  const arma::uword position_per_charge = 3;
  const arma::uword sigma_per_charge = 3;
  const arma::uword rotation_per_charge = 3;
  const arma::uword fraction_per_charge = 1;

  const arma::uword position_offset = 2; // interfaces
  const arma::uword sigma_offset = position_offset + position_per_charge;
  const arma::uword rotation_offset = sigma_offset + sigma_per_charge;
  const arma::uword fraction_offset = rotation_offset + rotation_per_charge;

  const arma::uword variables_per_charge =
      position_per_charge + sigma_per_charge + rotation_per_charge +
      fraction_per_charge;
  const arma::uword n_charges =
      optimizer_vars_vec.size() / variables_per_charge;
  charge_fraction.zeros();
  interfaces = {optimizer_vars_vec.at(0), optimizer_vars_vec.at(1)};
  for (arma::uword i = 0; i < n_charges; ++i) {
    for (arma::uword j = 0; j < 3; ++j) {
      charge_position(i, j) =
          optimizer_vars_vec.at(position_offset + i * variables_per_charge + j);
      charge_sigma(i, j) =
          optimizer_vars_vec.at(sigma_offset + i * variables_per_charge + j);
      charge_rotations(i, j) =
          optimizer_vars_vec.at(rotation_offset + i * variables_per_charge + j);
    }
    if (i != n_charges - 1) {
      charge_fraction(i) =
          optimizer_vars_vec.at(fraction_offset + i * variables_per_charge);
    } else {
      charge_fraction(i) = 1.0 - accu(charge_fraction);
    }
  }
}

void slabcc_model::verify_interface_optimization(
    const arma::rowvec2 &initial_interfaces) const {

  const arma::rowvec2 unshifted_initial_interfaces =
      fmod_p(initial_interfaces - rounded_relative_shift(normal_direction), 1);
  const arma::rowvec2 unshifted__interfaces =
      fmod_p(interfaces - rounded_relative_shift(normal_direction), 1);

  const double max_total_safe_movement = 4; // Ang
  const double normal_length =
      cell_vectors_lengths(normal_direction) / ang_to_bohr;
  auto log = spdlog::get("loggers");

  const arma::rowvec3 mirroring_template = {-1, 0, 1};
  const arma::mat interface_all_mirrors =
      repmat(mirroring_template, 2, 1) +
      repmat(unshifted__interfaces.t(), 1, 3);
  const arma::mat interface_n0_shiftes =
      abs(interface_all_mirrors - unshifted_initial_interfaces(0));
  arma::vec interface_n0_distances = min(interface_n0_shiftes, 1);
  const arma::mat interface_n1_shiftes =
      abs(interface_all_mirrors - unshifted_initial_interfaces(1));
  arma::vec interface_n1_distances = min(interface_n1_shiftes, 1);
  const arma::rowvec2 total_changes = {
      interface_n0_distances(0) + interface_n1_distances(1),
      interface_n0_distances(1) + interface_n1_distances(0)};

  // convert to Ang
  interface_n0_distances *= normal_length;
  interface_n1_distances *= normal_length;

  if (total_changes(0) < total_changes(1)) {
    log->trace("The order of the interfaces has not been swapped!");
    log->trace("1st interface has moved {} angstrom from {} to {}",
               interface_n0_distances(0), unshifted_initial_interfaces(0),
               unshifted__interfaces(0));
    log->trace("2nd interface has moved {} angstrom from {} to {}",
               interface_n1_distances(1), unshifted_initial_interfaces(1),
               unshifted__interfaces(1));
  } else {
    log->trace("The order of the interfaces has been swapped!");
    log->trace("1st interface has moved {} angstrom from {} to {}",
               interface_n0_distances(1), unshifted_initial_interfaces(0),
               unshifted__interfaces(1));
    log->trace("2nd interface has moved {} angstrom from {} to {}",
               interface_n1_distances(0), unshifted_initial_interfaces(1),
               unshifted__interfaces(0));
  }

  if (total_changes(0) * normal_length > max_total_safe_movement) {
    log->warn("There is a relatively large change in the position of the "
              "interfaces after the optimization! Please ensure the "
              "correctness of the interfaces.");
  }
}

void slabcc_model::verify_charge_optimization() const {
  // maximum sigma for the localized charges
  const double max_sigma = 6.5;
  const double min_non_negligable_q = 0.01;
  auto log = spdlog::get("loggers");
  for (arma::uword i = 0; i < charge_fraction.n_elem; ++i) {
    if (abs(total_charge) * charge_fraction(i) > min_non_negligable_q) {
      if (charge_sigma.row(i).max() > max_sigma) {
        log->error("The model extra charge seems to be (at least partially) "
                   "delocalized. The current charge correction method is not "
                   "suitable for the delocalized charges.");
        break;
      }
    }
  }
}

bool slabcc_model::had_discretization_error() {
  if (in_optimization)
    return false;

  auto log = spdlog::get("loggers");
  const double tolerance = 1e-4; // minimum significant discretization error
  const double new_charge_error = abs(defect_charge - total_charge);
  if ((last_charge_error > tolerance) &&
      (new_charge_error > last_charge_error)) {
    // increasing the grid size is not helping
    log->debug("Model charge error on the new grid size: {}", new_charge_error);
    log->critical("Increasing the calculation grid size did not decrease the "
                  "discretization error. Most probably the model charge is "
                  "fairly delocalized!");

    if (is_active(verbosity::write_planarAvg_file)) {
      write_planar_avg(arma::real(POT) * Hartree_to_eV,
                       arma::real(CHG) * voxel_vol, "M", cell_vectors_lengths);
    } else if (is_active(verbosity::write_normal_planarAvg)) {
      write_planar_avg(arma::real(POT) * Hartree_to_eV,
                       arma::real(CHG) * voxel_vol, "M", cell_vectors_lengths,
                       normal_direction);
    }

    finalize_loggers();
    exit(1);
  } else {
    if (new_charge_error > tolerance) {
      log->debug("Model charge error on the former grid size: {}",
                 new_charge_error);
      last_charge_error = new_charge_error;
      const arma::urowvec3 new_grid = cell_grid + cell_grid / 2;
      change_grid(new_grid);
      log->debug("New model charge grid size: {}", to_string(cell_grid));
      return true;
    } else {
      last_charge_error = new_charge_error;
      return false;
    }
  }
}

void slabcc_model::update_V_target() {
  auto log = spdlog::get("loggers");
  if (as_size(cell_grid) != arma::size(POT_target)) {
    POT_target.set_size(as_size(cell_grid));
    const arma::rowvec new_grid_x = arma::linspace<arma::rowvec>(
        1.0, POT_target_on_input_grid.n_rows, cell_grid(0));
    const arma::rowvec new_grid_y = arma::linspace<arma::rowvec>(
        1.0, POT_target_on_input_grid.n_cols, cell_grid(1));
    const arma::rowvec new_grid_z = arma::linspace<arma::rowvec>(
        1.0, POT_target_on_input_grid.n_slices, cell_grid(2));

    POT_target =
        interp3(POT_target_on_input_grid, new_grid_x, new_grid_y, new_grid_z);
    POT_target -= arma::accu(POT_target) / POT_target.n_elem;
    log->debug("New potential grid size: " + to_string(SizeVec(POT_target)));
  }
}

void slabcc_model::adjust_extrapolation_grid(
    const int &extrapol_steps_num, const double &extrapol_steps_size) {

  auto log = spdlog::get("loggers");
  log->trace("Checking the extrapolation grid size");
  const arma::mat33 cell_vectors0 = cell_vectors;
  const double total_charge0 = total_charge;
  // Force discretization error checks
  for (auto step = extrapol_steps_num - 1; step > 0; --step) {
    const double extrapol_factor = extrapol_steps_size * step + 1;
    change_size(cell_vectors0 * extrapol_factor);
    gaussian_charges_gen();
  }
  change_size(cell_vectors0);
  total_charge = total_charge0;
}

std::tuple<arma::rowvec, arma::rowvec>
slabcc_model::extrapolate(int extrapol_steps_num, double extrapol_steps_size) {

  auto log = spdlog::get("loggers");
  const arma::mat33 cell_vectors0 = cell_vectors;
  const arma::rowvec3 cell_vectors_lengths0 = cell_vectors_lengths;
  const auto interfaces0 = interfaces;
  const auto charge_position0 = charge_position;
  const double slab_thickness = abs(interfaces(0) - interfaces(1));
  arma::rowvec Es = arma::zeros<arma::rowvec>(extrapol_steps_num - 1),
               sizes = Es;
  for (auto step = 1; step < extrapol_steps_num; ++step) {
    const double extrapol_factor = extrapol_steps_size * step + 1;
    change_size(cell_vectors0 * extrapol_factor);
    if (this->type == model_type::slab) {
      // increase the slab thickness
      const arma::uvec interface_sorted_i = sort_index(interfaces);
      interfaces(interface_sorted_i(1)) =
          interfaces(interface_sorted_i(0)) + slab_thickness;
      // move the charges to the same distance from their original nearest
      // interface
      for (arma::uword charge_i = 0; charge_i < charge_position.n_rows;
           ++charge_i) {
        const arma::rowvec2 initial_distance_to_interfaces =
            (charge_position0(charge_i, normal_direction) - interfaces0) *
            cell_vectors_lengths0(normal_direction);
        if (abs(initial_distance_to_interfaces(0)) <
            abs(initial_distance_to_interfaces(1))) {
          charge_position(charge_i, normal_direction) =
              interfaces(0) + initial_distance_to_interfaces(0) /
                                  cell_vectors_lengths(normal_direction);
        } else {
          charge_position(charge_i, normal_direction) =
              interfaces(1) + initial_distance_to_interfaces(1) /
                                  cell_vectors_lengths(normal_direction);
          ;
        }
      }
    }

    gaussian_charges_gen();
    dielectric_profiles_gen();

    // (only works for the orthogonal cells!)
    const auto CHG_normalized =
        CHG - total_charge / arma::prod(cell_vectors_lengths);
    const auto V = poisson_solver_3D(CHG_normalized, dielectric_profiles,
                                     cell_vectors_lengths, normal_direction);
    const auto EperModel = 0.5 * arma::accu(arma::real(V % CHG_normalized)) *
                           voxel_vol * Hartree_to_eV;
    const arma::rowvec2 interface_pos =
        interfaces * cell_vectors_lengths(normal_direction);
    std::string extrapolation_info =
        to_string(extrapol_factor) + "\t" + to_string(EperModel) + "\t" +
        ::to_string(total_charge) + "\t" + to_string(interface_pos);
    for (arma::uword i = 0; i < charge_position.n_rows; ++i) {
      extrapolation_info +=
          "\t" + to_string(charge_position(i, normal_direction) *
                           cell_vectors_lengths(normal_direction));
    }
    log->debug(extrapolation_info);
    Es(step - 1) = EperModel;
    sizes(step - 1) = 1.0 / extrapol_factor;
  }

  return std::make_tuple(Es, sizes);
}

double slabcc_model::Eiso_bessel() const {

  auto logger = spdlog::get("loggers");
  const double K_min = 0.00001;
  const double K_max = 100;
  logger->debug("Isolated energy integration limits in the k-space: {}:{}",
                K_min, K_max);
  // integration grid
  arma::rowvec K;
  // hand-tuned adaptive grid based on the derivative of the Uk ~~log(k)
  for (int ki = 0; ki < K_max; ++ki) {
    const int grid_density =
        static_cast<int>(pow(10, -log((ki + 1) / 100.0)) / 5.0) + 1;
    K = arma::join_horiz(
        K, arma::linspace<arma::rowvec>(ki + K_min, ki + 1, grid_density));
  }
  logger->debug("Number of k-space integration grid points: {}", K.n_elem);

  // eq. 7 in the SI (Supplementary Information for `First-principles
  // electrostatic potentials for reliable alignment at interfaces and defects`)
  const arma::rowvec integrand =
      K % exp(-square(K) * pow(charge_sigma(0, 0), 2)) % Uk(K);
  const arma::mat integral = trapz(K, integrand, 1);
  const double Q = charge_fraction(0) * total_charge;
  const double U_total = pow(Q, 2) * integral(0) * Hartree_to_eV;

  return U_total;
}

arma::rowvec slabcc_model::Uk(arma::rowvec k) const {
  const double z0 = charge_position(0, normal_direction) *
                    cell_vectors_lengths(normal_direction);
  const arma::rowvec3 length = cell_vectors_lengths;
  const arma::urowvec3 n_points = cell_grid;
  const arma::uword normal = normal_direction;
  const arma::rowvec Gs = 2.0 * PI / length;

  arma::rowvec Gz0 = arma::ceil(arma::regspace<arma::rowvec>(
                         -0.5 * n_points(normal), 0.5 * n_points(normal) - 1)) *
                     Gs(normal);
  Gz0 = ifftshift(Gz0);
  const arma::rowvec Gz02 = arma::square(Gz0);

  const arma::uword LGz = Gz0.n_elem;
  const double dielbulk = 1; // bulk response far away in vaccum!
  const arma::cx_mat dielsG = fft(dielectric_profiles - dielbulk);
  const arma::cx_mat dielGz = arma::circ_toeplitz(dielsG.col(normal)) / LGz;
  const arma::uword inplane_direction = (normal == 0) ? 1 : 0;
  const arma::cx_mat dielpGz =
      arma::circ_toeplitz(dielsG.col(inplane_direction)) / LGz;

  const arma::cx_mat Ag1 = dielGz;
  const arma::cx_mat Ag1p = dielpGz;
  const arma::mat Ag2 = Gz0.t() * Gz0;
  const arma::cx_double coef(0, -1); // compiler-specific problem!
  const arma::cx_mat rhok =
      exp(coef * Gz0 * z0 - Gz02 * pow(charge_sigma(0, 0), 2) / 2.0);
  const arma::cx_mat rhok_t = rhok.t();
  const arma::cx_mat rho = ifft(rhok_t);
  arma::rowvec Uk = arma::zeros(arma::size(k));

  const arma::cx_mat Ag12 = Ag1 % Ag2;
  const auto cosGL_2 = cos(Gz0 * length(normal) / 2.0);
  for (arma::uword i = 0; i < k.n_elem; ++i) {
    const arma::cx_mat Ag = Ag12 + Ag1p * k(i) * k(i);
    const double keff = k(i);
    const arma::mat Kinvg =
        arma::diagmat(dielbulk * length(normal) * (pow(keff, 2) + Gz02) /
                      (1 - exp(-keff * length(normal) / 2.0) * cosGL_2));
    const arma::cx_mat Dg = Kinvg + length(normal) * Ag;
    const auto VGz = arma::solve(Dg, rhok_t);
    const auto Vz = ifft(VGz) * LGz;
    Uk(i) = real(arma::accu(Vz % rho));
  }

  return Uk;
}

double potential_error(const std::vector<double> &x, std::vector<double> &grad,
                       void *model_ptr) {
  slabcc_model &model = *static_cast<slabcc_model *>(model_ptr);
  return model.potential_error(x, grad);
}

double slabcc_model::potential_error(const std::vector<double> &x,
                                     std::vector<double> &grad) {
  auto log = spdlog::get("loggers");

  data_unpacker(x);
  arma::rowvec normalized_charge_fraction = charge_fraction;

  // total charge error
  double bounds_factor = 0;

  if (charge_fraction(charge_fraction.n_elem - 1) < 0) {
    bounds_factor = -charge_fraction(charge_fraction.n_elem - 1);
    normalized_charge_fraction(normalized_charge_fraction.n_elem - 1) = 0;
    normalized_charge_fraction /= (bounds_factor + 1);
  }

  gaussian_charges_gen();
  dielectric_profiles_gen();

  POT = poisson_solver_3D(CHG, dielectric_profiles, cell_vectors_lengths,
                          normal_direction);
  POT_diff = real(POT) * Hartree_to_eV - POT_target;
  // bigger output for out-of-bounds input: quadratic penalty
  const double bounds_correction =
      bounds_factor + 10 * bounds_factor * bounds_factor;
  potential_RMSE = sqrt(arma::accu(arma::square(POT_diff)) / POT_diff.n_elem) +
                   bounds_correction;

  if (initial_potential_RMSE < 0) {
    initial_potential_RMSE = potential_RMSE;
  }

  log->debug("-----------------------------------------");
  if (this->type != model_type::bulk) {
    const arma::rowvec2 unshifted_interfaces =
        fmod_p(interfaces - rounded_relative_shift(normal_direction), 1);
    log->debug(" > interfaces={}", to_string(unshifted_interfaces));
  }

  const arma::mat unshifted_charge_position =
      fmod_p(charge_position - arma::repmat(rounded_relative_shift,
                                            charge_position.n_rows, 1),
             1);

  for (arma::uword i = 0; i < normalized_charge_fraction.n_elem; ++i) {
    log->debug("{}> charge_position={}", i + 1,
               to_string(unshifted_charge_position.row(i)));
    if (trivariate_charge) {
      log->debug("{}> charge_sigma={}", i + 1, to_string(charge_sigma.row(i)));
      if (abs(charge_rotations).max() > 0) {
        const arma::rowvec3 rotation = charge_rotations.row(i) * 180.0 / PI;
        log->debug("{}> charge_rotation={}", i + 1, to_string(rotation));
      }
    } else {
      log->debug("{}> charge_sigma={}", i + 1, to_string(charge_sigma(i, 0)));
    }
    if (charge_fraction.n_elem > 1) {
      log->debug("{}> charge_fraction={}", i + 1, charge_fraction(i));
    }
  }
  if (bounds_correction > 0) {
    log->debug("Out of the bounds correction to the RMSE: {}",
               bounds_correction);
  }
  log->debug("Potential Root Mean Square Error: {}", potential_RMSE);

  return potential_RMSE;
}

void slabcc_model::optimize(const std::string &opt_algo, const double &opt_tol,
                            const int &max_eval, const int &max_time,
                            const opt_switches &optimize) {

  auto log = spdlog::get("loggers");
  in_optimization = true;

  auto opt_algorithm = nlopt::LN_COBYLA;
  if (opt_algo == "BOBYQA") {
    opt_algorithm = nlopt::LN_BOBYQA;
  } else if (opt_algo == "SBPLX") {
    opt_algorithm = nlopt::LN_SBPLX;
  }

  std::vector<double> opt_param, low_b, upp_b, step_size;
  tie(opt_param, low_b, upp_b, step_size) = data_packer(optimize);
  nlopt::opt opt(opt_algorithm, opt_param.size());
  opt.set_lower_bounds(low_b);
  opt.set_upper_bounds(upp_b);
  opt.set_initial_step(step_size);
  opt.set_min_objective(::potential_error, this);
  opt.set_xtol_rel(opt_tol);
  if (max_eval > 0) {
    opt.set_maxeval(max_eval);
  }
  if (max_time > 0) {
    opt.set_maxtime(60.0 * max_time);
  }

  const int sigma_per_charge = trivariate_charge ? 3 : 1;
  const int var_per_charge =
      static_cast<int>(optimize.charge_position) * 3 +
      static_cast<int>(optimize.charge_rotation) * 3 +
      static_cast<int>(optimize.charge_sigma) * sigma_per_charge +
      static_cast<int>(optimize.charge_fraction) * 1;
  const arma::uword opt_parameters =
      charge_fraction.n_elem * var_per_charge + 2 * optimize.interfaces;
  log->trace("Started optimizing {} model parameters", opt_parameters);
  log->trace("Optimization algorithm: " +
             std::string(opt.get_algorithm_name()));
  try {
    const nlopt::result nlopt_final_result =
        opt.optimize(opt_param, potential_RMSE);
    log->debug("-----------------------------------------");
    if (nlopt_final_result == nlopt::MAXEVAL_REACHED) {
      log->warn("Optimization ended after {} steps before reaching the "
                "requested accuracy!",
                max_eval);
    } else if (nlopt_final_result == nlopt::MAXTIME_REACHED) {
      log->warn("Optimization ended after {} minutes before reaching the "
                "requested accuracy!",
                max_time);
    }
  } catch (const std::exception &e) {
    log->error("Optimization of the slabcc parameters failed: " +
               std::string(e.what()));
    log->error("Please start with better initial guess for the input "
               "parameters or use a different optimization algorithm.");
  }

  data_unpacker(opt_param);
  in_optimization = false;
  log->trace("Optimization ended.");
}

void slabcc_model::check_V_error() {
  auto log = spdlog::get("loggers");

  const bool isotropic_screening = accu(abs(diff(diel_in))) < 0.02;
  if (potential_RMSE > 0.1) {
    if (type == model_type::bulk && isotropic_screening) {
      log->debug("RMSE of the model charge potential is large but for the bulk "
                 "models with an isotropic screening (dielectric tensor) "
                 "this shouldn't make much difference in the total correction "
                 "energy!");
    } else {
      log->warn("RMSE of the model charge potential is large. The calculated "
                "correction energies may not be accurate!");
    }
  }

  const auto V_error_x = planar_average(0, POT_diff);
  const auto V_error_y = planar_average(1, POT_diff);
  const auto V_error_z = planar_average(2, POT_diff);

  // potential error in each direction
  arma::rowvec3 V_error_planars = {arma::accu(arma::square(V_error_x)),
                                   arma::accu(arma::square(V_error_y)),
                                   arma::accu(arma::square(V_error_z))};
  V_error_planars = arma::sqrt(V_error_planars / POT_diff.n_elem);
  log->debug("Directional RMSE: " + to_string(V_error_planars));
  if (arma::max(V_error_planars) / arma::min(V_error_planars) > 10) {
    log->warn("The potential error is highly anisotropic.");
    log->warn(
        "If the potential error is large, this usually means that either the "
        "extra charge is not properly described by the model Gaussian charge "
        "or the chosen dielectric tensor is not a good representation of the "
        "actual tensor! "
        "This can be fixed by properly optimizing the model parameters, using "
        "multiple Gaussian charges, using trivaritate Gaussians, or using the "
        "dielectric tensor calculated for the same VASP model/method.");
  }

  log->debug("Potential error anisotropy: {}",
             max(V_error_planars) / min(V_error_planars));
}

void slabcc_model::verify_CHG(const arma::cube &defect_charge) {
  auto log = spdlog::get("loggers");

  if (this->type != model_type::bulk) {

    // find the grid index of the interfaces
    arma::rowvec2 interfaces_index = cell_grid(normal_direction) * interfaces;
    arma::urowvec2 interfaces_grid_i = {
        static_cast<arma::uword>(interfaces_index(0)),
        static_cast<arma::uword>(interfaces_index(1))};
    interfaces_grid_i = arma::sort(interfaces_grid_i);
    std::vector<arma::span> spans = {
        arma::span(), arma::span(),
        arma::span(interfaces_grid_i(0), interfaces_grid_i(1))};
    std::swap(spans[normal_direction], spans[2]);

    const double model_total = arma::accu(arma::real(CHG)) * voxel_vol;
    const double model_in =
        arma::accu(arma::real(CHG(spans[0], spans[1], spans[2]))) * voxel_vol;
    const double model_out = model_total - model_in;

    // The original defect may have a different grid size than the final model!
    const arma::urowvec3 defect_grid = SizeVec(defect_charge);
    arma::rowvec2 defect_interfaces_index =
        defect_grid(normal_direction) * interfaces;
    arma::urowvec2 defect_interfaces_grid_i = {
        static_cast<arma::uword>(defect_interfaces_index(0)),
        static_cast<arma::uword>(defect_interfaces_index(1))};
    defect_interfaces_grid_i = sort(defect_interfaces_grid_i);
    std::vector<arma::span> defect_spans = {
        arma::span(), arma::span(),
        arma::span(defect_interfaces_grid_i(0), defect_interfaces_grid_i(1))};
    std::swap(defect_spans[normal_direction], defect_spans[2]);

    const double defect_voxel_vol =
        arma::prod(cell_vectors_lengths / defect_grid);
    const double defect_total = arma::accu(defect_charge) * defect_voxel_vol;
    const double defect_in =
        accu(defect_charge(defect_spans[0], defect_spans[1], defect_spans[2])) *
        defect_voxel_vol;
    const double defect_out = defect_total - defect_in;

    // these two do NOT need to closely agree with each other if the charge is
    // near the surface
    log->debug("Total charge of the model slab (inside, outside): {}, {}",
               model_in, model_out);
    log->debug("Total charge of the defect slab (inside, outside): {}, {}",
               defect_in, defect_out);

    const bool large_charge_difference = abs(model_in - defect_in) > 0.1;
    const arma::mat relative_charge_interface_distance =
        abs(arma::repmat(charge_position.col(normal_direction), 1, 2) -
            arma::repmat(interfaces, charge_position.n_rows, 1));
    const arma::mat charge_interface_distance =
        relative_charge_interface_distance *
        cell_vectors_lengths(normal_direction) / ang_to_bohr;
    const bool charge_near_surface = charge_interface_distance.min() < 1;

    if (large_charge_difference && !charge_near_surface) {
      log->debug("Minimum distance of the model charge centers to the surface "
                 "(Ang): {}",
                 charge_interface_distance.min());
      log->debug("Difference of the charge inside of the slab in the defect "
                 "and the model: {}",
                 abs(model_in - defect_in));
      log->warn(
          "The amount of the charge inside and outside of the slab seems to be "
          "different in the input files and the slabcc model. The extra charge "
          "may be partially localized on the slab surfaces. ");
    }
  }
}