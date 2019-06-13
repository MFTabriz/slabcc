// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_model.hpp"


void slabcc_model::set_input_variables(const input_data& inputfile_variables) {
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

void slabcc_model::set_model_type(const bool& model_2d, const rowvec& diel_in, const rowvec& diel_out) {
	if (approx_equal(diel_in, diel_out, "absdiff", 0.02)) {
		this->type = model_type::bulk;
	}
	else if (model_2d) {
		this->type = model_type::monolayer;
	}
	else {
		this->type = model_type::slab;
	}
}

void slabcc_model::init_supercell(const mat33& new_vectors, const urowvec3& new_grid) {
	cell_vectors = new_vectors;
	cell_grid = new_grid;
	update_cell_vectors_lengths();
	update_voxel_vol();
}

void slabcc_model::change_grid(const urowvec3& new_cell_grid) {
	cell_grid = new_cell_grid;
	update_voxel_vol();
}

void slabcc_model::change_size(const mat33& new_cell_vectors) {
	auto log = spdlog::get("loggers");
	const rowvec3 cell_vectors_lengths0 = cell_vectors_lengths;
	cell_vectors = new_cell_vectors;
	update_cell_vectors_lengths();
	update_voxel_vol();

	interfaces *= cell_vectors_lengths0(normal_direction) / cell_vectors_lengths(normal_direction);
	charge_position *= cell_vectors_lengths0(normal_direction) / cell_vectors_lengths(normal_direction);

	const rowvec3 scaling = cell_vectors_lengths / cell_vectors_lengths0;
	if (abs(max(scaling) - min(scaling)) > 0.00001) {
		log->debug("Model cell scaling is anisotropic!");
	}
}

void slabcc_model::update_voxel_vol() {
	voxel_vol = prod(cell_vectors_lengths / cell_grid);
}

void slabcc_model::update_cell_vectors_lengths() {
	for (uword i = 0; i < 3; ++i) {
		cell_vectors_lengths(i) = norm(cell_vectors.col(i));
	}
}

void slabcc_model::dielectric_profiles_gen() {
	const auto length = cell_vectors_lengths(normal_direction);
	const auto n_points = cell_grid(normal_direction);
	rowvec2 interfaces_cartesian = interfaces * length;
	interfaces_cartesian = sort(interfaces_cartesian);
	const auto positions = linspace<rowvec>(0, length, n_points + 1);
	dielectric_profiles = arma::zeros<mat>(n_points, 3);
	const rowvec3 diel_sum = diel_in + diel_out;
	const rowvec3 diel_diff = diel_out - diel_in;

	for (uword k = 0; k < n_points; ++k) {
		double min_distance = 0;

		double diel_side = 1;
		// min distance of the point positions(k) to the interfaces in PBC [-length/2  length/2]
		const rowvec2 distances = fmod_p(positions(k) - interfaces_cartesian + length / 2, length) - length / 2;
		if (abs(distances(0)) < abs(distances(1))) {
			min_distance = distances(0);
			diel_side = -1;
		}
		else {
			min_distance = distances(1);
		}

		const auto diel_edge = erf(min_distance / diel_erf_beta);
		dielectric_profiles.row(k) = (diel_diff * diel_side * diel_edge + diel_sum) / 2;
	}

}

void slabcc_model::gaussian_charges_gen() {

	do {
		rowvec x0 = linspace<rowvec>(0, cell_vectors_lengths(0) - cell_vectors_lengths(0) / cell_grid(0), cell_grid(0));
		rowvec y0 = linspace<rowvec>(0, cell_vectors_lengths(1) - cell_vectors_lengths(1) / cell_grid(1), cell_grid(1));
		rowvec z0 = linspace<rowvec>(0, cell_vectors_lengths(2) - cell_vectors_lengths(2) / cell_grid(2), cell_grid(2));

		CHG = arma::zeros<cx_cube>(as_size(cell_grid));

		for (uword i = 0; i < charge_fraction.n_elem; ++i) {
			// shift the axis reference to position of the Gaussian charge center
			rowvec x = x0 - accu(cell_vectors.col(0) * charge_position(i, 0));
			rowvec y = y0 - accu(cell_vectors.col(1) * charge_position(i, 1));
			rowvec z = z0 - accu(cell_vectors.col(2) * charge_position(i, 2));
			//handle the minimum distance from the mirror charges
			for (auto& pos : x) {
				if (abs(pos) > cell_vectors_lengths(0) / 2) {
					pos = cell_vectors_lengths(0) - abs(pos);
				}
			}
			for (auto& pos : y) {
				if (abs(pos) > cell_vectors_lengths(1) / 2) {
					pos = cell_vectors_lengths(1) - abs(pos);
				}
			}
			for (auto& pos : z) {
				if (abs(pos) > cell_vectors_lengths(2) / 2) {
					pos = cell_vectors_lengths(2) - abs(pos);
				}
			}

			cube xs, ys, zs;
			tie(xs, ys, zs) = ndgrid(x, y, z);

			const rowvec3 rotation_angle = charge_rotations.row(i);
			//rotate around xyz axis
			if (max(abs(rotation_angle)) > 0.002) {
				const mat33 rot_x = {
					{1, 0, 0},
					{0, cos(rotation_angle(0)), -sin(rotation_angle(0))},
					{0, sin(rotation_angle(0)), cos(rotation_angle(0))}
				};

				const mat33 rot_y = {
					{ cos(rotation_angle(1)), 0, sin(rotation_angle(1))},
					{0, 1, 0},
					{-sin(rotation_angle(1)), 0, cos(rotation_angle(1))}
				};

				const mat33 rot_z = {
					{cos(rotation_angle(2)), -sin(rotation_angle(2)), 0},
					{sin(rotation_angle(2)), cos(rotation_angle(2)), 0},
					{0, 0, 1}
				};

				const mat33 rotation_mat = rot_x * rot_y * rot_z;
				for (uword i = 0; i < xs.n_elem; ++i) {
					const vec3 old_coordinates = { xs(i), ys(i), zs(i) };
					const vec3 new_coordinates = rotation_mat * old_coordinates;
					xs(i) = new_coordinates(0);
					ys(i) = new_coordinates(1);
					zs(i) = new_coordinates(2);
				}
			}

			const cube r2 = square(xs) + square(ys) + square(zs);

			// this charge distribution is due to the 1st nearest gaussian image. 
			// In case of the very small supercells or very diffuse charges (large sigma), the higher order of the image charges must also be included.
			// But the validity of the correction method for these cases must be checked!	

			const double Q = charge_fraction(i) * defect_charge;
			if (trivariate_charge) {
				CHG += cx_cube(Q / (pow(2 * PI, 1.5) * prod(charge_sigma.row(i)))
					* exp(-square(xs) / (2 * square(charge_sigma(i, 0))) - square(ys) / (2 * square(charge_sigma(i, 1))) - square(zs) / (2 * square(charge_sigma(i, 2))))
					, arma::zeros(as_size(cell_grid)));
			}
			else {
				CHG += cx_cube(Q / pow((charge_sigma(i, 0) * sqrt(2 * PI)), 3) * exp(-r2 / (2 * square(charge_sigma(i, 0)))), arma::zeros(as_size(cell_grid)));
			}
		}

		total_charge = accu(real(CHG)) * voxel_vol;
	}while(had_discretization_error());
	update_V_target();
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>> slabcc_model::data_packer(opt_switches optimize) const {
	auto log = spdlog::get("loggers");
	//size of the first step for each parameter
	const double move_step = 1; //Ang
	const double sigma_step = 0.5;
	const double fraction_step = 0.2;
	const double rotation_step = 0.2; //rad
	const rowvec3 relative_move_step = move_step * ang_to_bohr / cell_vectors_lengths;
	//------interfaces----
	vector<double> step_size = { relative_move_step(normal_direction) , relative_move_step(normal_direction) };
	vector<double> opt_param = { interfaces(0), interfaces(1) };
	vector<double> low_b = { 0, 0 };		//lower bounds
	vector<double> upp_b = { 1, 1 };		//upper bounds
	if (!optimize.interfaces) {
		low_b = opt_param;
		upp_b = opt_param;
	}

	//------charge positions----
	for (uword i = 0; i < charge_position.n_rows; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_param.push_back(charge_position(i, j));
			if (optimize.charge_position) {
				low_b.push_back(0);
				upp_b.push_back(1);
			}
			else {
				low_b.push_back(charge_position(i, j));
				upp_b.push_back(charge_position(i, j));
			}
		}
		step_size.insert(step_size.end(), { relative_move_step(0), relative_move_step(1), relative_move_step(2) });
		//------charge_sigma----
		opt_param.insert(opt_param.end(), { charge_sigma(i, 0), charge_sigma(i, 1), charge_sigma(i, 2) });
		step_size.insert(step_size.end(), { sigma_step,sigma_step,sigma_step });
		if (optimize.charge_sigma) {
			if (trivariate_charge) {
				low_b.insert(low_b.end(), { 0.1, 0.1, 0.1 });
				upp_b.insert(upp_b.end(), { 7, 7, 7 });
			}
			else {
				low_b.insert(low_b.end(), { 0.1, charge_sigma(i, 1), charge_sigma(i, 2) });
				upp_b.insert(upp_b.end(), { 7, charge_sigma(i, 1), charge_sigma(i, 2) });
			}

		}
		else {
			low_b.insert(low_b.end(), { charge_sigma(i, 0), charge_sigma(i, 1), charge_sigma(i, 2) });
			upp_b.insert(upp_b.end(), { charge_sigma(i, 0), charge_sigma(i, 1), charge_sigma(i, 2) });
		}
		//-----charge_rotations-----
		opt_param.insert(opt_param.end(), { charge_rotations(i, 0), charge_rotations(i, 1), charge_rotations(i, 2) });
		step_size.insert(step_size.end(), { rotation_step,rotation_step,rotation_step });
		if (optimize.charge_rotation) {
			if (trivariate_charge) {
				low_b.insert(low_b.end(), { -0.5 * PI , -0.5 * PI , -0.5 * PI });
				upp_b.insert(upp_b.end(), { 0.5 * PI , 0.5 * PI , 0.5 * PI });
			}
			else {
				//this must be caught in the checker!
				log->warn("Optimizing the rotation angles for the simple Gaussian charges is not possible!");
				low_b.insert(low_b.end(), { charge_rotations(i, 0), charge_rotations(i, 1), charge_rotations(i, 2) });
				upp_b.insert(upp_b.end(), { charge_rotations(i, 0), charge_rotations(i, 1), charge_rotations(i, 2) });
			}
		}
		else {
			low_b.insert(low_b.end(), { charge_rotations(i, 0), charge_rotations(i, 1), charge_rotations(i, 2) });
			upp_b.insert(upp_b.end(), { charge_rotations(i, 0), charge_rotations(i, 1), charge_rotations(i, 2) });
		}
		//------charge fractions----
		opt_param.push_back(charge_fraction(i));
		step_size.push_back(fraction_step);
		if (optimize.charge_fraction) {
			low_b.push_back(0);
			upp_b.push_back(1);
		}
		else {
			low_b.push_back(charge_fraction(i));
			upp_b.push_back(charge_fraction(i));
		}
	}
	//remove the last charge_fraction from the parameters 
	//this gives a better chance to optimizer and also removes the charge_fraction if we have only one charge
	opt_param.pop_back();
	upp_b.pop_back();
	low_b.pop_back();
	step_size.pop_back();
	return make_tuple(opt_param, low_b, upp_b, step_size);
}

void slabcc_model::data_unpacker(const vector<double>& optimizer_vars_vec) {

	// optimizer_vars_vec is ordered as: 2x interface, :|[x, y, z, 3x sigma, 3x rotation, fraction]|:
	const uword position_per_charge = 3;
	const uword sigma_per_charge = 3;
	const uword rotation_per_charge = 3;
	const uword fraction_per_charge = 1;

	const uword position_offset = 2; //interfaces
	const uword sigma_offset = position_offset + position_per_charge;
	const uword rotation_offset = sigma_offset + sigma_per_charge;
	const uword fraction_offset = rotation_offset + rotation_per_charge;

	const uword variables_per_charge = position_per_charge + sigma_per_charge + rotation_per_charge + fraction_per_charge;
	const uword n_charges = optimizer_vars_vec.size() / variables_per_charge;
	charge_fraction.zeros();
	interfaces = { optimizer_vars_vec.at(0), optimizer_vars_vec.at(1) };
	for (uword i = 0; i < n_charges; ++i) {
		for (uword j = 0; j < 3; ++j) {
			charge_position(i, j) = optimizer_vars_vec.at(position_offset + i * variables_per_charge + j);
			charge_sigma(i, j) = optimizer_vars_vec.at(sigma_offset + i * variables_per_charge + j);
			charge_rotations(i, j) = optimizer_vars_vec.at(rotation_offset + i * variables_per_charge + j);
		}
		if (i != n_charges - 1) {
			charge_fraction(i) = optimizer_vars_vec.at(fraction_offset + i * variables_per_charge);
		}
		else {
			charge_fraction(i) = 1.0 - accu(charge_fraction);
		}
	}
}

void slabcc_model::verify_interface_optimization(const rowvec2& initial_interfaces) const {

	const rowvec2 unshifted_initial_interfaces = fmod_p(initial_interfaces - rounded_relative_shift(normal_direction), 1);
	const rowvec2 unshifted__interfaces = fmod_p(interfaces - rounded_relative_shift(normal_direction), 1);

	const double max_total_safe_movement = 4; //Ang
	const double normal_length = cell_vectors_lengths(normal_direction) / ang_to_bohr;
	auto log = spdlog::get("loggers");

	const rowvec3 mirroring_template = { -1,0,1 };
	const mat interface_all_mirrors = repmat(mirroring_template, 2, 1) + repmat(unshifted__interfaces.t(), 1, 3);
	const mat interface_n0_shiftes = abs(interface_all_mirrors - unshifted_initial_interfaces(0));
	vec interface_n0_distances = min(interface_n0_shiftes, 1);
	const mat interface_n1_shiftes = abs(interface_all_mirrors - unshifted_initial_interfaces(1));
	vec interface_n1_distances = min(interface_n1_shiftes, 1);
	const rowvec2 total_changes = { interface_n0_distances(0) + interface_n1_distances(1),  interface_n0_distances(1) + interface_n1_distances(0) };

	//convert to Ang
	interface_n0_distances *= normal_length;
	interface_n1_distances *= normal_length;

	if (total_changes(0) < total_changes(1)) {
		log->trace("The order of the interfaces has not been swapped!");
		log->trace("1st interface has moved {} angstrom from {} to {}", interface_n0_distances(0), unshifted_initial_interfaces(0), unshifted__interfaces(0));
		log->trace("2nd interface has moved {} angstrom from {} to {}", interface_n1_distances(1), unshifted_initial_interfaces(1), unshifted__interfaces(1));
	}
	else {
		log->trace("The order of the interfaces has been swapped!");
		log->trace("1st interface has moved {} angstrom from {} to {}", interface_n0_distances(1), unshifted_initial_interfaces(0), unshifted__interfaces(1));
		log->trace("2nd interface has moved {} angstrom from {} to {}", interface_n1_distances(0), unshifted_initial_interfaces(1), unshifted__interfaces(0));
	}

	if (total_changes(0) * normal_length > max_total_safe_movement) {
		log->warn("There is a relatively large change in the position of the interfaces after the optimization! Please ensure the correctness of the interfaces.");
	}

}

void slabcc_model::verify_charge_optimization() const {
	//maximum sigma for the localized charges
	const double max_sigma = 6.5;
	const double min_non_negligable_q = 0.01;
	auto log = spdlog::get("loggers");
	for (uword i = 0; i < charge_fraction.n_elem; ++i) {
		if (abs(total_charge) * charge_fraction(i) > min_non_negligable_q) {
			if (charge_sigma.row(i).max() > max_sigma) {
				log->error("The model extra charge seems to be (at least partially) delocalized. The current charge correction method is not suitable for the delocalized charges.");
			}
		}
	}
}

bool slabcc_model::had_discretization_error() {
	auto log = spdlog::get("loggers");
	if (abs(defect_charge - total_charge) > 0.00001) {
		log->debug("The amount of the extra charge in model is not the same as the extra charge in the VASP input files.");
		log->debug("Model charge error: {}", abs(defect_charge - total_charge));
		if (charge_sigma.max() > 6) {
			log->critical("charge_sigma value is too large, the model charge is fairly delocalized and the model supercell cannot contain the whole Gaussian charge and the present charge correction method is not suitable for this cases!");
			finalize_loggers();
			exit(1);
		} else if (charge_sigma.min() < 0.3) {
			log->debug("Possible discretization error was detected!");
		}
		const rowvec3 new_grid_size = 1.2 * conv_to<rowvec>::from(cell_grid);
		const urowvec3 new_grid = { (uword)new_grid_size(0), (uword)new_grid_size(1), (uword)new_grid_size(2) };
		change_grid(new_grid);
		log->debug("Model charge grid size was adjusted to: {}", to_string(cell_grid));
		return true;
	}
	return false;
}

void slabcc_model::update_V_target() {
	auto log = spdlog::get("loggers");
	if (as_size(cell_grid) != arma::size(V_target)) {
		V_target.set_size(as_size(cell_grid));
		V_target.zeros();
		const rowvec new_grid_x = linspace<rowvec>(1.0, V_target0.n_rows, cell_grid(0));
		const rowvec new_grid_y = linspace<rowvec>(1.0, V_target0.n_cols, cell_grid(1));
		const rowvec new_grid_z = linspace<rowvec>(1.0, V_target0.n_slices, cell_grid(2));

		V_target = interp3(V_target0, new_grid_x, new_grid_y, new_grid_z);
		V_target -= accu(V_target) / V_target.n_elem;
		log->debug("New potential grid size: " + to_string(SizeVec(V_target)));
	}
}

void slabcc_model::adjust_extrapolation_params(int &extrapol_steps_num, double &extrapol_steps_size, double &extrapol_grid_x) {
	auto log = spdlog::get("loggers");

	const mat33 cell_vectors0 = cell_vectors;
	const urowvec3 grid0 = cell_grid;

	//Force discretization error checks
	for (auto n = extrapol_steps_num - 1; n > 0; --n) {
		const double extrapol_factor = extrapol_steps_size * n + 1;
		change_size(cell_vectors0 * extrapol_factor);
		gaussian_charges_gen();
		//Adjust the parameters
		if (as_size(grid0)!= as_size(cell_grid)) {
			string adjusted_parameters = "";
			if (type!=model_type::monolayer) {
				if (extrapol_steps_num > 4) {
					extrapol_steps_num = 4;
					adjusted_parameters += " extrapolate_steps_number=" + to_string(extrapol_steps_num);
				}
				if (extrapol_steps_size > 0.25) {
					extrapol_steps_size = 0.25;
					adjusted_parameters += " extrapolate_steps_size=" + to_string(extrapol_steps_size);
				}
			}
			log->debug("Adjusted parameters:{}", adjusted_parameters);
		}
	}
	change_size(cell_vectors0);
	change_grid(grid0);
}

tuple <rowvec, rowvec> slabcc_model::extrapolate(int extrapol_steps_num, double extrapol_steps_size) {


	auto log = spdlog::get("loggers");
	const mat33 cell_vectors0 = cell_vectors;
	const rowvec3 cell_vectors_lengths0 = cell_vectors_lengths;
	const auto interfaces0 = interfaces;
	const auto charge_position0 = charge_position;
	const double slab_thickness = abs(interfaces(0) - interfaces(1));
	rowvec Es = arma::zeros<rowvec>(extrapol_steps_num - 1), sizes = Es;
	for (auto n = 0; n < extrapol_steps_num - 1; ++n) {
		const double extrapol_factor = extrapol_steps_size * (1.0 + n) + 1;
		change_size(cell_vectors0 * extrapol_factor);
		if (this->type == model_type::slab) {
			//increase the slab thickness
			const uvec interface_sorted_i = sort_index(interfaces);
			interfaces(interface_sorted_i(1)) = interfaces(interface_sorted_i(0)) + slab_thickness;
			//move the charges to the same distance from their original nearest interface
			for (uword charge_i = 0; charge_i < charge_position.n_rows; ++charge_i) {
				const rowvec2 initial_distance_to_interfaces = (charge_position0(charge_i, normal_direction) - interfaces0) * cell_vectors_lengths0(normal_direction);
				if (abs(initial_distance_to_interfaces(0)) < abs(initial_distance_to_interfaces(1))) {
					charge_position(charge_i, normal_direction) = interfaces(0) + initial_distance_to_interfaces(0) / cell_vectors_lengths(normal_direction);
				}
				else {
					charge_position(charge_i, normal_direction) = interfaces(1) + initial_distance_to_interfaces(1) / cell_vectors_lengths(normal_direction);;
				}
			}
		}


		gaussian_charges_gen();
		dielectric_profiles_gen();

		// (only works for the orthogonal cells!)
		const auto CHG_normalized = CHG - total_charge / prod(cell_vectors_lengths);
		const auto V = poisson_solver_3D(CHG_normalized, dielectric_profiles, cell_vectors_lengths, normal_direction);
		const auto EperModel = 0.5 * accu(real(V % CHG_normalized)) * voxel_vol * Hartree_to_eV;
		const rowvec2 interface_pos = interfaces * cell_vectors_lengths(normal_direction);
		string extrapolation_info = to_string(extrapol_factor) + "\t" + ::to_string(EperModel) + "\t" + ::to_string(total_charge) + "\t" + to_string(interface_pos);
		for (uword i = 0; i < charge_position.n_rows; ++i) {
			extrapolation_info += "\t" + to_string(charge_position(i, normal_direction) * cell_vectors_lengths(normal_direction));
		}
		log->debug(extrapolation_info);
		Es(n) = EperModel;
		sizes(n) = 1.0 / extrapol_factor;
	}

	return make_tuple(Es, sizes);
}

double slabcc_model::Eiso_bessel() const {
	
	auto logger = spdlog::get("loggers");
	const double K_min = 0.00001;
	const double K_max = 100;
	logger->debug("Isolated energy integration limits in the k-space: {}:{}", K_min, K_max);
	//integration grid
	rowvec K;
	//hand-tuned adaptive grid based on the derivative of the Uk ~~log(k)
	for (int ki = 0; ki < K_max; ++ki) {
		const int grid_density = static_cast<int>(pow(10, -log((ki + 1) / 100.0)) / 5.0) + 1;
		K = join_horiz(K, linspace<rowvec>(ki + K_min, ki + 1, grid_density));
	}
	logger->debug("Number of k-space integration grid points: {}", K.n_elem);

	// eq. 7 in the SI (Supplementary Information for `First-principles electrostatic potentials for reliable alignment at interfaces and defects`)
	const rowvec integrand = K % exp(-square(K) * pow(charge_sigma(0, 0), 2)) % Uk(K);
	const mat integral = trapz(K, integrand, 1);
	const double Q = charge_fraction(0) * total_charge;
	const double U_total = pow(Q, 2) * integral(0) * Hartree_to_eV;

	return U_total;

}

rowvec slabcc_model::Uk(rowvec k) const {
	const double z0 = charge_position(0, normal_direction) * cell_vectors_lengths(normal_direction);
	const rowvec3 length = cell_vectors_lengths;
	const urowvec3 n_points = cell_grid;
	const uword normal = normal_direction;
	const rowvec Gs = 2.0 * PI / length;

	rowvec Gz0 = ceil(regspace<rowvec>(-0.5 * n_points(normal), 0.5 * n_points(normal) - 1)) * Gs(normal);
	Gz0 = ifftshift(Gz0);
	const rowvec Gz02 = square(Gz0);

	const uword LGz = Gz0.n_elem;
	const double dielbulk = 1;				// bulk response far away in vaccum!
	const cx_mat dielsG = fft(dielectric_profiles - dielbulk);
	const cx_mat dielGz = circ_toeplitz(dielsG.col(normal)) / LGz;
	const uword inplane_direction = (normal == 0) ? 1 : 0;
	const cx_mat dielpGz = circ_toeplitz(dielsG.col(inplane_direction)) / LGz;

	const cx_mat Ag1 = dielGz;
	const cx_mat Ag1p = dielpGz;
	const mat Ag2 = Gz0.t() * Gz0;
	const cx_double coef(0, -1); //compiler-specific problem!
	const cx_mat rhok = exp(coef * Gz0 * z0 - Gz02 * pow(charge_sigma(0, 0), 2) / 2.0);
	const cx_mat rhok_t = rhok.t();
	const cx_mat rho = ifft(rhok_t);
	rowvec Uk = zeros(arma::size(k));

	const cx_mat Ag12 = Ag1 % Ag2;
	const auto cosGL_2 = cos(Gz0 * length(normal) / 2.0);
	for (uword i = 0; i < k.n_elem; ++i) {
		const cx_mat Ag = Ag12 + Ag1p * k(i) * k(i);
		const double keff = k(i);
		const mat Kinvg = diagmat(dielbulk * length(normal) * (pow(keff, 2) + Gz02) / (1 - exp(-keff * length(normal) / 2.0) * cosGL_2));
		const cx_mat Dg = Kinvg + length(normal) * Ag;
		const auto VGz = solve(Dg, rhok_t);
		const auto Vz = ifft(VGz) * LGz;
		Uk(i) = real(accu(Vz % rho));
	}

	return Uk;
}

double potential_error(const vector<double>& x, vector<double>& grad, void* model_ptr) {
	auto log = spdlog::get("loggers");

	//input data
	slabcc_model &model = *static_cast<slabcc_model*>(model_ptr);
	model.data_unpacker(x);
	rowvec normalized_charge_fraction = model.charge_fraction;

	//total charge error
	double bounds_factor = 0;

	if (model.charge_fraction(model.charge_fraction.n_elem - 1) < 0) {
		bounds_factor = -model.charge_fraction(model.charge_fraction.n_elem - 1);
		normalized_charge_fraction(normalized_charge_fraction.n_elem - 1) = 0;
		normalized_charge_fraction /= (bounds_factor + 1);
	}

	model.gaussian_charges_gen();
	model.dielectric_profiles_gen();

	model.V = poisson_solver_3D(model.CHG, model.dielectric_profiles, model.cell_vectors_lengths, model.normal_direction);
	model.V_diff = real(model.V) * Hartree_to_eV - model.V_target;
	//bigger output for out-of-bounds input: quadratic penalty
	const double bounds_correction = bounds_factor + 10 * bounds_factor * bounds_factor;
	model.potential_RMSE = sqrt(accu(square(model.V_diff)) / model.V_diff.n_elem) + bounds_correction;

	if (model.initial_potential_RMSE < 0) {
		model.initial_potential_RMSE = model.potential_RMSE;
	}


	log->debug("-----------------------------------------");
	if (!approx_equal(model.diel_in, model.diel_out, "absdiff", 0.02)) {
		const rowvec2 unshifted_interfaces = fmod_p(model.interfaces - model.rounded_relative_shift(model.normal_direction), 1);
		log->debug(" > interfaces={}", ::to_string(unshifted_interfaces));
	}

	const mat unshifted_charge_position = fmod_p(model.charge_position - repmat(model.rounded_relative_shift, model.charge_position.n_rows, 1), 1);

	for (uword i = 0; i < normalized_charge_fraction.n_elem; ++i) {
		log->debug("{}> charge_position={}", i + 1, ::to_string(unshifted_charge_position.row(i)));
		if (model.trivariate_charge) {
			log->debug("{}> charge_sigma={}", i + 1, ::to_string(model.charge_sigma.row(i)));
			if (abs(model.charge_rotations).max() > 0) {
				const rowvec3 rotation = model.charge_rotations.row(i) * 180.0 / PI;
				log->debug("{}> charge_rotation={}", i + 1, ::to_string(rotation));
			}
		}
		else {
			log->debug("{}> charge_sigma={}", i + 1, ::to_string(model.charge_sigma(i, 0)));
		}
		if (model.charge_fraction.n_elem > 1) {
			log->debug("{}> charge_fraction={}", i + 1, model.charge_fraction(i));
		}
	}
	if (bounds_correction > 0) {
		log->debug("Out of the bounds correction to the RMSE: {}", bounds_correction);
	}
	log->debug("Potential Root Mean Square Error: {}", model.potential_RMSE);

	return model.potential_RMSE;
}

void optimize(const string& opt_algo, const double& opt_tol, const int& max_eval, const int& max_time, slabcc_model& model, const opt_switches& optimize) {

	auto log = spdlog::get("loggers");
	auto opt_algorithm = nlopt::LN_COBYLA;
	if (opt_algo == "BOBYQA") {
		opt_algorithm = nlopt::LN_BOBYQA;
	}
	else if (opt_algo == "SBPLX") {
		opt_algorithm = nlopt::LN_SBPLX;
	}
	
	vector<double> opt_param, low_b, upp_b, step_size;
	tie(opt_param, low_b, upp_b, step_size) = model.data_packer(optimize);
	nlopt::opt opt(opt_algorithm, opt_param.size());
	opt.set_lower_bounds(low_b);
	opt.set_upper_bounds(upp_b);
	opt.set_initial_step(step_size);
	opt.set_min_objective(potential_error, &model);
	opt.set_xtol_rel(opt_tol);
	if (max_eval > 0) {
		opt.set_maxeval(max_eval);
	}
	if (max_time > 0) {
		opt.set_maxtime(60.0 * max_time);
	}

	const int sigma_per_charge = model.trivariate_charge ? 3 : 1;
	const int var_per_charge = static_cast<int>(optimize.charge_position) * 3
		+ static_cast<int>(optimize.charge_rotation) * 3
		+ static_cast<int>(optimize.charge_sigma) * sigma_per_charge
		+ static_cast<int>(optimize.charge_fraction) * 1;
	const uword opt_parameters = model.charge_fraction.n_elem * var_per_charge + 2 * optimize.interfaces;
	log->trace("Started optimizing {} model parameters", opt_parameters);
	log->trace("Optimization algorithm: " + string(opt.get_algorithm_name()));
	try {
		const nlopt::result nlopt_final_result = opt.optimize(opt_param, model.potential_RMSE);
		log->debug("-----------------------------------------");
		if (nlopt_final_result == nlopt::MAXEVAL_REACHED) {
			log->warn("Optimization ended after {} steps before reaching the requested accuracy!", max_eval);
		}
		else if (nlopt_final_result == nlopt::MAXTIME_REACHED) {
			log->warn("Optimization ended after {} minutes before reaching the requested accuracy!", max_time);
		}
	}
	catch (const exception & e) {
		log->error("Optimization of the slabcc parameters failed: " + string(e.what()));
		log->error("Please start with better initial guess for the input parameters or use a different optimization algorithm.");
	}

	model.data_unpacker(opt_param);

	log->trace("Optimization ended.");
}
