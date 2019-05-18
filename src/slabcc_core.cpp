// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_core.hpp"

extern slabcc_cell model_cell;

mat dielectric_profiles_gen(const rowvec2 &interfaces, const rowvec3 &diel_in, const rowvec3 &diel_out, const double &diel_erf_beta) {
	const auto length = model_cell.vec_lengths(model_cell.normal_direction);
	const auto n_points = model_cell.grid(model_cell.normal_direction);
	rowvec2 interfaces_cartesian = interfaces * length;
	interfaces_cartesian = sort(interfaces_cartesian);
	const auto positions = linspace<rowvec>(0, length, n_points + 1);
	mat dielectric_profiles = zeros(n_points, 3);
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

	return dielectric_profiles;
}

cx_cube gaussian_charge(const double& Q, const rowvec3& rel_pos, const rowvec3& sigma, const rowvec3& rotation_angle, const bool& trivariate) {

	rowvec x0 = linspace<rowvec>(0, model_cell.vec_lengths(0) - model_cell.vec_lengths(0) / model_cell.grid(0), model_cell.grid(0));
	rowvec y0 = linspace<rowvec>(0, model_cell.vec_lengths(1) - model_cell.vec_lengths(1) / model_cell.grid(1), model_cell.grid(1));
	rowvec z0 = linspace<rowvec>(0, model_cell.vec_lengths(2) - model_cell.vec_lengths(2) / model_cell.grid(2), model_cell.grid(2));

	// shift the axis reference to position of the Gaussian charge center
	x0 -= accu(model_cell.vectors.col(0) * rel_pos(0));
	y0 -= accu(model_cell.vectors.col(1) * rel_pos(1));
	z0 -= accu(model_cell.vectors.col(2) * rel_pos(2));

	//handle the minimum distance from the mirror charges
	for (auto &pos : x0) {
		if (abs(pos) > model_cell.vec_lengths(0) / 2) {
			pos = model_cell.vec_lengths(0) - abs(pos);
		}
	}
	for (auto &pos : y0) {
		if (abs(pos) > model_cell.vec_lengths(1) / 2) {
			pos = model_cell.vec_lengths(1) - abs(pos);
		}
	}
	for (auto &pos : z0) {
		if (abs(pos) > model_cell.vec_lengths(2) / 2) {
			pos = model_cell.vec_lengths(2) - abs(pos);
		}
	}

	cube x, y, z;
	tie(x, y, z) = ndgrid(x0, y0, z0);

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
		for (uword i = 0; i < x.n_elem; ++i) {
			const vec3 old_coordinates = { x(i), y(i), z(i) };
			const vec3 new_coordinates = rotation_mat * old_coordinates;
			x(i) = new_coordinates(0);
			y(i) = new_coordinates(1);
			z(i) = new_coordinates(2);
		}
	}

	const cube r2 = square(x) + square(y) + square(z);

	// this charge distribution is due to the 1st nearest gaussian image. 
	// In case of the very small supercells or very diffuse charges (large sigma), the higher order of the image charges must also be included.
	// But the validity of the correction method for these cases must be checked!	
	cx_cube charge_dist;

	if (trivariate) {
		charge_dist = cx_cube(Q / (pow(2 * PI, 1.5) * prod(sigma))
			* exp(-square(x) / (2 * square(sigma(0))) - square(y) / (2 * square(sigma(1))) - square(z) / (2 * square(sigma(2))))
			, zeros(as_size(model_cell.grid)));
	}
	else {
		charge_dist = cx_cube(Q / pow((sigma(0) * sqrt(2 * PI)), 3) * exp(-r2 / (2 * square(sigma(0)))), zeros(as_size(model_cell.grid)));
	}
	return charge_dist;
}

cx_cube poisson_solver_3D(const cx_cube &rho, mat diel) {
	auto length = model_cell.vec_lengths;
	auto n_points = model_cell.grid;

	if (model_cell.normal_direction != 2) {
		n_points.swap_cols(model_cell.normal_direction, 2);
		length.swap_cols(model_cell.normal_direction, 2);
		diel.swap_cols(model_cell.normal_direction, 2);
	}

	const rowvec Gs = 2.0 * PI / length;

	rowvec Gx0 = ceil(regspace<rowvec>(-0.5 * n_points(0), 0.5 * n_points(0) - 1)) * Gs(0);
	rowvec Gy0 = ceil(regspace<rowvec>(-0.5 * n_points(1), 0.5 * n_points(1) - 1)) * Gs(1);
	rowvec Gz0 = ceil(regspace<rowvec>(-0.5 * n_points(2), 0.5 * n_points(2) - 1)) * Gs(2);

	Gx0 = ifftshift(Gx0);
	Gy0 = ifftshift(Gy0);
	Gz0 = ifftshift(Gz0);

	// 4PI is for the atomic units
	const auto rhok = fft(cx_cube(4.0 * PI * rho));
	const cx_mat dielsG = fft(diel);
	const cx_mat eps11 = circ_toeplitz(dielsG.col(0)) / Gz0.n_elem;
	const cx_mat eps22 = circ_toeplitz(dielsG.col(1)) / Gz0.n_elem;
	const cx_mat eps33 = circ_toeplitz(dielsG.col(2)) / Gz0.n_elem;
	const mat GzGzp = Gz0.t() * Gz0;
	const cx_mat Az = eps33 % GzGzp;
	cx_cube Vk(arma::size(rhok));

	#pragma omp parallel for firstprivate(Az,eps11,eps22,rhok)
	for (uword k = 0; k < Gx0.n_elem; ++k) {
		const cx_mat eps11_Gx0k2 = eps11 * square(Gx0(k));
		for (uword m = 0; m < Gy0.n_elem; ++m) {
			vector<span> spans = { span(k), span(m), span() };
			swap(spans[model_cell.normal_direction], spans[2]);
			cx_mat AG = Az + eps11_Gx0k2 + eps22 * square(Gy0(m));
			if ((k == 0) && (m == 0)) { AG(0, 0) = 1; }
			Vk(spans[0], spans[1], spans[2]) = solve(AG, vectorise(rhok(spans[0], spans[1], spans[2])));
		}
	}
	// 0,0,0 in k-space corresponds to a constant in the real space: average potential over the supercell.
	Vk(0, 0, 0) = 0;
	const cx_cube V = ifft(Vk);

	return V;
}

double potential_eval(const vector<double> &x, vector<double> &grad, void *slabcc_data) {

	rowvec2 shifted_interfaces;
	mat charge_sigma;
	mat charge_rotations;
	rowvec input_charge_fraction;
	mat charge_position;
	opt_variable variables = { shifted_interfaces, charge_sigma, charge_rotations, input_charge_fraction, charge_position };
	optimizer_unpacker(x, variables);

	//input data
	const auto d = static_cast<opt_data *>(slabcc_data);
	const rowvec3 &diel_in = d->diel_in;
	const rowvec3 &diel_out = d->diel_out;
	const auto &diel_erf_beta = d->diel_erf_beta;
	const cube &defect_potential = d->defect_potential;
	const double &total_vasp_charge = d->total_vasp_charge;
	const bool &trivariate = d->trivariate;
	const rowvec3 &rounded_relative_shift = d->rounded_relative_shift;
	
	rowvec normalized_charge_fraction = input_charge_fraction;
	
	//total charge error
	double bounds_factor = 0;

	if(input_charge_fraction(input_charge_fraction.n_elem - 1) < 0){
		bounds_factor = -input_charge_fraction(input_charge_fraction.n_elem - 1);
		normalized_charge_fraction(normalized_charge_fraction.n_elem - 1) = 0;
		normalized_charge_fraction /= (bounds_factor + 1);
	}

	//output data
	mat& dielectric_profiles = d->dielectric_profiles;
	cx_cube& rhoM = d->rhoM;
	cx_cube& V = d->V;
	cube& V_diff = d->V_diff;
	auto& initial_potential_RMSE = d->initial_potential_RMSE;

	dielectric_profiles = dielectric_profiles_gen(shifted_interfaces, diel_in, diel_out, diel_erf_beta);

	rhoM.reset();
	rhoM = zeros<cx_cube>(as_size(model_cell.grid));
	for (uword i = 0; i < normalized_charge_fraction.n_elem; ++i) {
		rhoM += gaussian_charge(normalized_charge_fraction(i) * total_vasp_charge, charge_position.row(i), charge_sigma.row(i), charge_rotations.row(i), trivariate);
	}

	V = poisson_solver_3D(rhoM, dielectric_profiles);
	V_diff = real(V) * Hartree_to_eV - defect_potential;
	//bigger output for out-of-bounds input: quadratic penalty
	const double bounds_correction = bounds_factor + 10 * bounds_factor * bounds_factor;
	const double potential_RMSE = sqrt(accu(square(V_diff)) / V_diff.n_elem) + bounds_correction;

	if (initial_potential_RMSE < 0) {
		initial_potential_RMSE = potential_RMSE;
	}
	
	auto log = spdlog::get("loggers");

	log->debug("-----------------------------------------");
	if (!approx_equal(diel_in, diel_out, "absdiff", 0.02)) {
		const rowvec2 unshifted_interfaces = fmod_p(shifted_interfaces - rounded_relative_shift(model_cell.normal_direction), 1);
		log->debug(" > interfaces={}", ::to_string(unshifted_interfaces));
	}

	const mat unshifted_charge_position = fmod_p(charge_position - repmat(rounded_relative_shift, charge_position.n_rows, 1), 1);

	for (uword i = 0; i < normalized_charge_fraction.n_elem; ++i) {
		log->debug("{}> charge_position={}", i + 1, ::to_string(unshifted_charge_position.row(i)));
		if (trivariate) {
			log->debug("{}> charge_sigma={}", i + 1, ::to_string(charge_sigma.row(i)));
			if (abs(charge_rotations).max() > 0) {
				const rowvec3 rotation = charge_rotations.row(i) * 180.0 / PI;
				log->debug("{}> charge_rotation={}", i + 1, ::to_string(rotation));
			}
		}
		else {
			log->debug("{}> charge_sigma={}", i + 1, ::to_string(charge_sigma(i,0)));
		}
		if (input_charge_fraction.n_elem > 1) {
			log->debug("{}> charge_fraction={}", i + 1, input_charge_fraction(i));
		}		
	}
	if (bounds_correction > 0) {
		log->debug("Out of the bounds correction to the RMSE: {}", bounds_correction);
	}
	log->debug("Potential Root Mean Square Error: {}", potential_RMSE);

	return potential_RMSE;
}

double do_optimize(const string& opt_algo, const double& opt_tol, const int &max_eval, const int &max_time, opt_data& opt_data, opt_variable& opt_vars, opt_switch& optimize) {
	auto log = spdlog::get("loggers");
	double pot_MSE_min = 0;
	auto opt_algorithm = nlopt::LN_COBYLA;
	if (opt_algo == "BOBYQA") {
		opt_algorithm = nlopt::LN_BOBYQA;
	}
	else if (opt_algo == "SBPLX") {
		opt_algorithm = nlopt::LN_SBPLX;
	}
	
	vector<double> opt_param, low_b, upp_b, step_size;
	tie(opt_param, low_b, upp_b, step_size) = optimizer_packer(opt_vars, optimize, opt_data.trivariate);
	nlopt::opt opt(opt_algorithm, opt_param.size());
	opt.set_lower_bounds(low_b);
	opt.set_upper_bounds(upp_b);
	opt.set_initial_step(step_size);

	opt.set_min_objective(potential_eval, &opt_data);
	opt.set_xtol_rel(opt_tol); 
	if (max_eval > 0) {
		opt.set_maxeval(max_eval);
	}
	if (max_time > 0) {
		opt.set_maxtime(60.0 * max_time);
	}
	//this seems to improve the COBYLA's performance is some cases
	if (opt_vars.charge_fraction.n_elem > 2) {
		if (opt_algo == "COBYLA") {
			opt.add_inequality_constraint(opt_charge_constraint, &opt_data, 1e-6);
		}
	}

	const int sigma_per_charge = opt_data.trivariate ? 3 : 1;
	const int var_per_charge = static_cast<int>(optimize.charge_position) * 3
		+ static_cast<int>(optimize.charge_rotation) * 3
		+ static_cast<int>(optimize.charge_sigma) * sigma_per_charge
		+ static_cast<int>(optimize.charge_fraction) * 1;
	const uword opt_parameters = opt_vars.charge_fraction.n_elem * var_per_charge + 2 * optimize.interfaces;
	log->trace("Started optimizing {} model parameters", opt_parameters);
	log->trace("Optimization algorithm: " + string(opt.get_algorithm_name()));
	try {
		const nlopt::result nlopt_final_result = opt.optimize(opt_param, pot_MSE_min);
		log->debug("-----------------------------------------");
		if (nlopt_final_result == nlopt::MAXEVAL_REACHED) {
			log->warn("Optimization ended after {} steps before reaching the requested accuracy!", max_eval);
		}
		else if (nlopt_final_result == nlopt::MAXTIME_REACHED) {
			log->warn("Optimization ended after {} minutes before reaching the requested accuracy!", max_time);
		}
	}
	catch (const exception &e) {
		log->error("Optimization of the slabcc parameters failed: " + string(e.what()));
		log->error("Please start with better initial guess for the input parameters or use a different optimization algorithm.");
	}

	optimizer_unpacker(opt_param, opt_vars);
	log->trace("Optimization ended.");

	return pot_MSE_min;
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>> optimizer_packer(const opt_variable& opt_vars, opt_switch optimize, const bool trivariate) {
	auto log = spdlog::get("loggers");
	//size of the first step for each parameter
	const double move_step = 1; //Ang
	const double sigma_step = 0.5;
	const double fraction_step = 0.2;
	const double rotation_step = 0.2; //rad
	const rowvec3 relative_move_step = move_step * ang_to_bohr / model_cell.vec_lengths;
	
	//------interfaces----
	vector<double> step_size = { relative_move_step(model_cell.normal_direction) , relative_move_step(model_cell.normal_direction) };
	vector<double> opt_param = { opt_vars.interfaces(0), opt_vars.interfaces(1) };
	vector<double> low_b = { 0, 0 };		//lower bounds
	vector<double> upp_b = { 1, 1 };		//upper bounds
	if (!optimize.interfaces) {
		low_b = opt_param;
		upp_b = opt_param;
	}


	//------charge positions----
	for (uword i = 0; i < opt_vars.charge_position.n_rows; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_param.push_back(opt_vars.charge_position(i, j));
			if (optimize.charge_position) {
				low_b.push_back(0);
				upp_b.push_back(1);
			}
			else {
				low_b.push_back(opt_vars.charge_position(i, j));
				upp_b.push_back(opt_vars.charge_position(i, j));
			}
		}
		step_size.insert(step_size.end(), { relative_move_step(0), relative_move_step(1), relative_move_step(2) });

		//------charge_sigma----
		opt_param.insert(opt_param.end(), { opt_vars.charge_sigma(i, 0), opt_vars.charge_sigma(i, 1), opt_vars.charge_sigma(i, 2) });
		step_size.insert(step_size.end(), { sigma_step,sigma_step,sigma_step });

		if (optimize.charge_sigma) {
			if (trivariate) {
				low_b.insert(low_b.end(), { 0.1, 0.1, 0.1 });
				upp_b.insert(upp_b.end(), { 7, 7, 7 });
			}
			else {
				low_b.insert(low_b.end(), { 0.1, opt_vars.charge_sigma(i, 1), opt_vars.charge_sigma(i, 2) });
				upp_b.insert(upp_b.end(), { 7, opt_vars.charge_sigma(i, 1), opt_vars.charge_sigma(i, 2) });
			}

		}
		else {
			low_b.insert(low_b.end(), { opt_vars.charge_sigma(i, 0), opt_vars.charge_sigma(i, 1), opt_vars.charge_sigma(i, 2) });
			upp_b.insert(upp_b.end(), { opt_vars.charge_sigma(i, 0), opt_vars.charge_sigma(i, 1), opt_vars.charge_sigma(i, 2) });
		}

		//-----charge_rotations-----
		opt_param.insert(opt_param.end(), { opt_vars.charge_rotations(i, 0), opt_vars.charge_rotations(i, 1), opt_vars.charge_rotations(i, 2) });
		step_size.insert(step_size.end(), { rotation_step,rotation_step,rotation_step });

		if (optimize.charge_rotation) {
			if (trivariate) {
				low_b.insert(low_b.end(), { -0.5 * PI , -0.5 * PI , -0.5 * PI });
				upp_b.insert(upp_b.end(), { 0.5 * PI , 0.5 * PI , 0.5 * PI });
			}
			else {
				//this must be caught in the checker!
				log->warn("Optimizing the rotation angles for the simple Gaussian charges is not possible!");
				low_b.insert(low_b.end(), { opt_vars.charge_rotations(i, 0), opt_vars.charge_rotations(i, 1), opt_vars.charge_rotations(i, 2) });
				upp_b.insert(upp_b.end(), { opt_vars.charge_rotations(i, 0), opt_vars.charge_rotations(i, 1), opt_vars.charge_rotations(i, 2) });
			}
		}
		else {
			low_b.insert(low_b.end(), { opt_vars.charge_rotations(i, 0), opt_vars.charge_rotations(i, 1), opt_vars.charge_rotations(i, 2) });
			upp_b.insert(upp_b.end(), { opt_vars.charge_rotations(i, 0), opt_vars.charge_rotations(i, 1), opt_vars.charge_rotations(i, 2) });
		}

		//------charge fractions----
		opt_param.push_back(opt_vars.charge_fraction(i));
		step_size.push_back(fraction_step);

		if (optimize.charge_fraction) {
			low_b.push_back(0);
			upp_b.push_back(1);
		}
		else {
			low_b.push_back(opt_vars.charge_fraction(i));
			upp_b.push_back(opt_vars.charge_fraction(i));
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

void optimizer_unpacker(const vector<double> &optimizer_vars_vec, opt_variable &opt_vars) {

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
	opt_vars.charge_sigma = zeros(n_charges, 3);
	opt_vars.charge_rotations = zeros(n_charges, 3);
	opt_vars.charge_fraction = zeros<rowvec>(n_charges);
	opt_vars.charge_position = zeros(n_charges, 3);
	opt_vars.interfaces = { optimizer_vars_vec.at(0), optimizer_vars_vec.at(1) };

	for (uword i = 0; i < n_charges; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_vars.charge_position(i, j) = optimizer_vars_vec.at(position_offset + i * variables_per_charge + j);
			opt_vars.charge_sigma(i, j) = optimizer_vars_vec.at(sigma_offset + i * variables_per_charge + j);
			opt_vars.charge_rotations(i, j) = optimizer_vars_vec.at(rotation_offset + i * variables_per_charge + j);
		}
		if (i != n_charges - 1) {
			opt_vars.charge_fraction(i) = optimizer_vars_vec.at(fraction_offset + i * variables_per_charge);
		}
		else {
			opt_vars.charge_fraction(i) = 1.0 - accu(opt_vars.charge_fraction);
		}
	}
}



double opt_charge_constraint(const vector<double> &x, vector<double> &grad, void *data)
{
	auto log = spdlog::get("loggers");
	rowvec2 interfaces;
	rowvec charge_fraction;
	mat charge_sigma, charge_rotations, charge_position;
	opt_variable variables = { interfaces, charge_sigma, charge_rotations, charge_fraction, charge_position };
	optimizer_unpacker(x, variables);

	const auto d = static_cast<const opt_data *>(data);
	const auto total_vasp_charge = d->total_vasp_charge;
	const auto constraint = -charge_fraction(charge_fraction.n_elem - 1);
	const rowvec charge_q = charge_fraction * total_vasp_charge;
	log->debug("Charge in each Gaussian:" + ::to_string(charge_q));
	return constraint;
}

void verify_cells(const supercell& Neutral_supercell, const supercell& Charged_supercell) {

	auto log = spdlog::get("loggers");

	// equal size
	if (!approx_equal(Neutral_supercell.cell_vectors * Neutral_supercell.scaling,
		Charged_supercell.cell_vectors * Charged_supercell.scaling, "reldiff", 0.001)) {
		log->debug("Neutral supercell vectors: " + to_string(Neutral_supercell.cell_vectors * Neutral_supercell.scaling));
		log->debug("Charged supercell vectors: " + to_string(Charged_supercell.cell_vectors * Charged_supercell.scaling));
		log->critical("Cell vectors of the input files does not match!");
		finalize_loggers();
		exit(1);
	}

	//orthogonal
	auto cell = Neutral_supercell.cell_vectors;
	cell.each_col([](vec& vec) { vec /= norm(vec); });
	if (!approx_equal(cell.t() * cell, eye(3, 3), "reldiff", 0.001)) {
		log->debug("Cell vectors: " + to_string(Neutral_supercell.cell_vectors));
		log->debug("Cell vectors basis: " + to_string(cell));
		log->debug("Orthogonality criteria: " + to_string(cell.t() * cell));
		log->critical("Supercell vectors are not orthogonal!");
		finalize_loggers();
		exit(1);
	}

	// equal grid
	const SizeCube input_grid = arma::size(Neutral_supercell.potential);
	if ((input_grid != arma::size(Charged_supercell.potential)) ||
		(input_grid != arma::size(Charged_supercell.charge)) ||
		(input_grid != arma::size(Neutral_supercell.charge))) {
		log->debug("Neutral CHGCAR grid: " + to_string(arma::size(Neutral_supercell.charge)));
		log->debug("Neutral LOCPOT grid: " + to_string(arma::size(Neutral_supercell.potential)));
		log->debug("Charged CHGCAR grid: " + to_string(arma::size(Charged_supercell.charge)));
		log->debug("Charged LOCPOT grid: " + to_string(arma::size(Charged_supercell.potential)));
		log->critical("Data grids of the CHGCAR or LOCPOT files does not match!");
		finalize_loggers();
		exit(1);
	}

	log->trace("All files are loaded and cross-checked!");
}

void verify_interface_optimization(const rowvec2& initial_interfaces, const rowvec2& optimized_interfaces) {
	
	const double max_total_safe_movement = 4; //Ang
	const double normal_length = model_cell.vec_lengths(model_cell.normal_direction) / ang_to_bohr;
	auto log = spdlog::get("loggers");

	const rowvec3 mirroring_template = { -1,0,1 };
	const mat interface_all_mirrors = repmat(mirroring_template, 2, 1) + repmat(optimized_interfaces.t(), 1, 3);
	const mat interface_n0_shiftes = abs(interface_all_mirrors - initial_interfaces(0));
	vec interface_n0_distances = min(interface_n0_shiftes, 1);
	const mat interface_n1_shiftes = abs(interface_all_mirrors - initial_interfaces(1));
	vec interface_n1_distances = min(interface_n1_shiftes, 1);
	const rowvec2 total_changes = { interface_n0_distances(0) + interface_n1_distances(1),  interface_n0_distances(1) + interface_n1_distances(0) };
	
	//convert to Ang
	interface_n0_distances *= normal_length;
	interface_n1_distances *= normal_length;

	if (total_changes(0) < total_changes(1)) {
		log->trace("The order of the interfaces has not been swapped!");
		log->trace("1st interface has moved {} angstrom from {} to {}", interface_n0_distances(0), initial_interfaces(0), optimized_interfaces(0));
		log->trace("2nd interface has moved {} angstrom from {} to {}", interface_n1_distances(1), initial_interfaces(1), optimized_interfaces(1));
	}
	else {
		log->trace("The order of the interfaces has been swapped!");
		log->trace("1st interface has moved {} angstrom from {} to {}", interface_n0_distances(1), initial_interfaces(0), optimized_interfaces(1));
		log->trace("2nd interface has moved {} angstrom from {} to {}", interface_n1_distances(0), initial_interfaces(1), optimized_interfaces(0));
	}

	if (total_changes(0) * normal_length > max_total_safe_movement) {
		log->warn("There is a relatively large change in the position of the interfaces after the optimization! Please ensure the correctness of the interfaces.");
	}

}


void verify_charge_sigma_optimization(const rowvec& charge_q, const mat& charge_sigma) {
	//maximum sigma for the localized charges
	const double max_sigma = 6.5;
	const double min_non_negligable_q = 0.01;
	auto log = spdlog::get("loggers");
	for (uword i = 0; i < charge_q.n_elem; ++i) {
		if (abs(charge_q(i)) > min_non_negligable_q) {
			if (charge_sigma.row(i).max() > max_sigma) {
				log->error("The model extra charge seems to be (at least partially) delocalized. The current charge correction method is not suitable for the delocalized charges.");
			}
		}
	}
}