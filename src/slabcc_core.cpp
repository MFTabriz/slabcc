// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_core.hpp"

extern slabcc_cell_type slabcc_cell;

mat dielectric_profiles_gen(const rowvec2 &interfaces, const rowvec3 &diel_in, const rowvec3 &diel_out, const double &diel_erf_beta) {
	const auto length = slabcc_cell.vec_lengths(slabcc_cell.normal_direction);
	const auto n_points = slabcc_cell.grid(slabcc_cell.normal_direction);
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

void UpdateCell(const mat33& vectors, const urowvec3& grid) {
	slabcc_cell.vectors = vectors;
	for (uword i = 0; i < 3; ++i) {
		slabcc_cell.vec_lengths(i) = norm(vectors.col(i));
	}
	slabcc_cell.grid = grid;
	slabcc_cell.voxel_vol = prod(slabcc_cell.vec_lengths / grid);
}

cx_cube gaussian_charge(const double& Q, const vec3& rel_pos, const rowvec3& sigma, const bool& trivariate) {

	rowvec x0 = linspace<rowvec>(0, slabcc_cell.vec_lengths(0) - slabcc_cell.vec_lengths(0) / slabcc_cell.grid(0), slabcc_cell.grid(0));
	rowvec y0 = linspace<rowvec>(0, slabcc_cell.vec_lengths(1) - slabcc_cell.vec_lengths(1) / slabcc_cell.grid(1), slabcc_cell.grid(1));
	rowvec z0 = linspace<rowvec>(0, slabcc_cell.vec_lengths(2) - slabcc_cell.vec_lengths(2) / slabcc_cell.grid(2), slabcc_cell.grid(2));

	// shift the axis reference to position of the Gaussian charge center
	x0 -= accu(slabcc_cell.vectors.col(0) * rel_pos(0));
	y0 -= accu(slabcc_cell.vectors.col(1) * rel_pos(1));
	z0 -= accu(slabcc_cell.vectors.col(2) * rel_pos(2));

	//handle the minimum distance from the mirror charges
	for (auto &pos : x0) {
		if (abs(pos) > slabcc_cell.vec_lengths(0) / 2) {
			pos = slabcc_cell.vec_lengths(0) - abs(pos);
		}
	}
	for (auto &pos : y0) {
		if (abs(pos) > slabcc_cell.vec_lengths(1) / 2) {
			pos = slabcc_cell.vec_lengths(1) - abs(pos);
		}
	}
	for (auto &pos : z0) {
		if (abs(pos) > slabcc_cell.vec_lengths(2) / 2) {
			pos = slabcc_cell.vec_lengths(2) - abs(pos);
		}
	}

	cube x, y, z;
	tie(x, y, z) = ndgrid(x0, y0, z0);

	const cube r2 = square(x) + square(y) + square(z);
	// this charge distribution is due to the 1st nearest gaussian image. 
	// In case of the very small supercells or very diffuse charges (large sigma), the higher order of the image charges must also be included.
	// But the validity of the correction method for these cases must be checked!
	
	cx_cube charge_dist;

	if (trivariate) {
		charge_dist = cx_cube(Q / (pow(2 * PI, 1.5) * prod(sigma))
			* exp(-square(x) / (2 * square(sigma(0))) - square(y) / (2 * square(sigma(1))) - square(z) / (2 * square(sigma(2))))
			, zeros(as_size(slabcc_cell.grid)));
	}
	else {
		charge_dist = cx_cube(Q / pow((sigma(0) * sqrt(2 * PI)), 3) * exp(-r2 / (2 * square(sigma(0)))), zeros(as_size(slabcc_cell.grid)));
	}
	return charge_dist;
}

cx_cube poisson_solver_3D(const cx_cube &rho, mat diel) {
	auto length = slabcc_cell.vec_lengths;
	auto n_points = slabcc_cell.grid;

	if (slabcc_cell.normal_direction != 2) {
		n_points.swap_cols(slabcc_cell.normal_direction, 2);
		length.swap_cols(slabcc_cell.normal_direction, 2);
		diel.swap_cols(slabcc_cell.normal_direction, 2);
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
			swap(spans[slabcc_cell.normal_direction], spans[2]);
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
	rowvec charge_fraction;
	mat charge_position;
	opt_vars variables = { shifted_interfaces, charge_sigma, charge_fraction, charge_position };
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

	//rest of the charge goes to the last Gaussian
	charge_fraction(charge_fraction.n_elem - 1) = 1 - accu(charge_fraction);

	//output data
	mat& dielectric_profiles = d->dielectric_profiles;
	cx_cube& rhoM = d->rhoM;
	cx_cube& V = d->V;
	cube& V_diff = d->V_diff;
	auto& initial_potential_MSE = d->initial_potential_MSE;

	dielectric_profiles = dielectric_profiles_gen(shifted_interfaces, diel_in, diel_out, diel_erf_beta);

	rhoM.reset();
	rhoM = zeros<cx_cube>(as_size(slabcc_cell.grid));
	for (uword i = 0; i < charge_fraction.n_elem; ++i) {
		rhoM += gaussian_charge(charge_fraction(i) * total_vasp_charge, charge_position.row(i).t(), charge_sigma.row(i), trivariate);
	}

	V = poisson_solver_3D(rhoM, dielectric_profiles);
	V_diff = real(V) * Hartree_to_eV - defect_potential;

	const auto potential_MSE = accu(square(V_diff)) / V_diff.n_elem;

	if (initial_potential_MSE < 0) {
		initial_potential_MSE = potential_MSE;
	}

	auto log = spdlog::get("loggers");

	log->debug("-----------------------------------------");
	if (!approx_equal(diel_in, diel_out, "absdiff", 0.02)) {
		const rowvec2 unshifted_interfaces = fmod_p(shifted_interfaces - rounded_relative_shift(slabcc_cell.normal_direction), 1);
		log->debug("> interfaces=" + ::to_string(unshifted_interfaces));
	}

	const mat unshifted_charge_position = fmod_p(charge_position - repmat(rounded_relative_shift, charge_position.n_rows, 1), 1);

	for (uword i = 0; i < charge_fraction.n_elem; ++i) {
		log->debug("> charge_position=" + ::to_string(unshifted_charge_position.row(i)));
		if (trivariate) {
			log->debug("> charge_sigma=" + ::to_string(charge_sigma.row(i)));
		}
		else {
			log->debug("> charge_sigma=" + ::to_string(charge_sigma(i,0)));
		}
		if (charge_fraction.n_elem > 1) {
			log->debug("> charge_fraction={}", charge_fraction(i));
		}		
	}
	log->debug("Potential Root Mean Square Error: {}", sqrt(potential_MSE));

	return potential_MSE;
}

double do_optimize(const string& opt_algo, const double& opt_tol, const int &max_eval, const int &max_time, opt_data& opt_data, opt_vars& opt_vars, const bool &optimize_charge_position, const bool &optimize_charge_sigma, const bool &optimize_charge_fraction, const bool &optimize_interfaces) {
	auto log = spdlog::get("loggers");
	double pot_MSE_min = 0;
	auto opt_algorithm = nlopt::LN_COBYLA;
	
	if (opt_vars.charge_fraction.n_elem < 3) {
		if (opt_algo == "BOBYQA") {
			opt_algorithm = nlopt::LN_BOBYQA;
		} 
		else if (opt_algo == "SBPLX") {
			opt_algorithm = nlopt::LN_SBPLX;
		}
	}
	else {
		if (opt_algo != "COBYLA") {
			//SBPLX/BOBYQA in NLOPT 2.4.2 does not support the inequality constraints which is needed to preserve the total charge of the model!
			log->warn("SBPLX/BOBYQA optimization algorithms does not support the models with multiple charges! COBYLA will be used instead!");
		}
	}

	vector<double> opt_param, low_b, upp_b, step_size;
	tie(opt_param, low_b, upp_b, step_size) = optimizer_packer(opt_vars, optimize_charge_position, optimize_charge_sigma, optimize_charge_fraction, optimize_interfaces, opt_data.trivariate);
	nlopt::opt opt(opt_algorithm, opt_param.size());
	opt.set_lower_bounds(low_b);
	opt.set_upper_bounds(upp_b);
	opt.set_initial_step(step_size);

	opt.set_min_objective(potential_eval, &opt_data);
	opt.set_xtol_rel(square(opt_tol));   //tolerance is defined for RMSE in the input, optimizer is checking MSE
	if (max_eval > 0) {
		opt.set_maxeval(max_eval);
	}
	if (max_time > 0) {
		opt.set_maxtime(60.0 * max_time);
	}
	if (opt_vars.charge_fraction.n_elem > 2) {
		//add constraint to keep the total charge constant
		opt.add_inequality_constraint(opt_charge_constraint, &opt_data, 1e-6);
	}

	const int sigma_per_charge = opt_data.trivariate ? 3 : 1;
	const int var_per_charge = static_cast<int>(optimize_charge_position) * 3
		+ static_cast<int>(optimize_charge_sigma) * sigma_per_charge
		+ static_cast<int>(optimize_charge_fraction) * 1;
	const uword opt_parameters = opt_vars.charge_fraction.n_elem * var_per_charge + 2 * optimize_interfaces;
	log->trace("Started optimizing {} model parameters", opt_parameters);
	log->trace("Optimization algorithm: " + string(opt.get_algorithm_name()));
	try {
		const nlopt::result nlopt_final_result = opt.optimize(opt_param, pot_MSE_min);
		log->debug("-----------------------------------------");
		if (nlopt_final_result == nlopt::MAXEVAL_REACHED) {
			log->warn("Optimization ended after {} steps before reaching the requested accuracy!", max_eval);
		}
		else if (nlopt_final_result == nlopt::MAXTIME_REACHED) {
			log->warn("Optimization ended after {} minutes of search before reaching the requested accuracy!", max_time);
		}
	}
	catch (const exception &e) {
		log->error("Optimization of the slabcc parameters failed: " + string(e.what()));
		log->error("Please start with better initial guess for the input parameters or change the optimization algorithm.");
	}

	optimizer_unpacker(opt_param, opt_vars);
	log->trace("Optimization ended.");

	return pot_MSE_min;
}

tuple<vector<double>, vector<double>, vector<double>, vector<double>> optimizer_packer(const opt_vars& opt_vars, const bool optimize_charge_position, const bool optimize_charge_sigma, const bool optimize_charge_fraction, const bool optimize_interface, const bool trivariate) {
	//size of the first step for each parameter
	const double move_step = 1; //Ang
	const double sigma_step = 0.5;
	const double fraction_step = 0.2;
	const rowvec3 relative_move_step = move_step * ang_to_bohr / slabcc_cell.vec_lengths;
	
	//------interfaces----
	vector<double> step_size = { relative_move_step(slabcc_cell.normal_direction) , relative_move_step(slabcc_cell.normal_direction) };
	vector<double> opt_param = { opt_vars.interfaces(0), opt_vars.interfaces(1) };
	vector<double> low_b = { 0, 0 };		//lower bounds
	vector<double> upp_b = { 1, 1 };		//upper bounds
	if (!optimize_interface) {
		low_b = opt_param;
		upp_b = opt_param;
	}


	//------charge positions----
	for (uword i = 0; i < opt_vars.charge_position.n_rows; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_param.push_back(opt_vars.charge_position(i, j));
			if (optimize_charge_position) {
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

		if (optimize_charge_sigma) {
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

		//------charge fractions----
		opt_param.push_back(opt_vars.charge_fraction(i));
		step_size.push_back(fraction_step);

		if (optimize_charge_fraction) {
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

void optimizer_unpacker(const vector<double> &optimizer_vars_vec, opt_vars &opt_vars) {
	// the coefficients in this part depend on the set of our variables
	// they are ordered as: interfaces, :|[x, y, z, 3x sigma, q]|:
	// NOTE: the last "q" must be calculated from the total charge (which in not available here!)
	// NOTE2: this function is also called by opt_charge_constraint() with a reference to empty variable.
	const uword variable_per_charge = 7;
	const uword n_charges = optimizer_vars_vec.size() / variable_per_charge;
	opt_vars.charge_sigma = zeros(n_charges, 3);
	opt_vars.charge_fraction = zeros<rowvec>(n_charges);
	opt_vars.charge_position = zeros(n_charges, 3);
	opt_vars.interfaces = { optimizer_vars_vec.at(0), optimizer_vars_vec.at(1) };

	for (uword i = 0; i < n_charges; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_vars.charge_position(i, j) = optimizer_vars_vec.at(2 + i * variable_per_charge + j);
			opt_vars.charge_sigma(i, j) = optimizer_vars_vec.at(5 + i * variable_per_charge + j);
		}
		if (i != n_charges - 1) {
			opt_vars.charge_fraction(i) = optimizer_vars_vec.at(8 + i * variable_per_charge);
		}
	}
}


void check_inputs(input_data input_set) {
	auto log = spdlog::get("loggers");
	input_set.charge_sigma = abs(input_set.charge_sigma);
	input_set.max_eval = abs(input_set.max_eval);
	input_set.max_time = abs(input_set.max_time);
	input_set.interfaces = fmod_p(input_set.interfaces, 1);
	input_set.extrapol_grid_x = abs(input_set.extrapol_grid_x);
	input_set.opt_grid_x = abs(input_set.opt_grid_x);
	input_set.opt_tol = abs(input_set.opt_tol);

	if ((input_set.max_eval != 0) && (input_set.max_eval < 3)) {
		log->warn("Searching for the optimum model parameters only for {} steps most probably will not be any useful!", input_set.max_eval);
	}

	if (input_set.diel_in.n_elem == 1) {
		log->trace("Isotropic dielectric constant inside of the slab.");
		input_set.diel_in = repelem(input_set.diel_in, 1, 3);
	}

	if (input_set.diel_out.n_elem == 1) {
		log->trace("Isotropic dielectric constant outside of the slab.");
		input_set.diel_out = repelem(input_set.diel_out, 1, 3);
	}

	if ((min(input_set.diel_in) <= 0) || (min(input_set.diel_out) <= 0)) {
		log->debug("Minimum of the dielectric tensor inside the slab: {}", min(input_set.diel_in));
		log->debug("Minimum of the dielectric tensor outside the slab: {}", min(input_set.diel_out));
		log->critical("The dielectric tensor is not defined properly! None of the tensor elements should be negative!");
		exit(1);
	}
	if (approx_equal(input_set.diel_in, input_set.diel_out, "absdiff", 0.02)) {
		log->debug("Model type: Bulk");
		if (input_set.optimize_interface) {
			log->trace("The position of the interfaces for the bulk models is irrelevant! \"interfaces\" will not be optimized!");
			input_set.optimize_interface = false;
		}
	}
	else {
		if (input_set.model_2D) {
			log->debug("Model type: 2D");
		}
		else {
			log->debug("Model type: Slab");
		}
	}
	if (input_set.optimize_charge_position || input_set.optimize_charge_sigma || input_set.optimize_charge_fraction || input_set.optimize_interface) {
		if ((input_set.opt_algo != "BOBYQA") && (input_set.opt_algo != "COBYLA") && (input_set.opt_algo != "SBPLX")) {
			log->debug("Optimization algorithm: {}", input_set.opt_algo);
			log->warn("Unsupported optimization algorithm is selected!");
			if (input_set.charge_position.n_rows < 3) {
				input_set.opt_algo = "BOBYQA";
			}
			else {
				input_set.opt_algo = "COBYLA";
			}
			log->warn("{} will be used instead!", input_set.opt_algo);
		}

		if ((input_set.optimize_charge_fraction) && (input_set.charge_fraction.n_elem == 1)) {
			log->debug("There is only 1 Gaussian charge in your slabcc model. The charge_fraction will not be optimized!");
			input_set.optimize_charge_fraction = false;
		}

		if ((input_set.opt_tol > 1) || (input_set.opt_tol < 0)) {
			log->debug("Requested optimization tolerance: {}", input_set.opt_tol);
			log->warn("The relative optimization tolerance is unacceptable! It must be in choosen in (0-1) range.");
			input_set.opt_tol = 0.01;
			log->warn("optimize_tolerance = {} will be used!", input_set.opt_tol);
		}
	}


	if (input_set.charge_position.n_cols != 3) {
		log->debug("Number of the parameters defined for the position of a charge: {}", input_set.charge_position.n_cols);
		log->critical("Incorrect definition of charge positions!");
		log->critical("Positions should be defined as: charge_position = 0.1 0.2 0.3; 0.1 0.2 0.4;");
		exit(1);
	}

	const uword charge_number = input_set.charge_position.n_rows;
	const uword sigma_rows = input_set.charge_sigma.n_rows;
	const uword sigma_cols = input_set.charge_sigma.n_cols;

	if (input_set.charge_sigma.n_elem == 1) {
		input_set.charge_sigma = input_set.charge_sigma(0) * ones(charge_number, 3);
		log->debug("Only one charge_sigma is defined!");

	}else if (sigma_rows == charge_number) {
		if (sigma_cols == 3) {
			if (input_set.trivariate) {
				log->debug("All the charge_sigma values are properly defined!");
			}
			else {
				log->warn("charge_sigma is not defined properly! charge_sigma={} will be used.", to_string(input_set.charge_sigma.col(0)));
			}
		}
		else if (sigma_cols == 1) {
			input_set.charge_sigma = repmat(input_set.charge_sigma, 1, 3);
			if (input_set.trivariate) {
				log->debug("Equal values will be assumed for the charge_sigma in all directions!");
			}
			else {
				log->debug("charge_sigma is properly defined!");
			}
		}
	}
	else if (sigma_cols == charge_number) {
		if (sigma_rows == 1) {
			input_set.charge_sigma = repmat(input_set.charge_sigma.t(), 1, 3);
			if (input_set.trivariate) {	
				log->debug("Equal values will be assumed for the charge_sigma in all directions!");
			}
			else {
				log->debug("charge_sigma is properly defined!");
			}
		}
	}
	else if ((sigma_cols == 3) && (sigma_rows == 1)) {
		if (input_set.trivariate) {
			input_set.charge_sigma = repmat(input_set.charge_sigma, charge_number, 1);
			log->debug("Equal charge_sigma will be assumed for all the Gaussian charges!");
		}
	}
	
	
	//if all the sensible checks and remedies have been failed!
	if(arma::size(input_set.charge_sigma) != arma::size(input_set.charge_position)){
		log->warn("Number of the defined Gaussian charges and the charge_sigma sets does not match!");
		input_set.charge_sigma = mat(arma::size(input_set.charge_position), fill::ones);
		if (input_set.trivariate) {
			log->warn("charge_sigma = {} will be used!", to_string(input_set.charge_sigma));
		}
		else {
			log->warn("charge_sigma = {} will be used!", to_string(input_set.charge_sigma.col(0)));
		}
	}

	log->debug("charge_sigma after the checks: {}", to_string(input_set.charge_sigma));

	//charge_fraction
	if (input_set.charge_fraction.n_elem != input_set.charge_position.n_rows) {
		log->debug("Number of charge position values: {}", input_set.charge_position.n_rows);
		log->debug("Number of charge fraction values: {}", input_set.charge_fraction.n_elem);
		input_set.charge_fraction = rowvec(input_set.charge_position.n_rows, fill::ones);
		log->warn("Number of the charge_fraction and charge_position sets does not match!");
		log->warn("Equal charge fractions will be assumed!");

	}


	if (!file_exists(input_set.CHGCAR_neutral) || !file_exists(input_set.CHGCAR_charged)
		|| !file_exists(input_set.LOCPOT_neutral) || !file_exists(input_set.LOCPOT_charged)) {
		log->debug("CHGCAR_neutral: '{}' found: {}", input_set.CHGCAR_neutral, to_string(file_exists(input_set.CHGCAR_neutral)));
		log->debug("CHGCAR_charged: '{}' found: {}", input_set.CHGCAR_charged, to_string(file_exists(input_set.CHGCAR_charged)));
		log->debug("LOCPOT_neutral: '{}' found: {}", input_set.LOCPOT_neutral, to_string(file_exists(input_set.LOCPOT_neutral)));
		log->debug("LOCPOT_charged: '{}' found: {}", input_set.LOCPOT_charged, to_string(file_exists(input_set.LOCPOT_charged)));
		log->critical("One or more of the input files could not be found!");
		exit(1);
	}
	
	if (input_set.extrapolate) {
		if (input_set.extrapol_steps_num < 3) {
			log->debug("Requested extrapolation steps: {}", input_set.extrapol_steps_num);
			log->warn("Minimum 3 steps are needed for the extrapolation algorithm to work reliably!");
			input_set.extrapol_steps_num = 3;
			log->warn("{} steps will be used instead!", input_set.extrapol_steps_num);
		}
	}
	else {
		if (input_set.model_2D) {
			// Bessel expansion limitations
			if (input_set.charge_fraction.n_elem > 1) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support multiple Gaussian charges. The extrapolation method \"extrapolate=yes\" will be used for this calculation.");
				input_set.extrapolate = true;
			}
			if (any(input_set.diel_out > 1)) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support the models embedded in any dielectric medium other than the vacuum. The extrapolation method \"extrapolate=yes\" will be used for this calculation.");
				input_set.extrapolate = true;
			}
			if (input_set.trivariate) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support the trivariate Gaussian model charges. The extrapolation method \"extrapolate=yes\" will be used for this calculation.");
				input_set.extrapolate = true;
			}

			vec2 inplane_diels;
			for (uword i = 0; i < 3; ++i) {
				if (i < input_set.normal_direction) {
					inplane_diels(i) = input_set.diel_in(i);
				}
				if (i > input_set.normal_direction) {
					inplane_diels(i - 1) = input_set.diel_in(i);
				}
			}

			if (abs(inplane_diels(0) - inplane_diels(1)) > 0.01) {
				log->debug("In-plane dielectric constants: {}", ::to_string(inplane_diels));
				log->error("In-plane dielectric constants are not equal. The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support this model. Will use the extrapolation method \"extrapolate=yes\" for this calculation.");
				input_set.extrapolate = true;
			}
		}
	}

	log->trace("Input parameters verified!");

}

void parse_input_params(const string& input_file, ofstream& output_fstream, const input_data& input_set) {
	auto log = spdlog::get("loggers");
	INIReader reader(input_file);
	if (reader.ParseError() < 0) {
		log->critical("Cannot load the input file: {}", input_file);
		exit(1);
	}

	verbosity_level = reader.GetInteger("verbosity", 1);

	input_set.CHGCAR_neutral = reader.GetStr("CHGCAR_neutral", "CHGCAR.N");
	input_set.LOCPOT_neutral = reader.GetStr("LOCPOT_neutral", "LOCPOT.N");
	input_set.LOCPOT_charged = reader.GetStr("LOCPOT_charged", "LOCPOT.C");
	input_set.CHGCAR_charged = reader.GetStr("CHGCAR_charged", "CHGCAR.C");
	input_set.charge_position = reader.GetMat("charge_position", {});
	input_set.charge_fraction = reader.GetVec("charge_fraction", rowvec(input_set.charge_position.n_rows, fill::ones) / input_set.charge_position.n_rows);
	input_set.trivariate = reader.GetBoolean("charge_trivariate", false);
	input_set.charge_sigma = reader.GetMat("charge_sigma", mat(input_set.charge_position.n_rows, 1, fill::ones));
	input_set.slabcenter = reader.GetVec("slab_center", { 0.5, 0.5, 0.5 });
	input_set.normal_direction = xyz2int(reader.GetStr("normal_direction", "z"));
	input_set.interfaces = reader.GetVec("interfaces", { 0.25, 0.75 });
	input_set.diel_in = reader.GetVec("diel_in", { 1 });
	input_set.diel_out = reader.GetVec("diel_out", { 1 });
	input_set.diel_erf_beta = reader.GetReal("diel_taper", 1);
	input_set.optimize_charge_position = reader.GetBoolean("optimize_charge_position", true);
	input_set.optimize_charge_sigma = reader.GetBoolean("optimize_charge_sigma", true);
	input_set.optimize_charge_fraction = reader.GetBoolean("optimize_charge_fraction", true);
	input_set.optimize_interface = reader.GetBoolean("optimize_interfaces", true);
	input_set.model_2D = reader.GetBoolean("2d_model", false);
	input_set.opt_algo = reader.GetStr("optimize_algorithm", (input_set.charge_position.n_rows == 1) ? "BOBYQA" : "COBYLA");
	input_set.opt_tol = reader.GetReal("optimize_tolerance", (input_set.opt_algo == "COBYLA") ? 5e-2 : 1e-2);
	input_set.max_eval = reader.GetInteger("optimize_maxsteps", 0);
	input_set.max_time = reader.GetInteger("optimize_maxtime", 0);
	input_set.opt_grid_x = reader.GetReal("optimize_grid_x", 0.8);
	input_set.extrapolate = reader.GetBoolean("extrapolate", input_set.model_2D ? false : true);
	input_set.extrapol_grid_x = reader.GetReal("extrapolate_grid_x", 1);
	input_set.extrapol_steps_num = reader.GetInteger("extrapolate_steps_number", 4);
	input_set.extrapol_steps_size = reader.GetReal("extrapolate_steps_size", 0.5);

	reader.dump_all(output_fstream);

	slabcc_cell.normal_direction = input_set.normal_direction;
}

double opt_charge_constraint(const vector<double> &x, vector<double> &grad, void *data)
{
	auto log = spdlog::get("loggers");
	rowvec2 interfaces;
	rowvec charge_fraction;
	mat charge_sigma;
	mat defcenter;
	opt_vars variables = { interfaces, charge_sigma, charge_fraction, defcenter };
	optimizer_unpacker(x, variables);

	const auto d = static_cast<const opt_data *>(data);
	const auto total_vasp_charge = d->total_vasp_charge;
	charge_fraction(charge_fraction.n_elem - 1) = 1 - accu(charge_fraction);
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
		exit(1);
	}

	log->trace("All files are loaded and cross-checked!");
}

void verify_interface_optimization(const rowvec2& initial_interfaces, const rowvec2& optimized_interfaces) {
	
	const double max_total_safe_movement = 4; //Ang
	const double normal_length = slabcc_cell.vec_lengths(slabcc_cell.normal_direction) / ang_to_bohr;
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