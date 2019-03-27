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


cx_cube gaussian_charge(const double& Q, const vec3& rel_pos, const double& sigma) {

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
	const cx_cube charge_dist = cx_cube(Q / pow((sigma * sqrt(2 * PI)), 3) * exp(-r2 / (2 * square(sigma))), zeros(as_size(slabcc_cell.grid)));
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

	rowvec2 interfaces;
	rowvec sigma, Qd;
	mat defcenter;
	opt_vars variables = { interfaces, sigma, Qd, defcenter };
	optimizer_unpacker(x, variables);

	//input data
	const auto d = static_cast<opt_data *>(slabcc_data);
	const rowvec3 &diel_in = d->diel_in;
	const rowvec3 &diel_out = d->diel_out;
	const auto &diel_erf_beta = d->diel_erf_beta;
	const cube &defect_potential = d->defect_potential;
	const auto &Q0 = d->Q0;

	//rest of the charge goes to the last Gaussian
	Qd(Qd.n_elem - 1) = Q0 - accu(Qd);

	//output data
	mat& dielectric_profiles = d->dielectric_profiles;
	cx_cube& rhoM = d->rhoM;
	cx_cube& V = d->V;
	cube& V_diff = d->V_diff;
	auto& initial_potential_MSE = d->initial_potential_MSE;

	dielectric_profiles = dielectric_profiles_gen(interfaces, diel_in, diel_out, diel_erf_beta);

	rhoM.reset();
	rhoM = zeros<cx_cube>(as_size(slabcc_cell.grid));
	for (uword i = 0; i < Qd.n_elem; ++i) {
		rhoM += gaussian_charge(Qd(i), defcenter.row(i).t(), sigma(i));
	}

	V = poisson_solver_3D(rhoM, dielectric_profiles);
	V_diff = real(V) * Hartree_to_eV - defect_potential;

	const auto potential_MSE = accu(square(V_diff)) / V_diff.n_elem;

	if (initial_potential_MSE < 0) {
		initial_potential_MSE = potential_MSE;
	}

	auto log = spdlog::get("loggers");

	log->debug("-----------------------------------------");
	log->debug("> shifted_interfaces=" + ::to_string(x.at(0)) + " " + ::to_string(x.at(1)));

	for (uword i = 0; i < Qd.n_elem; ++i) {
		log->debug("> shifted_charge_position=" + ::to_string(defcenter.row(i)));
		log->debug("> charge_sigma=" + ::to_string(sigma(i)));
		if (Qd.n_elem > 1) {
			log->debug("> charge_fraction={}", abs(Qd(i) / accu(Qd)));
		}		
	}
	log->debug("Potential Root Mean Square Error: {}", sqrt(potential_MSE));

	return potential_MSE;
}

double do_optimize(const string& opt_algo, const double& opt_tol, const int &max_eval, const int &max_time, opt_data& opt_data, opt_vars& opt_vars, const bool &optimize_charge, const bool &optimize_interfaces) {
	auto log = spdlog::get("loggers");
	double pot_MSE_min = 0;
	auto opt_algorithm = nlopt::LN_COBYLA;
	if (opt_algo == "BOBYQA") {
		if (opt_vars.Qd.n_elem == 1) {
			opt_algorithm = nlopt::LN_BOBYQA;
		}
		else {
			//BOBYQA in NLOPT 2.4.2 does not support the constraints!
			log->warn("BOBYQA optimization algorithm does not support the models with multiple charges! COBYLA will be used instead!");
		}
	}

	vector<double> opt_param, low_b, upp_b;
	tie(opt_param, low_b, upp_b) = optimizer_packer(opt_vars, optimize_charge, optimize_interfaces);
	nlopt::opt opt(opt_algorithm, opt_param.size());

	opt.set_lower_bounds(low_b);
	opt.set_upper_bounds(upp_b);

	opt.set_min_objective(potential_eval, &opt_data);
	opt.set_xtol_rel(square(opt_tol));   //tolerance is defined for RMSE in the input, optimizer is checking MSE
	if (max_eval > 0) {
		opt.set_maxeval(max_eval);
	}
	if (max_time > 0) {
		opt.set_maxtime(60.0 * max_time);
	}
	if (opt_vars.Qd.n_elem > 1) {
		//add constraint to keep the total charge constant
		opt.add_inequality_constraint(opt_charge_constraint, &opt_data, 1e-8);
	}

	const int var_per_charge = (opt_vars.Qd.n_elem == 1) ? 4 : 5;
	const uword opt_parameters = opt_vars.Qd.n_elem * var_per_charge * optimize_charge + 2 * optimize_interfaces;
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
	catch (exception &e) {
		log->error("Parameters optimization failed: " + string(e.what()));
	}

	optimizer_unpacker(opt_param, opt_vars);
	log->trace("Optimization ended.");

	return pot_MSE_min;
}

tuple<vector<double>, vector<double>, vector<double>> optimizer_packer(const opt_vars& opt_vars, const bool optimize_charge, const bool optimize_interface) {

	vector<double> opt_param = { opt_vars.interfaces(0), opt_vars.interfaces(1) };
	vector<double> low_b = { 0, 0 };		//lower bounds
	vector<double> upp_b = { 1, 1 };		//upper bounds
	if (!optimize_interface) {
		low_b = opt_param;
		upp_b = opt_param;
	}

	for (uword i = 0; i < opt_vars.charge_position.n_rows; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_param.push_back(opt_vars.charge_position(i, j));
			if (optimize_charge) {
				low_b.push_back(0);
				upp_b.push_back(1);
			}
			else {
				low_b.push_back(opt_vars.charge_position(i, j));
				upp_b.push_back(opt_vars.charge_position(i, j));
			}
		}
		opt_param.push_back(opt_vars.sigma(i));
		opt_param.push_back(opt_vars.Qd(i));
		const auto min_charge = min(0.0, accu(opt_vars.Qd));
		const auto max_charge = max(0.0, accu(opt_vars.Qd));
		if (optimize_charge) {
			low_b.insert(low_b.end(), { 0.1, min_charge });	//sigma, q
			upp_b.insert(upp_b.end(), { 7, max_charge });
		}
		else {
			low_b.insert(low_b.end(), { opt_vars.sigma(i), opt_vars.Qd(i) });
			upp_b.insert(upp_b.end(), { opt_vars.sigma(i), opt_vars.Qd(i) });
		}
	}

	//remove the last Qd from the parameters 
	//this gives a better chance to optimizer and also removes the Qd if we have only one charge
	opt_param.pop_back();
	upp_b.pop_back();
	low_b.pop_back();

	return make_tuple(opt_param, low_b, upp_b);
}

void optimizer_unpacker(const vector<double> &optimizer_vars_vec, opt_vars &opt_vars) {
	// the coefficients in this part depend on the set of our variables
	// they are ordered as: interfaces, :|[x, y, z, sigma, q]|:
	// NOTE: the last "q" must be calculated from the total charge (which in not available here!)
	// NOTE2: this function is also called by opt_charge_constraint() with a reference to empty variable.
	const uword variable_per_charge = 5;
	const uword n_charges = optimizer_vars_vec.size() / variable_per_charge;
	opt_vars.sigma = zeros<rowvec>(n_charges);
	opt_vars.Qd = zeros<rowvec>(n_charges);
	opt_vars.charge_position = zeros(n_charges, 3);
	opt_vars.interfaces = { optimizer_vars_vec.at(0), optimizer_vars_vec.at(1) };

	for (uword i = 0; i < n_charges; ++i) {
		for (uword j = 0; j < 3; ++j) {
			opt_vars.charge_position(i, j) = optimizer_vars_vec.at(2 + i * variable_per_charge + j);
		}
		opt_vars.sigma(i) = optimizer_vars_vec.at(5 + variable_per_charge * i);
		if (i != n_charges - 1) {
			opt_vars.Qd(i) = optimizer_vars_vec.at(6 + variable_per_charge * i);
		}
	}
}


void check_inputs(input_data input_set) {
	auto log = spdlog::get("loggers");
	input_set.sigma = abs(input_set.sigma);
	input_set.max_eval = abs(input_set.max_eval);
	input_set.max_time = abs(input_set.max_time);
	input_set.interfaces = fmod_p(input_set.interfaces, 1);
	input_set.extrapol_grid_x = abs(input_set.extrapol_grid_x);
	input_set.opt_grid_x = abs(input_set.opt_grid_x);
	input_set.opt_tol = abs(input_set.opt_tol);

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
		log->critical("The dielectric tensor is not defined properly!");
		exit(1);
	}
	if (input_set.optimize_charge || input_set.optimize_interface) {
		if ((input_set.opt_algo != "BOBYQA") && (input_set.opt_algo != "COBYLA")) {

			log->debug("Optimization algorithm: {}", input_set.opt_algo);
			log->warn("Unknown optimization algorithm is selected!");
			input_set.opt_algo = "COBYLA";
			log->warn("{} will be used instead!", input_set.opt_algo);

		}

		if ((input_set.opt_tol > 1) || (input_set.opt_tol < 0)) {
			log->debug("Requested optimization tolerance: {}", input_set.opt_tol);
			log->warn("The optimization tolerance is unacceptable!");
			input_set.opt_tol = 0.001;
			log->warn("Will use optimize_tolerance = {}", input_set.opt_tol);
		}
	}


	if (input_set.sigma.n_elem != input_set.charge_position.n_rows) {
		input_set.sigma = rowvec(input_set.charge_position.n_rows, fill::ones);
		log->debug("Number of sigma values: {}", input_set.sigma.n_elem);
		log->debug("Number of Gaussian charges: {}", input_set.charge_position.n_rows);
		log->warn("Number of the defined Gaussian charges and the sigma values does not match!");
		log->warn("Will use charge_sigma = {}", to_string(input_set.sigma));
	}
	if (input_set.Qd.n_elem != input_set.sigma.n_elem) {
		log->debug("Number of sigma values: {}", input_set.sigma.n_elem);
		log->debug("Number of charge fraction values: {}", input_set.Qd.n_elem);
		input_set.Qd = rowvec(input_set.charge_position.n_rows, fill::ones);
		log->warn("Number of the charge_fraction and charges_sigma does not match!");
		log->warn("Equal charge fractions will be assumed!");

	}
	if (input_set.charge_position.n_cols != 3) {
		log->debug("Number of the parameters for position of each charge: {}", input_set.charge_position.n_cols);
		log->critical("Incorrect definition of charge positions!");
		log->critical("Positions should be defined as: charge_position = 0.1 0.2 0.3; 0.1 0.2 0.4;" );
		exit(1);
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
			log->warn("Minimum 3 steps are needed for extrapolation!");
			input_set.extrapol_steps_num = 3;
			log->warn("{} steps will be used instead!", input_set.extrapol_steps_num);
		}
	}
	else {
		if (input_set.model_2D) {
			// Bessel expansion limitations
			if (input_set.Qd.n_elem > 1) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support multiple Gaussian charges. Will use the extrapolation method \"extrapolate=yes\" for this calculation.");
				input_set.extrapolate = true;
			}
			if (any(input_set.diel_out > 1)) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support the models embedded in any dielectric medium other than the vacuum. Will use the extrapolation method \"extrapolate=yes\" for this calculation.");
				input_set.extrapolate = true;
			}

			vec2 inplace_diels;
			for (uword i = 0; i < 3; ++i) {
				if (i < input_set.normal_direction) {
					inplace_diels(i) = input_set.diel_in(i);
				}
				if (i > input_set.normal_direction) {
					inplace_diels(i - 1) = input_set.diel_in(i);
				}
			}

			if (abs(inplace_diels(0) - inplace_diels(1)) > 0.01) {
				log->debug("In-plane dielectric constants: {}", ::to_string(inplace_diels));
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
		log->critical("Cannot load '{}'", input_file);
		exit(1);
	}

	verbos = reader.GetInteger("verbosity", 1);

	input_set.CHGCAR_neutral = reader.GetStr("CHGCAR_neutral", "CHGCAR.N");
	input_set.LOCPOT_neutral = reader.GetStr("LOCPOT_neutral", "LOCPOT.N");
	input_set.LOCPOT_charged = reader.GetStr("LOCPOT_charged", "LOCPOT.C");
	input_set.CHGCAR_charged = reader.GetStr("CHGCAR_charged", "CHGCAR.C");
	input_set.charge_position = reader.GetMat("charge_position", {});
	input_set.Qd = reader.GetVec("charge_fraction", rowvec(input_set.charge_position.n_rows, fill::ones) / input_set.charge_position.n_rows);
	input_set.sigma = reader.GetVec("charge_sigma", rowvec(input_set.charge_position.n_rows, fill::ones));
	input_set.slabcenter = reader.GetVec("slab_center", { 0.5, 0.5, 0.5 });
	input_set.normal_direction = xyz2int(reader.GetStr("normal_direction", "z"));
	input_set.interfaces = reader.GetVec("interfaces", { 0.25, 0.75 });
	input_set.diel_in = reader.GetVec("diel_in", { 1 });
	input_set.diel_out = reader.GetVec("diel_out", { 1 });
	input_set.diel_erf_beta = reader.GetReal("diel_taper", 1);
	input_set.optimize_charge = reader.GetBoolean("optimize_charge", true);
	input_set.optimize_interface = reader.GetBoolean("optimize_interfaces", true);
	input_set.model_2D = reader.GetBoolean("2d_model", false);
	input_set.opt_algo = reader.GetStr("optimize_algorithm", "COBYLA");
	input_set.opt_tol = reader.GetReal("optimize_tolerance", 1e-3);
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
	rowvec sigma, Qd;
	mat defcenter;
	opt_vars variables = { interfaces, sigma, Qd, defcenter };
	optimizer_unpacker(x, variables);

	const auto d = static_cast<const opt_data *>(data);
	const auto Q0 = d->Q0;
	Qd(Qd.n_elem - 1) = Q0 - accu(Qd);
	const auto constraint = -Qd(Qd.n_elem - 1)  * Q0;
	log->debug("Charge in each Gaussian:" + ::to_string(Qd));
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

void verify_optimization(const rowvec2& initial_interfaces, const rowvec2& optimized_interfaces) {
	const double max_total_safe_movement = 4;
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