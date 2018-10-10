// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_core.hpp"

extern slabcc_cell_type slabcc_cell;

mat dielectric_profiles(const rowvec2 &interfaces, const rowvec3 &diel_in, const rowvec3 &diel_out, const double &diel_erf_beta) {
	const auto length = slabcc_cell.vec_lengths(slabcc_cell.normal_direction);
	const auto n_points = slabcc_cell.grid(slabcc_cell.normal_direction);
	rowvec2 interfaces_cartesian = interfaces * length;
	sort(interfaces_cartesian);
	const auto positions = linspace<rowvec>(0, length, n_points + 1);
	mat diels = zeros(n_points, 3);
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
		diels.row(k) = (diel_diff * diel_side * diel_edge + diel_sum) / 2;
	}

	return diels;
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
		cx_mat eps11_Gx0k2 = eps11 * square(Gx0(k));
		for (uword m = 0; m < Gy0.n_elem; ++m) {
			vector<span> spans = { span(k), span(m), span() };
			swap(spans.at(slabcc_cell.normal_direction), spans.at(2));
			cx_mat AG = Az + eps11_Gx0k2 + eps22 * square(Gy0(m));
			if ((k == 0) && (m == 0)) { AG(0, 0) = 1; }
			Vk(spans.at(0), spans.at(1), spans.at(2)) = solve(AG, vectorise(rhok(spans.at(0), spans.at(1), spans.at(2))));
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
	mat& diels = d->diels;
	cx_cube& rhoM = d->rhoM;
	cx_cube& V = d->V;
	cube& V_diff = d->V_diff;
	auto& initial_pot_MSE = d->initial_pot_MSE;

	diels = dielectric_profiles(interfaces, diel_in, diel_out, diel_erf_beta);

	rhoM.reset();
	rhoM = zeros<cx_cube>(as_size(slabcc_cell.grid));
	for (uword i = 0; i < Qd.n_elem; ++i) {
		rhoM += gaussian_charge(Qd(i), defcenter.row(i).t(), sigma(i));
	}

	V = poisson_solver_3D(rhoM, diels);
	V_diff = real(V) * Hartree_to_eV - defect_potential;

	const auto pot_MSE = accu(square(V_diff)) / V_diff.n_elem * 100;

	//if this is the 1st step, save the error for opt success checking
	if (initial_pot_MSE < 0) {
		initial_pot_MSE = pot_MSE;
	}

	auto log = spdlog::get("loggers");

	log->debug("-----------------------------------------");
	log->debug("> shifted_interfaces=" + ::to_string(x.at(0)) + " " + ::to_string(x.at(1)));

	for (uword i = 0; i < Qd.n_elem; ++i) {
		log->debug("> shifted_charge_position=" + ::to_string(defcenter.row(i)));
		log->debug("> charge_sigma=" + ::to_string(sigma(i)));
		if (Qd.n_elem > 1) {
			log->debug("> charge_fraction=" + to_string(abs(Qd(i) / accu(Qd))));
		}		
	}
	log->debug("Potential Root Mean Square Error: " + to_string(sqrt(pot_MSE)));

	return pot_MSE;
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
			log->warn("BOBYQA does not support the models with multiple charges! COBYLA will be used instead!");
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
	const int opt_parameters = opt_vars.Qd.n_elem * var_per_charge * optimize_charge + 2 * optimize_interfaces;
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
		input_set.diel_in = rowvec{ input_set.diel_in(0), input_set.diel_in(0), input_set.diel_in(0) };
	}

	if (input_set.diel_out.n_elem == 1) {
		input_set.diel_out = rowvec{ input_set.diel_out(0), input_set.diel_out(0), input_set.diel_out(0) };
	}

	if (input_set.opt_tol > 1) {
		input_set.opt_tol = 0.001;
		log->warn("The optimization tolerance is not defined properly");
		log->warn("Will use optimize_tolerance=" + to_string(input_set.opt_tol));
	}

	if (input_set.sigma.n_elem != input_set.charge_position.n_rows) {
		input_set.sigma = rowvec(input_set.charge_position.n_rows, fill::ones);
		log->warn("Number of the defined Gaussian charges and the sigma values does not match!");
		log->warn("Will use charge_sigma=" + to_string(input_set.sigma));
	}
	if (input_set.Qd.n_elem != input_set.sigma.n_elem) {
		input_set.Qd = rowvec(input_set.charge_position.n_rows, fill::ones);
		log->warn("Number of the charge_fraction and charges_sigma does not match!");
		log->warn("Equal charge fractions will be assumed!");

	}
	if (input_set.charge_position.n_cols != 3) {
		log->critical("Incorrect definition of charge positions!");
		log->critical("Positions should be defined as: charge_position = 0.1 0.2 0.3; 0.1 0.2 0.4;" );
		exit(1);
	}

	if ((max(input_set.diel_in < 0)) || (max(input_set.diel_out < 0))) {
		log->critical("The dielectric tensor is not defined properly!" );
		exit(1);
	}

	if (!file_exists(input_set.CHGCAR_NEU) || !file_exists(input_set.CHGCAR_CHG)
		|| !file_exists(input_set.LOCPOT_NEU) || !file_exists(input_set.LOCPOT_CHG)) {
		log->critical("One or more of the input files could not be found!");
		exit(1);
	}

	if ((input_set.opt_algo != "BOBYQA") && (input_set.opt_algo != "COBYLA")) {
		input_set.opt_algo = "COBYLA";
		log->warn("Unknown optimization algorithm is selected! {} will be used instead!", input_set.opt_algo);

	}

	if (input_set.extrapol_steps_num < 3) {
		log->warn("Extrapolation cannot be done with steps < 3" );
		input_set.extrapol_steps_num = 3;
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

	input_set.CHGCAR_NEU = reader.GetStr("CHGCAR_neutral", "CHGCAR.N");
	input_set.LOCPOT_NEU = reader.GetStr("LOCPOT_neutral", "LOCPOT.N");
	input_set.LOCPOT_CHG = reader.GetStr("LOCPOT_charged", "LOCPOT.C");
	input_set.CHGCAR_CHG = reader.GetStr("CHGCAR_charged", "CHGCAR.C");
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
	input_set.opt_algo = reader.GetStr("optimize_algorithm", "COBYLA");
	input_set.opt_tol = reader.GetReal("optimize_tolerance", 1e-3);
	input_set.max_eval = reader.GetInteger("optimize_maxsteps", 0);
	input_set.max_time = reader.GetInteger("optimize_maxtime", 0);
	input_set.opt_grid_x = reader.GetReal("optimize_grid_x", 0.8);
	input_set.extrapol_grid_x = reader.GetReal("extrapolate_grid_x", 1);
	input_set.extrapol_steps_num = reader.GetInteger("extrapolate_steps_number", 4);
	input_set.extrapol_steps_size = reader.GetReal("extrapolate_steps_size", 0.5);
	//input_set.extrapol_slab = reader.GetBoolean("extrapolate_slab", true);

	reader.dump_parsed(output_fstream);

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

tuple <rowvec, rowvec> extrapolate_3D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &Qd, const rowvec &sigma, const double &grid_multiplier) {
	auto log = spdlog::get("loggers");
	const uword normal_direction = slabcc_cell.normal_direction;
	rowvec Es = zeros<rowvec>(extrapol_steps_num - 1), sizes = Es;
	const mat33 cell_vectors0 = slabcc_cell.vectors;
	const urowvec grid0 = slabcc_cell.grid;
	const rowvec3 extrapolation_grid_size = grid_multiplier * conv_to<rowvec>::from(grid0);
	const urowvec3 extrapolation_grid = { (uword)extrapolation_grid_size(0),(uword)extrapolation_grid_size(1),(uword)extrapolation_grid_size(2) };
	for (auto n = 0; n < extrapol_steps_num - 1; ++n) {

		const auto extrapol_factor = extrapol_steps_size * (1.0 + n) + 1;

		UpdateCell(cell_vectors0 * extrapol_factor, extrapolation_grid);

		//extrapolated interfaces
		rowvec2 interfaces_ext = interfaces;
		const uvec interface_sorted_i = sort_index(interfaces);
		interfaces_ext(interface_sorted_i(1)) += abs(interfaces(0) - interfaces(1)) * (extrapol_factor - 1);
		interfaces_ext /= extrapol_factor;

		//charges moved to the same distance from their original nearest interface
		mat charge_position_shifted = charge_position / extrapol_factor;

		for (auto charge = 0; charge < charge_position.n_rows; ++charge) {
			const rowvec2 distance_to_interfaces = abs(charge_position(charge, normal_direction) - interfaces);
			if (distance_to_interfaces(0) < distance_to_interfaces(1)) {
				charge_position_shifted(charge, normal_direction) += interfaces_ext(0) - interfaces(0) / extrapol_factor;
			}
			else {
				charge_position_shifted(charge, normal_direction) += interfaces_ext(1) - interfaces(1) / extrapol_factor;
			}
		}

		const mat diels = dielectric_profiles(interfaces_ext, diel_in, diel_out, diel_erf_beta);

		cx_cube rhoM(as_size(slabcc_cell.grid), fill::zeros);
		for (uword i = 0; i < charge_position.n_rows; ++i) {
			rhoM += gaussian_charge(Qd(i), charge_position_shifted.row(i).t(), sigma(i));
		}
		const auto Q = accu(real(rhoM)) * slabcc_cell.voxel_vol;
		// (only works for the orthogonal cells!)
		rhoM -= Q / prod(slabcc_cell.vec_lengths);
		const auto V = poisson_solver_3D(rhoM, diels);
		const auto EperModel = 0.5 * accu(real(V % rhoM)) * slabcc_cell.voxel_vol * Hartree_to_eV;
		const rowvec2 interface_pos = interfaces_ext * slabcc_cell.vec_lengths(normal_direction);
		string extrapolation_info = to_string(extrapol_factor) + "\t" + ::to_string(EperModel) + "\t" + to_string(interface_pos);
		for (auto i = 0; i < charge_position_shifted.n_rows; ++i) {
			extrapolation_info +=  "\t" + to_string(charge_position_shifted(i, slabcc_cell.normal_direction) * slabcc_cell.vec_lengths(slabcc_cell.normal_direction));
		}
		log->debug(extrapolation_info);
		Es(n) = EperModel;
		sizes(n) = 1.0 / extrapol_factor;

	}

	return make_tuple(Es, sizes);
}

tuple <rowvec, rowvec> extrapolate_2D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &Qd, const rowvec &sigma, const double &grid_multiplier) {
	auto log = spdlog::get("loggers");
	const uword normal_direction = slabcc_cell.normal_direction;
	rowvec Es = zeros<rowvec>(extrapol_steps_num - 1), sizes = Es;
	const mat33 cell_vectors0 = slabcc_cell.vectors;
	const urowvec grid0 = slabcc_cell.grid;
	const rowvec3 extrapolation_grid_size = grid_multiplier * conv_to<rowvec>::from(grid0);
	const urowvec3 extrapolation_grid = { (uword)extrapolation_grid_size(0), (uword)extrapolation_grid_size(1), (uword)extrapolation_grid_size(2) };
	UpdateCell(cell_vectors0, extrapolation_grid);
	const mat diels0 = dielectric_profiles(interfaces, diel_in, diel_out, diel_erf_beta);

	for (auto n = 0; n < extrapol_steps_num - 1; ++n) {

		const auto extrapol_factor = extrapol_steps_size * (1.0 + n) + 1;
		UpdateCell(cell_vectors0 * extrapol_factor, extrapolation_grid);
		//extrapolated interfaces
		const rowvec2 interfaces_ext = interfaces / extrapol_factor;

		const mat charge_position_ext = charge_position / extrapol_factor;
		const mat diels = dielectric_profiles(interfaces_ext, diel_in, diel_out, diel_erf_beta);

		cx_cube rhoM(as_size(slabcc_cell.grid), fill::zeros);
		for (uword i = 0; i < charge_position.n_rows; ++i) {
			rhoM += gaussian_charge(Qd(i), charge_position_ext.row(i).t(), sigma(i));
		}
		const auto Q = accu(real(rhoM)) * slabcc_cell.voxel_vol;
		// (only works for the orthogonal cells!)
		rhoM -= Q / prod(slabcc_cell.vec_lengths);
		const auto V = poisson_solver_3D(rhoM, diels);
		const auto EperModel = 0.5 * accu(real(V % rhoM)) * slabcc_cell.voxel_vol * Hartree_to_eV;

		const rowvec2 interface_pos = interfaces_ext * slabcc_cell.vec_lengths(normal_direction);
		string extrapolation_info = to_string(extrapol_factor) + "\t" + ::to_string(EperModel) + "\t" + to_string(interface_pos);
		for (auto i = 0; i < charge_position_ext.n_rows; ++i) {
			extrapolation_info += "\t" + to_string(charge_position_ext(i, slabcc_cell.normal_direction) * slabcc_cell.vec_lengths(slabcc_cell.normal_direction));
		}
		log->debug(extrapolation_info);
		Es(n) = EperModel;
		sizes(n) = 1.0 / extrapol_factor;

	}

	return make_tuple(Es, sizes);
}


vector<double> nonlinear_fit(const double& opt_tol, nonlinear_fit_data& fit_data) {
	auto log = spdlog::get("loggers");
	double fit_MSE = 0;
	vector<double> fit_parameters = { 1, 1, 1, 1 };
	const auto opt_algorithm = nlopt::LN_COBYLA;

	nlopt::opt opt(opt_algorithm, 4);

	opt.set_min_objective(fit_eval, &fit_data);
	opt.set_xtol_rel(opt_tol);   //tolerance for error value

	try {
		const nlopt::result nlopt_final_result = opt.optimize(fit_parameters, fit_MSE);
	}
	catch (exception &e) {
		log->error("Nonlinear fitting failed: " + string(e.what()));
	}

	return fit_parameters;
}

double fit_eval(const vector<double> &c, vector<double> &grad, void *data)
{
	const auto d = static_cast<const nonlinear_fit_data *>(data);
	const rowvec scales = d->sizes;
	const rowvec energies = d->energies;
	const auto madelung_term = d->madelung_term;
	const rowvec model_energies = c.at(0) + c.at(1) * scales + c.at(2) * square(scales) + (c.at(1) - madelung_term) / c.at(3) * exp(-c.at(3) * scales);
	const auto fit_MSE = accu(square(energies - model_energies));
	return fit_MSE;
}

void check_cells(const supercell& Neutral_supercell, const supercell& Charged_supercell, const input_data& input_set) {

	auto log = spdlog::get("loggers");

	// equal cell size
	if (!approx_equal(Neutral_supercell.cell_vectors * Neutral_supercell.scaling,
		Charged_supercell.cell_vectors * Charged_supercell.scaling, "reldiff", 0.001)) {
		log->critical("Size vectors of the input files does not match!");
		exit(1);
	}

	// cell needs rotation or is not orthogonal
	const vec cellvec_nonzeros = nonzeros(Neutral_supercell.cell_vectors);
	if (cellvec_nonzeros.n_elem != 3) {
		log->critical("unsupported supercell shape!");
		exit(1);
	}

	// equal grid
	const SizeCube input_grid = arma::size(Neutral_supercell.potential);
	if ((input_grid != arma::size(Charged_supercell.potential)) ||
		(input_grid != arma::size(Charged_supercell.charge)) ||
		(input_grid != arma::size(Neutral_supercell.charge))) {

		log->critical("Data grids of the CHGCAR or LOCPOT files does not match!");
		exit(1);
	}

	log->trace("All files are loaded and cross-checked!");
}
