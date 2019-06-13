// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "stdafx.h"
#include "isolated.hpp"
#include "vasp.hpp"
#include "slabcc_model.hpp"
using namespace std;

int verbosity_level = 0;
const double ang_to_bohr = 1e-10 / datum::a_0;  //1.88972612546
const double Hartree_to_eV = datum::R_inf * datum::h * datum::c_0 / datum::eV * 2; //27.2113860193

int main(int argc, char *argv[]){

	slabcc_model model;
	string input_file = "slabcc.in";
	string output_file = "slabcc.out";
	string log_file = "slabcc.log";
	bool output_diffs_only = false;
	cli_params parameters_list = { input_file, output_file, log_file, output_diffs_only };
	parameters_list.parse(argc, argv);
	initialize_loggers(log_file, output_file);
	auto log = spdlog::get("loggers");
	auto output_log = spdlog::get("output");

	// default values should be defined in the parser function

	string CHGCAR_neutral = "";
	string LOCPOT_neutral = "";
	string LOCPOT_charged = "";
	string CHGCAR_charged = "";
	string opt_algo = "";			//optimization algorithm
	mat charge_position;			//center of each Gaussian model charge
	rowvec charge_fraction;			//charge fraction in each Gaussian
	mat charge_sigma;				//width of each Gaussian model charges
	mat charge_rotations;			//rotation angles along each axis for the trivariate Gaussians
	bool charge_trivariate = false; //use trivariate Gaussians
	rowvec diel_in;					//diagonal elements of slab dielectric tensor
	rowvec diel_out;				//diagonal elements of enviroment dielectric tensor
	rowvec3 slabcenter;
	uword normal_direction = 0;		//index of the normal direction (0/1/2)
	rowvec2 interfaces;				//interfaces in relative coordinates, ordered as the user input 
	double diel_erf_beta = 0;		//beta value of the erf for dielectric profile generation
	double opt_tol = 0;				//relative optimization tolerance
	double extrapol_grid_x = 0;		//extrapolation grid size multiplier
	double opt_grid_x = 0;			//optimization grid size multiplier
	int max_eval = 0;				//maximum number of steps for the optimization function evaluation
	int max_time = 0;				//maximum time for the optimization in minutes
	int extrapol_steps_num = 0;		//number of extrapolation steps for E_isolated calculation
	double extrapol_steps_size = 0; //size of each extrapolation step with respect to the initial supercell size
	bool optimize_charge_position = false;	//optimize the charge_position 
	bool optimize_charge_sigma = false;		//optimize the charge_sigma
	bool optimize_charge_rotation = false;	//optimize the charge_rotation 
	bool optimize_charge_fraction = false;	//optimize the charge_fraction
	bool optimize_interfaces = false;		//optimize the position of interfaces
	bool extrapolate = false;	//use the extrapolation for E-isolated calculations
	bool model_2D = false;		//the model is 2D
	
	// parameters read from the input file
	const input_data inputfile_variables = {
		CHGCAR_neutral, LOCPOT_charged, LOCPOT_neutral, CHGCAR_charged,
		opt_algo, charge_position, charge_fraction, charge_sigma, charge_rotations, slabcenter, diel_in, diel_out,
		normal_direction, interfaces, diel_erf_beta,
		opt_tol, optimize_charge_position, optimize_charge_sigma, optimize_charge_rotation, optimize_charge_fraction, optimize_interfaces, extrapolate, model_2D, charge_trivariate, opt_grid_x,
		extrapol_grid_x, max_eval, max_time, extrapol_steps_num, extrapol_steps_size };

	inputfile_variables.parse(input_file);
	if (!output_diffs_only) {
		inputfile_variables.verify();
	}
	model.set_input_variables(inputfile_variables);

	log->debug("SLABCC: version {}.{}.{}", SLABCC_VERSION_MAJOR, SLABCC_VERSION_MINOR, SLABCC_VERSION_PATCH);
	log->debug("Armadillo library: version {}.{}.{}", ARMA_VERSION_MAJOR, ARMA_VERSION_MINOR, ARMA_VERSION_PATCH);
	log->debug("NLOPT library: version {}.{}.{}", nlopt::version_major(), nlopt::version_minor(), nlopt::version_bugfix());
	log->debug("SPDLOG library: version {}.{}.{}", SPDLOG_VER_MAJOR, SPDLOG_VER_MINOR, SPDLOG_VER_PATCH);
	log->debug("SLABCC input file: {}", input_file);
	log->debug("SLABCC output file: {}", output_file);
	log->debug("SLABCC log file: {}", log_file);

	vector<pair<string, string>> calculation_results;


	//promises for async read of CHGCAR and POTCAR files
	vector<future<cube>> future_cells;

	//promises for async file writing (can be replaced by a deque if the number of files increases)
	vector<future<void>> future_files;

	future_cells.push_back(async(launch::async, read_CHGPOT, CHGCAR_neutral));
	future_cells.push_back(async(launch::async, read_CHGPOT, CHGCAR_charged));
	future_cells.push_back(async(launch::async, read_CHGPOT, LOCPOT_neutral));
	future_cells.push_back(async(launch::async, read_CHGPOT, LOCPOT_charged));

	supercell Neutral_supercell = read_POSCAR(CHGCAR_neutral);
	supercell Charged_supercell = read_POSCAR(CHGCAR_charged);

	Neutral_supercell.charge = future_cells.at(0).get();
	Charged_supercell.charge = future_cells.at(1).get();
	Neutral_supercell.potential = future_cells.at(2).get();
	Charged_supercell.potential = future_cells.at(3).get();

	check_slabcc_compatiblity(Neutral_supercell, Charged_supercell);

	//cell vectors of the CHGCAR and LOCPOT files (bohr)
	const mat33 input_cell_vectors = abs(Neutral_supercell.cell_vectors) * Neutral_supercell.scaling * ang_to_bohr;
	const urowvec3 input_grid_size = SizeVec(Neutral_supercell.charge);
	model.init_supercell(input_cell_vectors, input_grid_size);

	//(only works in the orthogonal case!)
	model.cell_volume = prod(model.cell_vectors_lengths);

	const rowvec3 relative_shift = 0.5 - slabcenter;

	model.rounded_relative_shift = round(model.cell_grid % relative_shift) / model.cell_grid;

	model.interfaces = fmod(model.interfaces + model.rounded_relative_shift(normal_direction), 1);

	if (!output_diffs_only) {
		shift_structure(Neutral_supercell, model.rounded_relative_shift);
		shift_structure(Charged_supercell, model.rounded_relative_shift);
		model.charge_position += repmat(model.rounded_relative_shift, model.charge_position.n_rows, 1);
		model.charge_position = fmod_p(model.charge_position, 1);
	}
	log->debug("Slab normal direction index (0-2): {}", model.normal_direction);
	log->trace("Shift to center done!");

	supercell Defect_supercell = Neutral_supercell;
	Defect_supercell.potential = Charged_supercell.potential - Neutral_supercell.potential;
	Defect_supercell.charge = Charged_supercell.charge - Neutral_supercell.charge;

	if (is_active(verbosity::write_defect_file) || output_diffs_only) {
		future_files.push_back(async(launch::async, write_CHGPOT, "LOCPOT", "slabcc_D.LOCPOT", Defect_supercell));
		future_files.push_back(async(launch::async, write_CHGPOT, "CHGCAR", "slabcc_D.CHGCAR", Defect_supercell));
	}

	if (output_diffs_only) {
		log->debug("Only the extra charge and the potential difference calculation have been requested!");
		write_planar_avg(Defect_supercell.potential, Defect_supercell.charge * model.voxel_vol, "D");
		for (auto &promise : future_files) { promise.get(); }
		finalize_loggers();
		exit(0);
	}

	//normalize the charges and potentials
	Neutral_supercell.charge *= -1.0 / model.cell_volume;
	Charged_supercell.charge *= -1.0 / model.cell_volume;
	Defect_supercell.charge *= -1.0 / model.cell_volume;
	Defect_supercell.potential *= -1.0;
	model.V_target0 = Defect_supercell.potential;

	// total extra charge of the VASP calculation
	model.defect_charge = accu(Defect_supercell.charge) * model.voxel_vol;

	if (abs(model.defect_charge) < 0.001) {
		log->debug("Total extra charge: {}", model.defect_charge);
		log->warn("Total extra charge seems to be very small. Please make sure the path to the input CHGCAR files are set properly!");
	}
	const opt_switches optimizer_activation_switches{ optimize_charge_position, optimize_charge_sigma, optimize_charge_rotation, optimize_charge_fraction, optimize_interfaces };
	const bool optimize_any = optimize_charge_position || optimize_charge_sigma || optimize_charge_rotation || optimize_charge_fraction || optimize_interfaces;

	if (optimize_any) {
		const rowvec2 shifted_interfaces0 = model.interfaces;
		const mat charge_position0 = model.charge_position;
		const rowvec3 optimization_grid_size = opt_grid_x * conv_to<rowvec>::from(model.cell_grid);
		const urowvec3 optimization_grid = { (uword)optimization_grid_size(0), (uword)optimization_grid_size(1), (uword)optimization_grid_size(2) };
		model.change_grid(optimization_grid);
		model.update_V_target();
		optimize(opt_algo, opt_tol, max_eval, max_time, model, optimizer_activation_switches);
		model.change_grid(input_grid_size);

		//write the unshifted optimized values to the file
		output_log->info("\n[Optimized_model_parameters]");
		if (optimize_interfaces) {
			const rowvec2 optimized_interfaces = fmod_p(model.interfaces - model.rounded_relative_shift(model.normal_direction), 1);
			output_log->info("interfaces_optimized = {}", to_string(optimized_interfaces));
		}

		if (optimize_charge_fraction) {
			output_log->info("charge_fraction_optimized = {}", to_string(model.charge_fraction));
		}
		if (optimize_charge_sigma) {
			const mat opt_charge_sigma = model.trivariate_charge ? model.charge_sigma : model.charge_sigma.col(0);
			output_log->info("charge_sigma_optimized = {}", to_string(opt_charge_sigma));
			model.verify_charge_optimization();
		
		}
		if (optimize_charge_rotation) {
			const mat rotations = model.charge_rotations * 180.0 / PI;
			output_log->info("charge_rotation_optimized = {}", to_string(rotations));
		}
		if (optimize_charge_position) {
			const mat optimized_charge_position = fmod_p(model.charge_position - repmat(model.rounded_relative_shift, model.charge_position.n_rows, 1), 1);
			output_log->info("charge_position_optimized = {}", to_string(optimized_charge_position));

			// TODO: need a better algorithm to handle the swaps and PBCs similar to verify_interface_optimization()
			const mat charge_position_change = abs(charge_position0 - model.charge_position);
			if (charge_position_change.max() > 0.1) {
				log->warn("The optimized position for the extra charge is significantly different from the initial value. "
							"Please make sure that the final position of the extra charge have been estimated correctly!");
				log->debug("Charge position changes: ", to_string(charge_position_change));
			}
		}

		if (optimize_interfaces) {
			model.verify_interface_optimization(shifted_interfaces0);
		}

		if (model.initial_potential_RMSE * (opt_tol + 1) < model.potential_RMSE) {
			// Don't panic! either NLOPT seems to be malfunctioning
			// or we are not correctly logging/checking the result
			log->critical("Optimization failed!");
			log->critical("Potential error of the initial parameters seems to be smaller than the optimized parameters! "
							"You may want to change the initial guess for charge_position, change the optimization algorithm, or turn off the optimization.");
			log->debug("Initial model potential RMSE: {}", model.initial_potential_RMSE);
			log->debug("Optimized model potential RMSE: {}", model.potential_RMSE);
			finalize_loggers();
			exit(1);
		}
	}
	auto local_param = model.data_packer();
	vector<double> gradients = {};
	model.potential_RMSE = potential_error(get<0>(local_param), gradients, &model);
	const bool bulk_model = approx_equal(diel_in, diel_out, "absdiff", 0.02);
	const bool isotropic_screening = (abs(diel_in(0) - diel_in(1)) < 0.02) && (abs(diel_in(0) - diel_in(2)) < 0.02);

	if (model.potential_RMSE > 0.1) {
		if (bulk_model && isotropic_screening) {
			log->debug("RMSE of the model charge potential is large but for the bulk models with an isotropic screening (dielectric tensor) "
			"this shouldn't make much difference in the total correction energy!");
		}
		else {
			log->warn("RMSE of the model charge potential is large. The calculated correction energies may not be accurate!");
		}
	}

	const auto V_error_x = conv_to<rowvec>::from(planar_average(0, model.V_diff));
	const auto V_error_y = conv_to<rowvec>::from(planar_average(1, model.V_diff));
	const auto V_error_z = conv_to<rowvec>::from(planar_average(2, model.V_diff));

	//potential error in each direction
	rowvec3 V_error_planars = { accu(square(V_error_x)), accu(square(V_error_y)), accu(square(V_error_z)) };
	V_error_planars = sqrt(V_error_planars/ model.V_diff.n_elem);
	log->debug("Directional RMSE: " + to_string(V_error_planars));
	if (max(V_error_planars) / min(V_error_planars) > 10) {
		log->warn("The potential error is highly anisotropic.");
		log->warn("If the potential error is large, this usually means that either the extra charge is not properly described by the model Gaussian charge "
			"or the chosen dielectric tensor is not a good representation of the actual tensor! "
			"This can be fixed by properly optimizing the model parameters, using multiple Gaussian charges, using trivaritate Gaussians, or using the dielectric tensor calculated for the same VASP model/method.");
	}

	log->debug("Potential error anisotropy: {}", max(V_error_planars) / min(V_error_planars));
	log->debug("Cell dimensions (bohr): " + to_string(model.cell_vectors_lengths));
	log->debug("Volume (bohr^3): {}", model.cell_volume);


	if (is_active(verbosity::write_defect_file)) {
		supercell Model_supercell = Neutral_supercell;
		//charge is normalized to the VASP CHGCAR convention (rho * Vol)
		//Also, positive value for the electron charge! (the probability of finding an electron)
		Model_supercell.charge = -real(model.CHG) * model.voxel_vol * model.CHG.n_elem;
		Model_supercell.potential = -real(model.V) * Hartree_to_eV;
		future_files.push_back(async(launch::async, write_CHGPOT, "CHGCAR", "slabcc_M.CHGCAR", Model_supercell));
		future_files.push_back(async(launch::async, write_CHGPOT, "LOCPOT", "slabcc_M.LOCPOT", Model_supercell));
	}

	if (is_active(verbosity::write_dielectric_file)) {
		write_mat2file(model.dielectric_profiles, "slabcc_DIEL.dat");
	}
	if (is_active(verbosity::write_planarAvg_file)) {
		write_planar_avg(Neutral_supercell.potential, Neutral_supercell.charge * model.voxel_vol, "N");
		write_planar_avg(Charged_supercell.potential, Charged_supercell.charge * model.voxel_vol, "C");
		write_planar_avg(model.V_target0, Defect_supercell.charge * model.voxel_vol, "D");
		write_planar_avg(real(model.V) * Hartree_to_eV, real(model.CHG) * model.voxel_vol, "M");
	}
	else if (is_active(verbosity::write_normal_planarAvg)) {
		write_planar_avg(model.V_target0, Defect_supercell.charge * model.voxel_vol, "D", model.normal_direction);
		write_planar_avg(real(model.V) * Hartree_to_eV, real(model.CHG) * model.voxel_vol, "M", model.normal_direction);
	}

	//add jellium to the charge (Because the V is normalized, it is not needed in solving the Poisson eq. but it is needed in the energy calculations)
	model.CHG -= model.total_charge / model.cell_volume;
	const uword farthest_element_index = model.total_charge < 0 ? real(model.V).index_max() : real(model.V).index_min();

	const auto dV = model.V_diff(farthest_element_index);
	log->info("Potential alignment (dV=): {}", ::to_string(dV));
	calculation_results.emplace_back("dV", ::to_string(dV));

	if (abs(dV) > 0.05) {
		if (bulk_model && isotropic_screening) {
			log->debug("The potential alignment term (dV) is relatively large. But in the isotropic bulk models "
			"this should not make much difference in the total energy correction value!");
		}
		else {
			log->warn("The potential alignment term (dV) is relatively large. The constructed model may not be accurate!");
		}
	}

	log->debug("Calculation grid point for the potential alignment term: {}", to_string(ind2sub(as_size(model.cell_grid), farthest_element_index)));

	const double EperModel0 = 0.5 * accu(real(model.V) % real(model.CHG)) * model.voxel_vol * Hartree_to_eV;
	log->info("E_periodic of the model charge: {}", ::to_string(EperModel0));
	calculation_results.emplace_back("E_periodic of the model charge", ::to_string(EperModel0));

	log->debug("Difference of the charge in the input files: {}", ::to_string(model.defect_charge));
	log->debug("Total charge of the model: {}", ::to_string(model.total_charge));

	double E_isolated = 0;
	double E_correction = 0;
	if (extrapolate) {

		const rowvec3 extrapolation_grid_size = extrapol_grid_x * conv_to<rowvec>::from(model.cell_grid);
		const urowvec3 extrapolation_grid = { (uword)extrapolation_grid_size(0), (uword)extrapolation_grid_size(1), (uword)extrapolation_grid_size(2) };
		model.change_grid(extrapolation_grid);
		model.update_V_target();
		model.adjust_extrapolation_params(extrapol_steps_num, extrapol_steps_size, extrapol_grid_x);
		
		log->debug("--------------------------------------------------------");
		log->debug("Scaling\tE_periodic\t\tmodel charge\t\tinterfaces\t\tcharge position");
		const rowvec2 interface_pos = model.interfaces * model.cell_vectors_lengths(model.normal_direction);
		string extrapolation_info = to_string(1.0) + "\t" + ::to_string(EperModel0) + "\t" + ::to_string(model.total_charge) + "\t" + to_string(interface_pos);
		for (uword i = 0; i < model.charge_position.n_rows; ++i) {
			extrapolation_info += "\t" + to_string(model.charge_position(i, model.normal_direction) * model.cell_vectors_lengths(model.normal_direction));
		}
		log->debug(extrapolation_info);
		rowvec Es = zeros<rowvec>(extrapol_steps_num - 1), sizes = Es;
		tie(Es, sizes) = model.extrapolate(extrapol_steps_num, extrapol_steps_size);

		if (model.type == model_type::monolayer) {
			const rowvec3 unit_cell = model.cell_vectors_lengths / max(model.cell_vectors_lengths);
			const auto radius = 10.0;
			const auto ewald_shells = generate_shells(unit_cell, radius);
			const auto madelung_const = jellium_madelung_constant(ewald_shells, unit_cell, 1);
			auto madelung_term = -pow(model.total_charge, 2) * madelung_const / 2;
			nonlinear_fit_data fit_data = { Es ,sizes, madelung_term };
			const auto cs = nonlinear_fit(1e-10, fit_data);

			log->info("Madelung constant = " + ::to_string(madelung_const));
			const string fit_params = "c0= " + ::to_string(cs.at(0)) +
				", c1=" + ::to_string(cs.at(1)) +
				", c2=" + ::to_string(cs.at(2)) +
				", c3=" + ::to_string(cs.at(3));
			log->info("Non-linear fit parameters:" + fit_params);
			calculation_results.emplace_back("Non-linear fit parameters", fit_params);
			calculation_results.emplace_back("Madelung constant", ::to_string(madelung_const));

			E_isolated = cs.at(0) + (cs.at(1) - madelung_term) / cs.at(3);
			E_correction = E_isolated - EperModel0 - model.total_charge * dV;
		}
		else { // bulk and slab models
			const colvec pols = polyfit(sizes, Es, 1);
			const colvec evals = polyval(pols, sizes.t());
			const auto linearfit_MSE = accu(square(evals.t() - Es)) / Es.n_elem * 100;
			const rowvec slopes = diff(Es) / diff(sizes);
			const auto extrapol_error_periodic = abs(slopes(0) - slopes(slopes.n_elem - 1));
			log->debug("--------------------------------------------------------");
			log->debug("Linear fit: Eper(Model) = {}/scaling + {}", ::to_string(pols(0)), ::to_string(pols(1)));
			log->debug("Linear fit Root Mean Square Error: {}", ::to_string(sqrt(linearfit_MSE)));
			log->debug("Polyfit evaluated energies: {}", ::to_string(evals));
			log->debug("Linear fit error for the periodic model: {}", ::to_string(extrapol_error_periodic));

			if (extrapol_error_periodic > 0.05) {
				log->debug("Extrapolation energy slopes: {}", to_string(slopes));
				log->critical("The extrapolated energies are not scaling linearly as expected! "
				"The slab thickness may be too small for this extrapolation algorithm. "
				"For calculating the charge correction energy for the 2D models use \"2D_model = yes\" in the input file.");
				finalize_loggers();
				exit(1);
			}

			E_isolated = EperModel0 - pols(0);
			E_correction = -pols(0) - model.total_charge * dV;

		}
		log->info("E_isolated from extrapolation with {}x{} steps: {}", to_string(extrapol_steps_num), to_string(extrapol_steps_size), ::to_string(E_isolated));
	}
	else {
		if (model.type == model_type::monolayer) {
			E_isolated = model.Eiso_bessel();
			E_correction = E_isolated - EperModel0 - model.total_charge * dV;
			log->info("E_isolated from the Bessel expansion of the Poisson equation: {}", ::to_string(E_isolated));
		}
		else {
			// input parameter checking function must prevent this from happening!
			log->critical("There is no algorithm other than the extrapolation for E_isolated calculation of the slab models in this version of the slabcc!");
			finalize_loggers();
			exit(1);
		}

	}
	calculation_results.emplace_back("E_isolated of the model charge", ::to_string(E_isolated));

	log->info("Energy correction for the model charge (E_iso-E_per-q*dV=): {}", ::to_string(E_correction) );
	calculation_results.emplace_back("Energy correction for the model charge (E_iso-E_per-q*dV)", ::to_string(E_correction));
	log->flush();
	
	finalize_loggers();
	output_log->info("\n[Results]");
	for (const auto &i : calculation_results) { output_log->info("{} = {}", i.first, i.second); }
	output_log->flush();
	
	//making sure all the files are written
	for (auto &promise : future_files) { promise.get(); }

	log->trace("Calculations successfully ended!");
	return 0;
}