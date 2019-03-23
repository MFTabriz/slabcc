// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "stdafx.h"
#include "slabcc_core.hpp"
#include "isolated.hpp"
#include "vasp.hpp"
using namespace std;

int verbos = 0;
slabcc_cell_type slabcc_cell;
const double ang_to_bohr = 1e-10 / datum::a_0;  //1.88972612546
const double Hartree_to_eV = datum::R_inf * datum::h * datum::c_0 / datum::eV * 2; //27.2113860193

int main(int argc, char *argv[])
{

	string input_file = "slabcc.in";
	string output_file = "slabcc.out";
	string log_file = "slabcc.log";

	parse_cli(argc, argv, input_file, output_file, log_file);
	initialize_logger(log_file);
	ofstream output_fstream(output_file);
	output_fstream << setprecision(12);
	
	// default values should be defined in the parser function

	string CHGCAR_neutral = "";
	string LOCPOT_neutral = "";
	string LOCPOT_charged = "";
	string CHGCAR_charged = "";
	string opt_algo = "";			//optimization algorithm
	mat charge_position;			//center of each Gaussian model charge
	rowvec Qd;						//charge of each Gaussian model charge
	rowvec sigma;					//width of the model charges
	rowvec diel_in;					//diagonal elements of slab dielectric tensor
	rowvec diel_out;				//diagonal elements of vacuum dielectric tensor
	rowvec3 slabcenter;
	uword normal_direction = 0;		//index of the normal direction (0/1/2)
	rowvec2 interfaces;				//interfaces in relative coordinates
	double diel_erf_beta = 0;		//beta value of the erf for dielectric profile generation
	double opt_tol = 0;				//relative optimization tolerance
	double extrapol_grid_x = 0;		//extrapolation grid size multiplier
	double opt_grid_x = 0;			//optimization grid size multiplier
	int max_eval = 0;				//maximum number of steps for the optimization function evaluation
	int max_time = 0;				//maximum time for the optimization in minutes
	int extrapol_steps_num = 0;		//number of extrapolation steps for E_isolated calculation
	double extrapol_steps_size = 0; //size of each extrapolation step with respect to the initial supercell size
	bool optimize_charge = false;	//optimize the charge_position and sigma
	bool optimize_interfaces = false;//optimize the position of interfaces
	bool extrapolate = false;		//use the extrapolation for E-isolated calculations
	bool model_2D = false;		//the model is 2D
	
	// parameters read from the input file
	const input_data input_file_variables = {
		CHGCAR_neutral, LOCPOT_charged, LOCPOT_neutral, CHGCAR_charged,
		opt_algo, charge_position, Qd, sigma, slabcenter, diel_in, diel_out,
		normal_direction, interfaces, diel_erf_beta,
		opt_tol, optimize_charge, optimize_interfaces, extrapolate, model_2D, opt_grid_x,
		extrapol_grid_x, max_eval, max_time, extrapol_steps_num, extrapol_steps_size };

	parse_input_params(input_file, output_fstream, input_file_variables);
	auto log = spdlog::get("loggers");
	log->debug("SLABCC: version {}.{}.{}", SLABCC_VERSION_MAJOR, SLABCC_VERSION_MINOR, SLABCC_VERSION_PATCH);
	log->debug("Armadillo library: version {}.{}.{}", ARMA_VERSION_MAJOR, ARMA_VERSION_MINOR, ARMA_VERSION_PATCH);
	log->debug("NLOPT library: version {}.{}.{}", nlopt::version_major(), nlopt::version_minor(), nlopt::version_bugfix());
	log->debug("SPDLOG library: version {}.{}.{}", SPDLOG_VER_MAJOR, SPDLOG_VER_MINOR, SPDLOG_VER_PATCH);
	log->debug("SLABCC input file: {}", input_file);
	log->debug("SLABCC output file: {}", output_file);
	log->debug("SLABCC log file: {}", log_file);

	vector<pair<string, string>> calculation_results;

	check_inputs(input_file_variables);

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

	verify_cells(Neutral_supercell, Charged_supercell);

	//cell vectors of the CHGCAR and LOCPOT files (bohr)
	const mat33 cell_vectors = abs(Neutral_supercell.cell_vectors) * Neutral_supercell.scaling * ang_to_bohr;

	const urowvec3 input_grid_size = SizeVec(Neutral_supercell.charge);
	UpdateCell(cell_vectors, input_grid_size);

	//slabcc_cell volume in bohr^3 (only works in the orthogonal case!)
	const auto volume = prod(slabcc_cell.vec_lengths);

	const rowvec3 relative_shift = 0.5 - slabcenter;

	const rowvec3 rounded_relative_shift = round(input_grid_size % relative_shift) / input_grid_size;

	interfaces = fmod(interfaces + rounded_relative_shift(normal_direction), 1);
	shift_structure(Neutral_supercell, rounded_relative_shift);
	shift_structure(Charged_supercell, rounded_relative_shift);
	charge_position += repmat(rounded_relative_shift, charge_position.n_rows, 1);
	charge_position = fmod_p(charge_position, 1);

	log->debug("Slab normal direction index (0-2): {}", normal_direction);
	log->trace("Shift to center done!");

	supercell Defect_supercell = Neutral_supercell;
	Defect_supercell.potential = Charged_supercell.potential - Neutral_supercell.potential;
	Defect_supercell.charge = Charged_supercell.charge - Neutral_supercell.charge;

	if (is_active(verbosity::write_defect_file)) {
		future_files.push_back(async(launch::async, write_CHGPOT, "LOCPOT", "slabcc_D.LOCPOT", Defect_supercell));
		future_files.push_back(async(launch::async, write_CHGPOT, "CHGCAR", "slabcc_D.CHGCAR", Defect_supercell));
	}

	//normalize the charges and potentials
	Neutral_supercell.charge *= -1.0 / volume;
	Charged_supercell.charge *= -1.0 / volume;
	Defect_supercell.charge *= -1.0 / volume;
	Defect_supercell.potential *= -1.0;

	// total extra charge of the VASP calculation
	const auto Q0 = accu(Defect_supercell.charge) * slabcc_cell.voxel_vol;
	if (abs(Q0) < 0.001) {
		log->debug("Total extra charge: {}", Q0);
		log->warn("Total extra charge seems to be very small. Please make sure the path to the input CHGCAR files are set properly!");
	}
	// convert charge fractions to actual charges
	Qd *= Q0 / accu(Qd);

	mat dielectric_profiles = zeros<mat>(slabcc_cell.grid(normal_direction), 3);

	// model charge distribution (e/bohr^3), negative for presence of the electron 
	cx_cube rhoM(as_size(slabcc_cell.grid), fill::zeros);

	//potential resulted from the model charge (Hartree)
	cx_cube V(arma::size(rhoM));
	//difference of the potential resulted from the model charge and the QM calculation (VASP) results (eV)
	cube V_diff(arma::size(V));

	double potential_MSE = 0;

	double initial_potential_MSE = -1;

	// variables to optimize
	opt_vars opt_vars = { interfaces, sigma, Qd, charge_position };

	if (optimize_charge || optimize_interfaces) {
		const auto interfaces0 = interfaces;
		const auto charge_position0 = charge_position;

		const rowvec optimization_grid_x = regspace<rowvec>(1, 1.0 / opt_grid_x, slabcc_cell.grid(0));
		const rowvec optimization_grid_y = regspace<rowvec>(1, 1.0 / opt_grid_x, slabcc_cell.grid(1));
		const rowvec optimization_grid_z = regspace<rowvec>(1, 1.0 / opt_grid_x, slabcc_cell.grid(2));

		cube interpolated_potential = interp3(Defect_supercell.potential, optimization_grid_x, optimization_grid_y, optimization_grid_z);
		interpolated_potential -= accu(interpolated_potential) / interpolated_potential.n_elem;
		log->debug("Optimization grid size: " + to_string(SizeVec(interpolated_potential)));

		//data needed for potential error calculation
		opt_data optimize_data = { Q0, diel_erf_beta, diel_in, diel_out, interpolated_potential, dielectric_profiles, rhoM, V, V_diff, initial_potential_MSE };
		UpdateCell(cell_vectors, SizeVec(interpolated_potential));
		potential_MSE = do_optimize(opt_algo, opt_tol, max_eval, max_time, optimize_data, opt_vars, optimize_charge, optimize_interfaces);
		UpdateCell(cell_vectors, input_grid_size);
		//add back the last Gaussian charge (removed in the optimization)
		Qd(Qd.n_elem - 1) = Q0 - accu(Qd);

		//write the unshifted optimized values to the file
		output_fstream << "\n[Optimized_model_parameters]\n";
		if (optimize_interfaces) {
			const rowvec2 original_interfaces = fmod_p(interfaces - rounded_relative_shift(normal_direction), 1);
			output_fstream << "interfaces_optimized = " << original_interfaces << '\n';
			// need a clever way to check the change in the interface positions
			// the order may change and also the interfaces may move from one supercell to the next
			// if there is too much change, warn the user
		}

		if (optimize_charge) {
			const mat orig_charge_position = fmod_p(charge_position - repmat(rounded_relative_shift, charge_position.n_rows, 1), 1);
			if (Qd.n_elem > 1) {
				output_fstream << "charge_fraction_optimized =" << abs(Qd / accu(Qd));
			}
			output_fstream << "charge_sigma_optimized = " << sigma << '\n';
			output_fstream << "charge_position_optimized = " << orig_charge_position << '\n';


			const mat charge_position_change = abs(charge_position0 - charge_position);
			if (charge_position_change.max() > 0.1) {
				log->warn("The optimized position for the extra charge is significantly different from the initial value.");
				log->warn("Please make sure that the final position of the extra charge have been estimated correctly!");
				log->trace("Charge position changes: ", to_string(charge_position_change));
			}
		}

		output_fstream.flush();

		check_inputs(input_file_variables);

		if (initial_potential_MSE * (opt_tol + 1) < potential_MSE) {
			// Don't panic! either NLOPT seems to be malfunctioning
			// or we are not correctly logging/checking the result
			log->critical("Optimization failed!");
			log->critical("Potential error of the initial parameters seems to be smaller than the optimized parameters!");
			log->critical("You may want to change the initial guess for charge_position, change the optimization algorithm, or turn off the optimization.");
			exit(1);
		}
	}
	opt_data optimized_data = { Q0, diel_erf_beta, diel_in, diel_out, Defect_supercell.potential, dielectric_profiles, rhoM, V, V_diff, initial_potential_MSE };
	auto local_param = optimizer_packer(opt_vars);
	vector<double> gradients = {};
	potential_MSE = potential_eval(get<0>(local_param), gradients, &optimized_data);

	if (potential_MSE > 0.01) {
		log->warn("MSE of the model charge potential is large. The calculated correction energies may not be accurate!");
	}

	const auto V_error_x = conv_to<rowvec>::from(planar_average(0, V_diff));
	const auto V_error_y = conv_to<rowvec>::from(planar_average(1, V_diff));
	const auto V_error_z = conv_to<rowvec>::from(planar_average(2, V_diff));

	//potential error in each direction
	rowvec3 V_error_planars = { accu(square(V_error_x)), accu(square(V_error_y)), accu(square(V_error_z)) };
	V_error_planars = sqrt(V_error_planars/V_diff.n_elem);

	log->debug("Directional RMSE: " + to_string(V_error_planars));
	if (max(V_error_planars) / min(V_error_planars) > 10) {
		log->warn("The potential error is highly anisotropic.");
		log->warn("If the potential error is large, this usually means that either the extra charge is not properly described by the model Gaussian charge "
			"or the chosen dielectric tensor is not a good representation of the actual tensor!");
	}

	log->debug("Potential error anisotropy: {}", max(V_error_planars) / min(V_error_planars));
	log->debug("Cell dimensions (bohr): " + to_string(slabcc_cell.vec_lengths));
	log->debug("Grid size: " + to_string(input_grid_size));
	log->debug("Volume (bohr^3): {}", volume);


	if (is_active(verbosity::write_defect_file)) {
		supercell Model_supercell = Neutral_supercell;
		//charge is normalized to the VASP CHGCAR convention (rho * Vol)
		//Also, positive value for the electron charge! (the probability of finding an electron)
		Model_supercell.charge = -real(rhoM) * slabcc_cell.voxel_vol * rhoM.n_elem;
		Model_supercell.potential = -real(V) * Hartree_to_eV;
		future_files.push_back(async(launch::async, write_CHGPOT, "CHGCAR", "slabcc_M.CHGCAR", Model_supercell));
		future_files.push_back(async(launch::async, write_CHGPOT, "LOCPOT", "slabcc_M.LOCPOT", Model_supercell));
	}

	if (is_active(verbosity::write_dielectric_file)) {
		write_mat2file(dielectric_profiles, "slabcc_DIEL.dat");
	}
	if (is_active(verbosity::write_planarAvg_file)) {
		write_planar_avg(Neutral_supercell.potential, Neutral_supercell.charge * slabcc_cell.voxel_vol, "N");
		write_planar_avg(Charged_supercell.potential, Charged_supercell.charge * slabcc_cell.voxel_vol, "C");
		write_planar_avg(Defect_supercell.potential, Defect_supercell.charge * slabcc_cell.voxel_vol, "D");
		write_planar_avg(real(V) * Hartree_to_eV, real(rhoM) * slabcc_cell.voxel_vol, "M");
	}
	else if (is_active(verbosity::write_normal_planarAvg)) {
		write_planar_avg(Defect_supercell.potential, Defect_supercell.charge * slabcc_cell.voxel_vol, "D", slabcc_cell.normal_direction);
		write_planar_avg(real(V) * Hartree_to_eV, real(rhoM) * slabcc_cell.voxel_vol, "M", slabcc_cell.normal_direction);
	}
	//total charge of the model
	const auto Q = accu(real(rhoM)) * slabcc_cell.voxel_vol;
	//add jellium to the charge (Because the V is normalized, it is not needed in solving the Poisson eq. but it is needed in the energy calculations)
	rhoM -= Q / volume;
	uword farthest_element_index = Q < 0 ? real(V).index_max() : real(V).index_min();

	const auto dV = V_diff(farthest_element_index);
	log->info("Potential alignment (dV=): " + ::to_string(dV));
	calculation_results.emplace_back("dV", ::to_string(dV));

	if (abs(dV) > 0.05) {
		log->warn("The potential alignment term (dV) is relatively large." );
		log->warn("The constructed model may not be accurate!");
	}

	log->debug("Alignment term calculation grid point: " + to_string(ind2sub(as_size(slabcc_cell.grid), farthest_element_index)));

	const auto EperModel0 = 0.5 * accu(real(V) % real(rhoM)) * slabcc_cell.voxel_vol * Hartree_to_eV;
	log->info("E_periodic of the model charge: " + ::to_string(EperModel0));
	calculation_results.emplace_back("E_periodic of the model charge", ::to_string(EperModel0));

	const rowvec max_size = slabcc_cell.vec_lengths * (1 + extrapol_steps_size * (extrapol_steps_num - 1));
	if (min(slabcc_cell.grid * extrapol_grid_x / max_size) < 1) {
		log->warn("The extrapolation grid is very coarse! The extrapolation energies for the large model charges may not be accurate.");
		log->warn("The energy of the largest extrapolated model will be calculated with {} points/bohr grid", min(slabcc_cell.grid * extrapol_grid_x / max_size));
		log->warn("You should increase the extrapolation grid multiplier or decrease the number/size of extrapolation steps.");
	}

	log->debug("Difference of the charge in the input files: " + ::to_string( Q0));
	log->debug("Total charge of the model: " + ::to_string(Q));
	if (abs(Q - Q0) > 0.05) {
		log->warn("Part of the extra charge is missing from the model. "
			"This usually happens when the size of the supercell is too small and cannot contain the whole "
			"Gaussian charge. Otherwise, the width of the Gaussian charge may be too large. "
			"The present charge correction method is not suitable for these cases!");
	}
	double E_isolated = 0;
	double E_correction = 0;
	if (extrapolate) {
		const rowvec3 extrapolation_grid_size = extrapol_grid_x * conv_to<rowvec>::from(slabcc_cell.grid);
		const urowvec3 extrapolation_grid = { (uword)extrapolation_grid_size(0), (uword)extrapolation_grid_size(1), (uword)extrapolation_grid_size(2) };
		log->debug("Extrapolation grid size: " + to_string(extrapolation_grid));
		log->debug("--------------------------------------------------------");
		log->debug("Scaling\tE_periodic\t\tinterfaces\t\tcharge position");
		const rowvec2 interface_pos = interfaces * slabcc_cell.vec_lengths(slabcc_cell.normal_direction);
		string extrapolation_info = to_string(1.0) + "\t" + ::to_string(EperModel0) + "\t" + to_string(interface_pos);
		for (auto i = 0; i < charge_position.n_rows; ++i) {
			extrapolation_info += "\t" + to_string(charge_position(i, slabcc_cell.normal_direction) * slabcc_cell.vec_lengths(slabcc_cell.normal_direction));
		}
		log->debug(extrapolation_info);
		rowvec Es = zeros<rowvec>(extrapol_steps_num), sizes = Es;

		if (model_2D) {
			tie(Es, sizes) = extrapolate_2D(extrapol_steps_num, extrapol_steps_size, diel_in, diel_out,
				interfaces, diel_erf_beta, charge_position, Qd, sigma, extrapol_grid_x);
			const rowvec3 unit_cell = slabcc_cell.vec_lengths / max(slabcc_cell.vec_lengths);
			const auto radius = 10.0;
			const auto ewald_shells = generate_shells(unit_cell, radius);
			const auto madelung_const = jellium_madelung_constant(ewald_shells, unit_cell, 1);
			auto madelung_term = -pow(Q, 2) * madelung_const / 2;
			nonlinear_fit_data fit_data = { Es ,sizes, madelung_term };
			const auto cs = nonlinear_fit(1e-12, fit_data);

			log->info("Madelung constant = " + ::to_string(madelung_const));
			const string fit_params = "c0= " + ::to_string(cs.at(0)) +
				", c1=" + ::to_string(cs.at(1)) +
				", c2=" + ::to_string(cs.at(2)) +
				", c3=" + ::to_string(cs.at(3));
			log->info("Non-linear fit parameters:" + fit_params);
			calculation_results.emplace_back("Non-linear fit parameters", fit_params);
			calculation_results.emplace_back("Madelung constant", ::to_string(madelung_const));

			E_isolated = cs.at(0) + (cs.at(1) - madelung_term) / cs.at(3);
			E_correction = E_isolated - EperModel0 - Q * dV;
		}
		else {
			tie(Es, sizes) = extrapolate_3D(extrapol_steps_num, extrapol_steps_size, diel_in, diel_out,
				interfaces, diel_erf_beta, charge_position, Qd, sigma, extrapol_grid_x);

			const colvec pols = polyfit(sizes, Es, 1);
			const colvec evals = polyval(pols, sizes.t());
			const auto linearfit_MSE = accu(square(evals.t() - Es)) / Es.n_elem * 100;
			const rowvec slopes = diff(Es) / diff(sizes);
			const auto extrapol_error_periodic = abs(slopes(0) - slopes(slopes.n_elem - 1));
			log->debug("--------------------------------------------------------");
			log->debug("Linear fit: Eper(Model) = " + ::to_string(pols(0)) + "/scaling + " + ::to_string(pols(1)));
			log->debug("Linear fit Root Mean Square Error: " + ::to_string(sqrt(linearfit_MSE)));
			log->debug("Polyfit evaluated energies: " + ::to_string(evals));
			log->debug("Linear fit error for the periodic model: " + ::to_string(extrapol_error_periodic));

			if (extrapol_error_periodic > 0.05) {
				log->error("The extrapolated energies are not scaling linearly!");
				log->error("The extrapolation steps/size may be too large for the grid size or the slab thickness in too small for this extrapolation algorithm");
				log->error("You may need to use larger \"extrapolate_grid_x\" or smaller \"extrapolate_steps_size\"/\"extrapolate_steps_number\"");
				log->error("For calculating the charge correction energy for the 2D models use \"2D_model = yes\" in the input file");
				log->error("Extrapolation energy slopes: " + to_string(slopes));
			}

			E_isolated = EperModel0 - pols(0);
			E_correction = -pols(0) - Q * dV;

		}
		log->info("E_isolated from extrapolation with " + to_string(extrapol_steps_num) + "x" + to_string(extrapol_steps_size) + " steps: " + ::to_string(E_isolated));
		calculation_results.emplace_back("E_isolated of the model charge", ::to_string(E_isolated));

	}
	else {
		if (model_2D) {
			// TODO: implement multiple charges
			const mat dielectric_profiles = dielectric_profiles_gen(interfaces, diel_in, diel_out, diel_erf_beta);
			E_isolated = Eiso_bessel(Qd(0), charge_position(0, normal_direction) * slabcc_cell.vec_lengths(normal_direction), sigma(0), dielectric_profiles);
			E_correction = E_isolated - EperModel0 - Q * dV;
			log->info("E_isolated from the Bessel expansion of the Poisson equation: " + ::to_string(E_isolated));
			calculation_results.emplace_back("E_isolated of the model charge", ::to_string(E_isolated));

		}
		else {
			log->error("There is no algorithm other than extrapolation for E_isolated calculation of the slab models in this version of the slabcc!");
		}

	}



	log->info("Energy correction for the model charge (E_iso-E_per-q*dV=): " + ::to_string(E_correction) );
	calculation_results.emplace_back("Energy correction for the model charge (E_iso-E_per-q*dV)", ::to_string(E_correction));
	log->flush();
	const string msgs_file = "MSG";
	ifstream messages_list(msgs_file);
	if (!file_is_empty(messages_list)) {
		output_fstream << "\n[Messages]\n";
		output_fstream << messages_list.rdbuf();
	}
	messages_list.close();
	
	output_fstream << "\n[Results]\n";
	for (const auto &i : calculation_results) { output_fstream << i.first << " = " << i.second << '\n'; }
	output_fstream.close();
	
	//making sure all the files are written
	for (auto &promise : future_files) { promise.get(); }

	log->trace("Calculations successfully ended!");
	log->~logger();
	remove(msgs_file.c_str());
	return 0;
}
