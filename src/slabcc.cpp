// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "stdafx.h"
#include "slabcc_core.hpp"
#include "vasp.hpp"
using namespace std;

//verbosity level to be compared with "enum verbosity" by is_active()
int verbos = 0;
slabcc_cell_type slabcc_cell;			//All the parameters various functions need to be aware of
const double ang_to_bohr = 1.889725989;
const double Hartree_to_eV = 27.2113714880369;
const float version = 0.3F;
chrono::time_point<chrono::steady_clock> t0; //initial start time of the program


int main(int argc, char *argv[])
{

	t0 = chrono::steady_clock::now();

	string input_file = "slabcc.in";
	string output_file = "slabcc.out";

	parse_cli(argc, argv, input_file, output_file);
	
	ofstream output_fstream;
	output_fstream.open(output_file);
	output_fstream << setprecision(12);

	// default values should be defined in the parser function

	string CHGCAR_NEU = "";			//charge density of neutral system
	string LOCPOT_NEU = "";			//potential of neutral system
	string LOCPOT_CHG = "";			//potential of charged system
	string CHGCAR_CHG = "";			//charge density of charged system
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
	int max_eval = 0;				//maximum number of steps for the optimization function evaluattion
	int max_time = 0;				//maximum time for the optimization in minutes
	int extrapol_steps_num = 0;		//number of extrapolation steps for E_isolated calculation
	double extrapol_steps_size = 0; //size of each extrapolation step with respect to the initial supercell size
	bool optimize_charge = false;	//optimize the charge_position and sigma
	bool optimize_interfaces = false;//optimize the position of interfaces
	bool extrapol_slab = true;		//extrapolate the slab thickness or not!

	// parameters read from the input file
	const input_data input_file_variables = {
		CHGCAR_NEU, LOCPOT_CHG, LOCPOT_NEU, CHGCAR_CHG, 
		opt_algo, charge_position, Qd, sigma, slabcenter, diel_in, diel_out,
		normal_direction, interfaces, diel_erf_beta,
		opt_tol, optimize_charge, optimize_interfaces, extrapol_slab, opt_grid_x, extrapol_grid_x, max_eval, max_time, extrapol_steps_num, extrapol_steps_size };

	parse_input_params(input_file, output_fstream, input_file_variables);
	
	if (is_active(verbosity::more_digits)) {
		cout << setprecision(12);
	}

	//calculation results to be written to the slabcc output file
	vector<pair<string, string>> results;

	check_inputs(input_file_variables);

	//promises for async read of CHGCAR and POTCAR files
	vector<future<cube>> future_cells;

	//promises for async file writing (can be replaced by a deque if the number of files increases)
	vector<future<void>> future_files;

	future_cells.push_back(async(launch::async, read_CHGPOT, CHGCAR_NEU));
	future_cells.push_back(async(launch::async, read_CHGPOT, CHGCAR_CHG));
	future_cells.push_back(async(launch::async, read_CHGPOT, LOCPOT_NEU));
	future_cells.push_back(async(launch::async, read_CHGPOT, LOCPOT_CHG));

	if (is_active(verbosity::detailed_progress)) {
		cout << timing() << "Slab normal direction index: " << normal_direction << endl;
		cout << timing() << "Started reading CHGCAR and LOCPOT files" << endl;
	}

	supercell Neutral_supercell = read_POSCAR(CHGCAR_NEU);
	supercell Charged_supercell = read_POSCAR(CHGCAR_CHG);

	Neutral_supercell.charge = future_cells.at(0).get();
	Charged_supercell.charge = future_cells.at(1).get();
	Neutral_supercell.potential = future_cells.at(2).get();
	Charged_supercell.potential = future_cells.at(3).get();

	check_cells(Neutral_supercell, Charged_supercell, input_file_variables);

	//volume of supercell models (A^3)
	const double volume_a3 = det(Neutral_supercell.cell_vectors * Neutral_supercell.scaling);

	Neutral_supercell.charge /= volume_a3;
	Charged_supercell.charge /= volume_a3;

	//lengths of the cell vectors of the CHGCAR and LOCPOT files in bohr
	const rowvec3 cell = diagvec(Neutral_supercell.cell_vectors).t() * Neutral_supercell.scaling * ang_to_bohr;

	//grid density of the CHGCAR and LOCPOT files
	const urowvec3 grid = SizeVec(Neutral_supercell.charge);
	UpdateCell(cell, grid);

	const rowvec3 relative_shift = 0.5 - slabcenter;

	//relative shifts rounded to the nearest grid point
	const rowvec3 rounded_relative_shift = round(grid % relative_shift) / grid;

	interfaces = fmod(interfaces + rounded_relative_shift(normal_direction), 1);
	shift_structure(Neutral_supercell, rounded_relative_shift);
	shift_structure(Charged_supercell, rounded_relative_shift);
	charge_position += repmat(rounded_relative_shift, charge_position.n_rows, 1);
	charge_position = fmod_p(charge_position, 1);

	if (is_active(verbosity::detailed_progress)) {
		cout << timing() << "Shift to center done!" << endl;
	}

	supercell Defect_supercell = Neutral_supercell;
	Defect_supercell.potential = Charged_supercell.potential - Neutral_supercell.potential;
	Defect_supercell.charge = Charged_supercell.charge - Neutral_supercell.charge;

	if (is_active(verbosity::write_defect_file)) {
		future_files.push_back(async(launch::async, write_CHGPOT, "LOCPOT", "slabcc_D.LOCPOT", Defect_supercell));
		cout << timing() << "Started writing slabcc_D.LOCPOT" << endl;
		future_files.push_back(async(launch::async, write_CHGPOT, "CHGCAR", "slabcc_D.CHGCAR", Defect_supercell));
		cout << timing() << "Started writing slabcc_D.CHGCAR" << endl;
	}

	//normalize the charges and potentials
	Neutral_supercell.charge *= -volume_a3 / Neutral_supercell.charge.n_elem;
	Charged_supercell.charge *= -volume_a3 / Charged_supercell.charge.n_elem;
	Defect_supercell.charge *= -volume_a3 / Defect_supercell.charge.n_elem;
	Defect_supercell.potential *= -1.0;

	// total extra charge of the VASP calculation
	const double Q0 = accu(Defect_supercell.charge);
	// normalize the value of the charges (is necessary if not provided or relatively defined)
	Qd *= Q0 / accu(Qd);

	if (is_active(verbosity::detailed_progress)) {
		cout << timing() << "Extra charge evaluated!" << endl;
	}
	//dielectric profile
	mat diels = zeros<mat>(slabcc_cell.grid(normal_direction), 3);

	// model charge distribution
	cx_cube rhoM(as_size(slabcc_cell.grid), fill::zeros);

	//potential resulted from the model charge
	cx_cube V(arma::size(rhoM));
	//difference of the potential resulted from the model charge and the QM calculation (VASP) results
	cube V_diff(arma::size(V));
	//mean squared potential error in 3D
	double pot_MSE = 0;

	//potential MSE returned from the 1st step of evaluation
	double initial_pot_MSE = -1;
		
	// variables to optimize
	opt_vars opt_vars = { interfaces, sigma, Qd, charge_position };

	if (optimize_charge || optimize_interfaces) {
		//interpolation grid
		rowvec vx = regspace<rowvec>(1, 1.0 / opt_grid_x, slabcc_cell.grid(0));
		rowvec vy = regspace<rowvec>(1, 1.0 / opt_grid_x, slabcc_cell.grid(1));
		rowvec vz = regspace<rowvec>(1, 1.0 / opt_grid_x, slabcc_cell.grid(2));
		//interpolated defect potential for the optimization process
		cube interpolated_potential = interp3(Defect_supercell.potential, vx, vy, vz);
		interpolated_potential -= accu(interpolated_potential) / interpolated_potential.n_elem;
		if (is_active(verbosity::detailed_progress)) {
			cout << timing() << "Optimization grid size:" << SizeVec(interpolated_potential) << endl;
		}
		
		//data needed for potential error calculation
		opt_data optimize_data = { Q0, diel_erf_beta, diel_in, diel_out, interpolated_potential, diels, rhoM, V, V_diff, initial_pot_MSE };
		UpdateCell(cell, SizeVec(interpolated_potential));
		pot_MSE = do_optimize(opt_algo, opt_tol, max_eval, max_time, optimize_data, opt_vars, optimize_charge, optimize_interfaces);
		UpdateCell(cell, grid);

		//write the unshifted optimized values to the file
		output_fstream << endl << "[Optimized_parameters]" << endl;
		if (optimize_interfaces) {
			const rowvec2 original_interfaces = fmod_p(interfaces - rounded_relative_shift(normal_direction), 1);
			output_fstream << "interfaces_optimized = " << original_interfaces << endl;
		}
		
		if (optimize_charge) {
			const mat orig_charge_position = fmod_p(charge_position - repmat(rounded_relative_shift, charge_position.n_rows, 1), 1);
			if (Qd.n_elem > 1) {
				output_fstream << "charge_fraction_optimized =" << abs(Qd / accu(Qd));
			}
			output_fstream << "charge_sigma_optimized =" << sigma << endl;
			output_fstream << "charge_position_optimized =" << orig_charge_position << endl;
		}

		output_fstream.flush();

		//make sure the parameters after the optimization still make sense
		check_inputs(input_file_variables);

		if (initial_pot_MSE * (opt_tol + 1) < pot_MSE) {
			// Don't panic! either NLOPT seems to be malfunctioning
			// or we are not correctly logging/checking the result
			cout << endl << "ERROR : Optimization failed!" << endl;
			cout << "Potential error of the initial parameters seems to be smaller than the optimized parameters!" << endl;
			cout << "You may want to change the initial guess for charge_position, change the optimization algorithm, or turn off the optimization." << endl << endl;
			exit(1);
		}
	}
	
	opt_data optimized_data = { Q0, diel_erf_beta, diel_in, diel_out, Defect_supercell.potential, diels, rhoM, V, V_diff, initial_pot_MSE };
	auto local_param = optimizer_packer(opt_vars, true, true);
	vector<double> grads = {};
	pot_MSE = potential_eval(get<0>(local_param), grads, &optimized_data);

	if (pot_MSE > 1) {
		cout << endl << timing() << ">> WARNING <<: MSE of the model charge potential is large. The calculated correction energies may not be accurate!" << endl << endl;
	}

	const auto V_error_x = conv_to<rowvec>::from(planar_average(0, V_diff));
	const auto V_error_y = conv_to<rowvec>::from(planar_average(1, V_diff));
	const auto V_error_z = conv_to<rowvec>::from(planar_average(2, V_diff));

	//potential error in each direction
	rowvec3 V_error_planars = { accu(square(V_error_x)), accu(square(V_error_y)), accu(square(V_error_z)) };
	V_error_planars /= V_diff.n_elem;

	if (is_active(verbosity::detailed_progress)) {
		cout << timing() << "Directional relative potential errors: " << V_error_planars << endl;
	}
	if (max(V_error_planars) / min(V_error_planars) > 10) {
		cout << endl << timing() << ">> WARNING <<: The potential error is highly anisotropic. Either the extra charge is not properly described by the model Gaussian charge or the choosen dielectric tensor is not a good representation of the actual tensor!" << endl << endl;
	}

	//slabcc_cell volume in bohr^3
	const double volume = prod(slabcc_cell.lengths);

	if (is_active(verbosity::steps)) {
		cout << timing() << "Potential error anisotropy: " << max(V_error_planars) / min(V_error_planars) <<endl;
		cout << timing() << "Shifted slab interfaces:" << interfaces << endl;
		cout << timing() << "Shifted charges position (bohr):" << mat(charge_position.each_row() % cell) << endl;
		cout << timing() << "Cell dimensions (bohr):" << cell << endl;
		cout << timing() << "Grid size:" << grid << endl;
		cout << timing() << "Volume (bohr^3): " << volume << endl;
	}

	if (is_active(verbosity::write_defect_file)) {
		supercell MDout = Neutral_supercell;
		MDout.charge = real(rhoM);
		MDout.potential = real(V) * Hartree_to_eV;
		future_files.push_back(async(launch::async, write_CHGPOT, "CHGCAR", "slabcc_M.CHGCAR", MDout));
		future_files.push_back(async(launch::async, write_CHGPOT, "LOCPOT", "slabcc_M.LOCPOT", MDout));
		cout << timing() << "Started writing slabcc_M files" << endl;
	}

	if (is_active(verbosity::write_dielectric_file)) {
		write_mat2file(diels, "slabcc_DIEL.dat");
	}
	if (is_active(verbosity::write_planarAvg_file)) {
		write_planar_avg(Neutral_supercell, "N");
		write_planar_avg(Charged_supercell, "C");
		write_planar_avg(Defect_supercell, "D");
		write_planar_avg(V * Hartree_to_eV, rhoM * slabcc_cell.voxel_vol, "M");
		if (is_active(verbosity::detailed_progress)) {
			cout << timing() << "Planar averages are written!" << endl;
		}
	}
	//total charge of the model
	double Q = accu(real(rhoM)) * slabcc_cell.voxel_vol;
	//add jelium to the charge (Because the V is normalized, it is not needed in solving the Poisson eq. but in energy calculation it is needed)
	rhoM -= Q / volume;
	//index of the element least affected by the extra charge
	uword index_far = 0;
	if (Q < 0) {
		index_far = real(V).index_max();
	}
	else {
		index_far = real(V).index_min();
	}
	auto dV = V_diff(index_far);
	cout << timing() << "Potential alignment (dV=): " << dV << endl;
	results.emplace_back("dV", INIReader::to_string(dV));

	if (dV > 0.05) {
		cout << endl << timing() << ">> WARNING <<: The potential alignment term (dV) is relatively large." << endl;
		cout << "The constructed model may not be accurate!" << endl << endl;
	}

	if (is_active(verbosity::detailed_progress)) {
		cout << timing() << "Alignment term calculated at grid point: " << ind2sub(as_size(slabcc_cell.grid), index_far) << endl;
	}

	//normalized charge of the original defect
	cube rho0 = Defect_supercell.charge / slabcc_cell.voxel_vol;
	double EperModel0 = 0.5 * accu(real(V) % real(rhoM)) * slabcc_cell.voxel_vol * Hartree_to_eV;
	cout << timing() << "E periodic of model charge: " << EperModel0 << endl;
	results.emplace_back("E periodic of model charge", INIReader::to_string(EperModel0));

	const rowvec max_size = slabcc_cell.lengths * (1 + extrapol_steps_size * (extrapol_steps_num - 1));
	if (min(slabcc_cell.grid * extrapol_grid_x / max_size) < 1) {
		cout << endl << timing() << ">> WARNING <<: The extrapolation grid is very coarse! The extrapolation energies for the large model charges may not be accurate." << endl;
		cout << timing() << "The energy of the largest extrapolated model will be calculated with "  << min(slabcc_cell.grid * extrapol_grid_x / max_size) << " points/bohr grid" << endl;
		cout << timing() << "You should increase the extrapolation grid multiplier or decrese the number/size of extrapolation steps." << endl << endl;
	}

	if (is_active(verbosity::steps)) {
		cout << timing() << "Total extra charge: " << Q0 << endl;
		cout << timing() << "Scaling\tE_periodic\t\tinterfaces\t\tcharge position" << endl;
		cout << timing() << "1\t\t" << EperModel0 << "\t" << interfaces(0) * slabcc_cell.lengths(slabcc_cell.normal_direction) << "\t" << interfaces(1) * slabcc_cell.lengths(slabcc_cell.normal_direction) << endl;
	}
	rowvec Es = zeros<rowvec>(extrapol_steps_num), sizes = Es;
	double E_isolated = 0;
	double E_correction = 0;
	if (extrapol_slab) {
		tie(Es, sizes) = extrapolate_3D(extrapol_steps_num, extrapol_steps_size, diel_in, diel_out, interfaces, diel_erf_beta, charge_position, Qd, sigma, extrapol_grid_x);

		const colvec pols = polyfit(sizes, Es, 1);
		const colvec evals = polyval(pols, sizes.t());
		const double linearfit_MSE = accu(square(evals.t() - Es)) / Es.n_elem * 100;
		const rowvec slopes = diff(Es) / diff(sizes);
		const double extrapol_error_periodic = abs(slopes(0) - slopes(slopes.n_elem - 1));

		if (is_active(verbosity::detailed_progress)) {
			cout << timing() << "Linear fit: Eper(Model) = " << pols(0) << "/scaling + " << pols(1) << endl;
			cout << timing() << "Linear fit Mean Squared Error: " << linearfit_MSE << " %" << endl;
			cout << timing() << "Polyfit evaluated results: " << evals.t();
			cout << timing() << "Linear fit error for the periodic model: " << extrapol_error_periodic << endl;
		}

		if (extrapol_error_periodic > 0.05) {
			cout << endl << timing() << ">> ERROR <<: The extrapolated energies are not scaling linearly!" << endl;
			cout << timing() << "The extrapolation steps/size may be too large for the grid size or the slab thickness in too small for this extrapolation algorithm" << endl;
			cout << timing() << "You may need to use larger \"extrapolate_grid_x\" or smaller \"extrapolate_steps_size\"/\"extrapolate_steps_number\"" << endl << endl;
			cout << timing() << "Extrapolation energy slopes:" << slopes << endl;
		}

		E_isolated = EperModel0 - pols(0);
		E_correction = -pols(0) - Q * dV;

	}
	else {
		tie(Es, sizes) = extrapolate_2D(extrapol_steps_num, extrapol_steps_size, diel_in, diel_out, interfaces, diel_erf_beta, charge_position, Qd, sigma, extrapol_grid_x);
		rowvec3 unit_cell = cell / max(cell);
		const double radius = 10;
		auto ewald_shells = generate_shells(unit_cell, radius);
		const double madelung_const = jellium_madelung_constant(ewald_shells, unit_cell, 1);
		double madelung_term = -pow(Q, 2) * madelung_const / 2;
		nonlinear_fit_data fit_data = { Es ,sizes, madelung_term };
		auto cs = nonlinear_fit(1e-12, fit_data);

		if (is_active(verbosity::detailed_progress)) {
			cout << timing() << "Madelung constant = " << madelung_const << endl;
			cout << timing() << "Non-linear fit parameters: c0=" << cs.at(0) << ", c1=" << cs.at(1) << ", c2=" << cs.at(2) << ", c3=" << cs.at(3) << endl;
			const string fit_params = "c0= " + INIReader::to_string(cs.at(0)) + ", c1=" + INIReader::to_string(cs.at(1)) + ", c2=" + INIReader::to_string(cs.at(2)) + ", c3=" + INIReader::to_string(cs.at(3));
			results.emplace_back("Non-linear fit parameters", fit_params);
			results.emplace_back("Madelung constant", INIReader::to_string(madelung_const));
		}

		E_isolated = cs.at(0) + (cs.at(1) - madelung_term) / cs.at(3);
		E_correction = E_isolated - EperModel0 - Q * dV;
	}
	

	cout << timing() << "E_iso from extrapolation with " << extrapol_steps_num << "x" << extrapol_steps_size << " steps: " << E_isolated << endl;
	results.emplace_back("E isolated of model charge", INIReader::to_string(E_isolated));

	cout << timing() << "Energy correction for model charge (E_iso-E_per-q*dV=): " << E_correction << endl;
	results.emplace_back("Energy correction for model charge (E_iso-E_per-q*dV)", INIReader::to_string(E_correction));

	//write the results to output file
	output_fstream << endl << "[Results]" << endl;
	for (const auto &i : results) { output_fstream << i.first << " = " << i.second << endl; }
	output_fstream.close();

	//making sure all the files are written
	for (auto &k : future_files) { k.get(); }

	if (is_active(verbosity::detailed_progress)) {
		cout << timing() << "slabcc version " << setprecision(3) << version << " calculations successfully ended!" << endl;
	}

	return 0;
}
