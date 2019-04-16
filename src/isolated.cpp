// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "isolated.hpp"

tuple <rowvec, rowvec> extrapolate_3D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &charge_q, const mat &charge_sigma, const double &grid_multiplier, const bool &trivariate) {
	auto log = spdlog::get("loggers");
	const uword normal_direction = slabcc_cell.normal_direction;
	rowvec Es = zeros<rowvec>(extrapol_steps_num - 1), sizes = Es;
	const mat33 cell_vectors0 = slabcc_cell.vectors;
	const urowvec grid0 = slabcc_cell.grid;
	const rowvec3 extrapolation_grid_size = grid_multiplier * conv_to<rowvec>::from(grid0);
	const urowvec3 extrapolation_grid = { (uword)extrapolation_grid_size(0),(uword)extrapolation_grid_size(1),(uword)extrapolation_grid_size(2) };

	for (auto n = 0; n < extrapol_steps_num - 1; ++n) {

		const double extrapol_factor = extrapol_steps_size * (1.0 + n) + 1;

		UpdateCell(cell_vectors0 * extrapol_factor, extrapolation_grid);
		rowvec2 interfaces_ext = interfaces;
		mat charge_position_shifted = charge_position / extrapol_factor;

		// not a bulk model!
		if (!approx_equal(diel_in, diel_out, "absdiff", 0.02)) {
			
			const uvec interface_sorted_i = sort_index(interfaces);
			interfaces_ext(interface_sorted_i(1)) += abs(interfaces(0) - interfaces(1)) * (extrapol_factor - 1);
			interfaces_ext /= extrapol_factor;

			//charges moved to the same distance from their original nearest interface
			for (auto charge = 0; charge < charge_position.n_rows; ++charge) {
				const rowvec2 distance_to_interfaces = abs(charge_position(charge, normal_direction) - interfaces);
				if (distance_to_interfaces(0) < distance_to_interfaces(1)) {
					charge_position_shifted(charge, normal_direction) += interfaces_ext(0) - interfaces(0) / extrapol_factor;
				}
				else {
					charge_position_shifted(charge, normal_direction) += interfaces_ext(1) - interfaces(1) / extrapol_factor;
				}
			}
		}

		const mat dielectric_profiles = dielectric_profiles_gen(interfaces_ext, diel_in, diel_out, diel_erf_beta);

		cx_cube rhoM(as_size(slabcc_cell.grid), fill::zeros);
		for (uword i = 0; i < charge_position.n_rows; ++i) {
			rhoM += gaussian_charge(charge_q(i), charge_position_shifted.row(i).t(), charge_sigma.row(i), trivariate);
		}
		const double Q = accu(real(rhoM)) * slabcc_cell.voxel_vol;
		// (only works for the orthogonal cells!)
		rhoM -= Q / prod(slabcc_cell.vec_lengths);
		const auto V = poisson_solver_3D(rhoM, dielectric_profiles);
		const auto EperModel = 0.5 * accu(real(V % rhoM)) * slabcc_cell.voxel_vol * Hartree_to_eV;
		const rowvec2 interface_pos = interfaces_ext * slabcc_cell.vec_lengths(normal_direction);
		string extrapolation_info = to_string(extrapol_factor) + "\t" + ::to_string(EperModel) + "\t" + to_string(interface_pos);
		for (auto i = 0; i < charge_position_shifted.n_rows; ++i) {
			extrapolation_info += "\t" + to_string(charge_position_shifted(i, slabcc_cell.normal_direction) * slabcc_cell.vec_lengths(slabcc_cell.normal_direction));
		}
		log->debug(extrapolation_info);
		Es(n) = EperModel;
		sizes(n) = 1.0 / extrapol_factor;

	}

	return make_tuple(Es, sizes);
}

tuple <rowvec, rowvec> extrapolate_2D(const int &extrapol_steps_num, const double &extrapol_steps_size, const rowvec3 &diel_in, const rowvec3 &diel_out, const rowvec2 &interfaces, const double &diel_erf_beta, const mat &charge_position, const rowvec &charge_q, const mat &charge_sigma, const double &grid_multiplier, const bool &trivariate) {
	auto log = spdlog::get("loggers");
	const uword normal_direction = slabcc_cell.normal_direction;
	rowvec Es = zeros<rowvec>(extrapol_steps_num - 1), sizes = Es;
	const mat33 cell_vectors0 = slabcc_cell.vectors;
	const urowvec grid0 = slabcc_cell.grid;
	const rowvec3 extrapolation_grid_size = grid_multiplier * conv_to<rowvec>::from(grid0);
	const urowvec3 extrapolation_grid = { (uword)extrapolation_grid_size(0), (uword)extrapolation_grid_size(1), (uword)extrapolation_grid_size(2) };
	UpdateCell(cell_vectors0, extrapolation_grid);

	for (auto n = 0; n < extrapol_steps_num - 1; ++n) {

		const auto extrapol_factor = extrapol_steps_size * (1.0 + n) + 1;
		UpdateCell(cell_vectors0 * extrapol_factor, extrapolation_grid);
		//extrapolated interfaces
		const rowvec2 interfaces_ext = interfaces / extrapol_factor;

		const mat charge_position_ext = charge_position / extrapol_factor;
		const mat dielectric_profiles = dielectric_profiles_gen(interfaces_ext, diel_in, diel_out, diel_erf_beta);

		cx_cube rhoM(as_size(slabcc_cell.grid), fill::zeros);
		for (uword i = 0; i < charge_position.n_rows; ++i) {
			rhoM += gaussian_charge(charge_q(i), charge_position_ext.row(i).t(), charge_sigma.row(i), trivariate);
		}
		const auto Q = accu(real(rhoM)) * slabcc_cell.voxel_vol;
		// (only works for the orthogonal cells!)
		rhoM -= Q / prod(slabcc_cell.vec_lengths);
		const auto V = poisson_solver_3D(rhoM, dielectric_profiles);
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
	catch (const exception &e) {
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


double Eiso_bessel(double Q, double z0, double sigma, mat diel) {
	auto logger = spdlog::get("loggers");
	const double dK = 0.005; //integration parameter
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
	const rowvec integrand = K % exp(-square(K) * pow(sigma, 2)) % Uk(K, z0, sigma, diel);
	const mat integral = trapz(K, integrand, 1);
 	const double U_total = pow(Q, 2) * integral(0) *  Hartree_to_eV;

	return U_total;

}



rowvec Uk(rowvec k, double z0, double sigma, mat diel) {
	const rowvec3 length = slabcc_cell.vec_lengths;
	const urowvec3 n_points = slabcc_cell.grid;
	const uword normal = slabcc_cell.normal_direction;
	const rowvec Gs = 2.0 * PI / length;

	rowvec Gz0 = ceil(regspace<rowvec>(-0.5 * n_points(normal), 0.5 * n_points(normal) - 1)) * Gs(normal);
	Gz0 = ifftshift(Gz0);
	const rowvec Gz02 = square(Gz0);

	const uword LGz = Gz0.n_elem;
	const double dielbulk = 1;				// bulk response far away in vaccum!
	const cx_mat dielsG = fft(diel - dielbulk);
	const cx_mat dielGz = circ_toeplitz(dielsG.col(normal)) / LGz;
	const uword inplane_direction = (normal == 0) ? 1 : 0;
	const cx_mat dielpGz = circ_toeplitz(dielsG.col(inplane_direction)) / LGz;

	const cx_mat Ag1 = dielGz;
	const cx_mat Ag1p = dielpGz;
	const mat Ag2 = Gz0.t() * Gz0;
	const cx_double coef(0,-1); //compiler-specific problem!
	const cx_mat rhok = exp(coef * Gz0 * z0 - Gz02 * pow(sigma, 2) / 2.0);
	const cx_mat rhok_t = rhok.t();
	const cx_mat rho = ifft(rhok_t);
	rowvec Uk = zeros(arma::size(k));

	const cx_mat Ag12 = Ag1 % Ag2;
	const auto cosGL_2 = cos(Gz0 * length(normal) / 2.0);
	for (int i = 0; i < k.n_elem; ++i) {
		const cx_mat Ag = Ag12 + Ag1p * pow(k(i), 2);
		const double keff = k(i);
		const mat Kinvg = diagmat(dielbulk * length(normal) * (pow(keff, 2) + Gz02) / (1 - exp(-keff * length(normal) / 2.0) * cosGL_2));
		const cx_mat Dg = Kinvg + length(normal) * Ag;
		const auto VGz = solve(Dg, rhok_t);
		const auto Vz = ifft(VGz) * LGz;
		Uk(i) = real(accu(Vz % rho));
	}

	return Uk;
}
