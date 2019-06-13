// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "slabcc_input.hpp"

void input_data::verify() const {
	auto log = spdlog::get("loggers");
	charge_sigma = abs(charge_sigma);
	max_eval = abs(max_eval);
	max_time = abs(max_time);
	interfaces = fmod_p(interfaces, 1);
	extrapol_grid_x = abs(extrapol_grid_x);
	opt_grid_x = abs(opt_grid_x);
	opt_tol = abs(opt_tol);
	charge_rotations = fmod_p(charge_rotations + 90, 180) - 90;
	charge_rotations *= PI / 180.0;

	if ((max_eval != 0) && (max_eval < 3)) {
		log->warn("Searching for the optimum model parameters only for {} steps most probably will not be any useful!", max_eval);
	}

	if (diel_in.n_elem == 1) {
		log->trace("Isotropic dielectric constant inside of the slab.");
		diel_in = repelem(diel_in, 1, 3);
	}

	if (diel_out.n_elem == 1) {
		log->trace("Isotropic dielectric constant outside of the slab.");
		diel_out = repelem(diel_out, 1, 3);
	}

	if ((min(diel_in) <= 0) || (min(diel_out) <= 0)) {
		log->debug("Minimum of the dielectric tensor inside the slab: {}", min(diel_in));
		log->debug("Minimum of the dielectric tensor outside the slab: {}", min(diel_out));
		log->critical("The dielectric tensor has not been defined properly! None of the tensor elements should be negative!");
		exit(1);
	}
	if (approx_equal(diel_in, diel_out, "absdiff", 0.02)) { //bulk
		if (optimize_interface) {
			log->trace("The position of the interfaces for the bulk models is irrelevant! \"interfaces\" will not be optimized!");
			optimize_interface = false;
		}
	}
	

	if (optimize_charge_position || optimize_charge_sigma || optimize_charge_fraction || optimize_interface) {
		if ((opt_algo != "BOBYQA") && (opt_algo != "COBYLA") && (opt_algo != "SBPLX")) {
			log->debug("Optimization algorithm: {}", opt_algo);
			log->warn("Unsupported optimization algorithm has been selected!");
			opt_algo = "BOBYQA";
			log->warn("{} will be used instead!", opt_algo);
		}

		if ((optimize_charge_fraction) && (charge_fraction.n_elem == 1)) {
			log->debug("There is only 1 Gaussian charge in your slabcc model. The charge_fraction will not be optimized!");
			optimize_charge_fraction = false;
		}

		if (opt_tol > 1) {
			log->debug("Requested optimization tolerance: {}", opt_tol);
			log->warn("The relative optimization tolerance is unacceptable! It must be in choosen in (0-1) range.");
			opt_tol = 0.01;
			log->warn("optimize_tolerance = {} will be used!", opt_tol);
		}
	}


	if (charge_position.n_cols != 3) {
		log->debug("Number of the parameters defined for the position of a charge: {}", charge_position.n_cols);
		log->critical("Incorrect definition of charge positions!");
		log->critical("Positions should be defined as: charge_position = 0.1 0.2 0.3; 0.1 0.2 0.4;");
		finalize_loggers();
		exit(1);
	}

	const uword charge_number = charge_position.n_rows;
	const uword sigma_rows = charge_sigma.n_rows;
	const uword sigma_cols = charge_sigma.n_cols;

	if (charge_sigma.n_elem == 1) {
		charge_sigma = charge_sigma(0) * ones(charge_number, 3);
		log->debug("Only one charge_sigma is defined!");

	}
	else if (sigma_rows == charge_number) {
		if (sigma_cols == 3) {
			if (trivariate) {
				log->debug("All the charge_sigma values are properly defined!");
			}
			else {	
				const mat isotropic_sig = repmat(charge_sigma.col(0), 1, 3);
				const mat sig_diff = abs(isotropic_sig - charge_sigma);
				if (any(vectorise(sig_diff) > 0.01)) {
					charge_sigma = isotropic_sig;
					log->warn("charge_sigma is not defined properly! charge_sigma={} will be used.", to_string(charge_sigma.col(0)));
				}		
			}
		}
		else if (sigma_cols == 1) {
			charge_sigma = repmat(charge_sigma, 1, 3);
			if (trivariate) {
				log->debug("Equal values will be assumed for the charge_sigma in all directions!");
			}
			else {
				log->debug("charge_sigma is properly defined!");
			}
		}
	}
	else if (sigma_cols == charge_number) {
		if (sigma_rows == 1) {
			charge_sigma = repmat(charge_sigma.t(), 1, 3);
			if (trivariate) {
				log->debug("Equal values will be assumed for the charge_sigma in all directions!");
			}
			else {
				log->debug("charge_sigma is properly defined!");
			}
		}
	}
	else if ((sigma_cols == 3) && (sigma_rows == 1)) {
		if (trivariate) {
			charge_sigma = repmat(charge_sigma, charge_number, 1);
			log->debug("Equal charge_sigma will be assumed for all the Gaussian charges!");
		}
	}

	if (charge_sigma.min() < 0.001) {
		log->warn("charge_sigma is too small. For accurate energy calculations with this parameter a very large grid is required!");
	}

	//if all the sensible checks and remedies have failed!
	if (arma::size(charge_sigma) != arma::size(charge_position)) {
		log->warn("Number of the defined Gaussian charges and the charge_sigma sets does not match!");
		charge_sigma = mat(arma::size(charge_position), fill::ones);
		if (trivariate) {
			log->warn("charge_sigma = {} will be used!", to_string(charge_sigma));
		}
		else {
			log->warn("charge_sigma = {} will be used!", to_string(charge_sigma.col(0)));
		}
	}

	log->debug("charge_sigma after the checks: {}", to_string(charge_sigma));

	if ((abs(charge_rotations)).max() > 0) {
		if (!trivariate) {
			log->warn("charge_rotation will be ignored for the simple Gaussian charges!");
			charge_rotations = zeros<mat>(arma::size(charge_position));
		}
		else {
			log->debug("charge_rotation (rad): {}", to_string(charge_rotations));
		}
	}

	if ((charge_rotations.n_rows == 1) && (charge_number > 1)) {
		charge_rotations = repmat(charge_rotations, charge_number, 1);
		log->debug("Equal charge_rotation will be assumed for all the Gaussian charges!");
	}

	if (arma::size(charge_rotations) != arma::size(charge_position)) {
		log->warn("charge_rotation is not defined properly!");
		charge_rotations = zeros<mat>(arma::size(charge_position));
		log->warn("charge_rotation = {} will be used!", to_string(charge_rotations));
	}


	//charge_fraction
	if (charge_fraction.n_elem != charge_position.n_rows) {
		log->debug("Number of charge position values: {}", charge_position.n_rows);
		log->debug("Number of charge fraction values: {}", charge_fraction.n_elem);
		charge_fraction = rowvec(charge_position.n_rows, fill::ones);
		log->warn("Number of the charge_fraction and charge_position sets does not match!");
		log->warn("Equal charge fractions will be assumed!");

	}


	if (!file_exists(CHGCAR_neutral) || !file_exists(CHGCAR_charged)
		|| !file_exists(LOCPOT_neutral) || !file_exists(LOCPOT_charged)) {
		log->debug("CHGCAR_neutral: '{}' found: {}", CHGCAR_neutral, to_string(file_exists(CHGCAR_neutral)));
		log->debug("CHGCAR_charged: '{}' found: {}", CHGCAR_charged, to_string(file_exists(CHGCAR_charged)));
		log->debug("LOCPOT_neutral: '{}' found: {}", LOCPOT_neutral, to_string(file_exists(LOCPOT_neutral)));
		log->debug("LOCPOT_charged: '{}' found: {}", LOCPOT_charged, to_string(file_exists(LOCPOT_charged)));
		log->critical("One or more of the input files could not be found!");
		finalize_loggers();
		exit(1);
	}

	if (extrapolate) {
		if (extrapol_steps_num < 3) {
			log->debug("Requested extrapolation steps: {}", extrapol_steps_num);
			log->warn("At least 3 steps are needed for verifiability of the extrapolation algorithm!");
			extrapol_steps_num = 3;
			log->warn("{} steps will be used instead!", extrapol_steps_num);
		}
	}
	else {
		if (model_2D) {
			// Bessel expansion limitations
			if (charge_fraction.n_elem > 1) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support multiple Gaussian charges. The extrapolation method \"extrapolate=yes\" will be used for this calculation.");
				extrapolate = true;
			}
			if (any(diel_out > 1)) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support the models embedded in any dielectric medium other than the vacuum. The extrapolation method \"extrapolate=yes\" will be used for this calculation.");
				extrapolate = true;
			}
			if (trivariate) {
				log->error("The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support the trivariate Gaussian model charges. The extrapolation method \"extrapolate=yes\" will be used for this calculation.");
				extrapolate = true;
			}

			vec2 inplane_diels;
			for (uword i = 0; i < 3; ++i) {
				if (i < normal_direction) {
					inplane_diels(i) = diel_in(i);
				}
				if (i > normal_direction) {
					inplane_diels(i - 1) = diel_in(i);
				}
			}

			if (abs(inplane_diels(0) - inplane_diels(1)) > 0.01) {
				log->debug("In-plane dielectric constants: {}", ::to_string(inplane_diels));
				log->error("In-plane dielectric constants are not equal. The current implementation of the E_isolated calculation algorithm using the Bessel expansion of the Poisson eq. does not support this model. Will use the extrapolation method \"extrapolate=yes\" for this calculation.");
				extrapolate = true;
			}
		}
	}

	log->trace("Input parameters verified!");

}

void input_data::parse(const string& input_file) const {
	auto log = spdlog::get("loggers");
	INIReader reader(input_file);
	if (reader.ParseError() < 0) {
		log->critical("Cannot load the input file: {}", input_file);
		finalize_loggers();
		exit(1);
	}

	verbosity_level = reader.GetInteger("verbosity", 1);

	CHGCAR_neutral = reader.GetStr("CHGCAR_neutral", "CHGCAR.N");
	LOCPOT_neutral = reader.GetStr("LOCPOT_neutral", "LOCPOT.N");
	LOCPOT_charged = reader.GetStr("LOCPOT_charged", "LOCPOT.C");
	CHGCAR_charged = reader.GetStr("CHGCAR_charged", "CHGCAR.C");
	charge_position = reader.GetMat("charge_position", {});
	charge_fraction = reader.GetVec("charge_fraction", rowvec(charge_position.n_rows, fill::ones) / charge_position.n_rows);
	trivariate = reader.GetBoolean("charge_trivariate", false);
	charge_rotations = reader.GetMat("charge_rotation", zeros<mat>(arma::size(charge_position)));
	charge_sigma = reader.GetMat("charge_sigma", ones<mat>(arma::size(charge_position)));
	slabcenter = reader.GetVec("slab_center", { 0.5, 0.5, 0.5 });
	normal_direction = xyz2int(reader.GetStr("normal_direction", "z"));
	interfaces = reader.GetVec("interfaces", { 0.25, 0.75 });
	diel_in = reader.GetVec("diel_in", { 1,1,1 });
	diel_out = reader.GetVec("diel_out", { 1,1,1 });
	diel_erf_beta = reader.GetReal("diel_taper", 1);
	optimize_charge_position = reader.GetBoolean("optimize_charge_position", true);
	optimize_charge_sigma = reader.GetBoolean("optimize_charge_sigma", true);
	optimize_charge_rotation = reader.GetBoolean("optimize_charge_rotation", false);
	optimize_charge_fraction = reader.GetBoolean("optimize_charge_fraction", true);
	optimize_interface = reader.GetBoolean("optimize_interfaces", true);
	model_2D = reader.GetBoolean("2d_model", false);
	opt_algo = reader.GetStr("optimize_algorithm", "BOBYQA");
	opt_tol = reader.GetReal("optimize_tolerance", 0.01);
	max_eval = reader.GetInteger("optimize_maxsteps", 0);
	max_time = reader.GetInteger("optimize_maxtime", 0);
	opt_grid_x = reader.GetReal("optimize_grid_x", 0.8);
	extrapolate = reader.GetBoolean("extrapolate", model_2D ? false : true);
	extrapol_grid_x = reader.GetReal("extrapolate_grid_x", 1);
	extrapol_steps_num = reader.GetInteger("extrapolate_steps_number", 4);
	extrapol_steps_size = reader.GetReal("extrapolate_steps_size", 0.5);

	reader.dump_parsed();

}