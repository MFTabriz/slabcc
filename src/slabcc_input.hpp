#pragma once
#include "slabcc_math.hpp"
#include "slabcc_cell.hpp"
#include "INIReader.h"

//references to the input data variables
struct input_data {
	string &CHGCAR_neutral, &LOCPOT_charged, &LOCPOT_neutral, &CHGCAR_charged, &opt_algo;
	mat &charge_position;
	rowvec &charge_fraction;
	mat &charge_sigma, &charge_rotations;
	rowvec3 &slabcenter;
	rowvec &diel_in, &diel_out;
	uword &normal_direction;
	rowvec2 &interfaces;
	double &diel_erf_beta, &opt_tol;
	bool &optimize_charge_position, &optimize_charge_sigma, &optimize_charge_rotation, &optimize_charge_fraction, &optimize_interface, &extrapolate, &model_2D, &trivariate;
	double &opt_grid_x, &extrapol_grid_x;
	int &max_eval, &max_time, &extrapol_steps_num;
	double &extrapol_steps_size;

	//read the input variables from the input_file
	void parse(const string& input_file) const;

	//sanity checks on the input parameters
	void verify() const;
};