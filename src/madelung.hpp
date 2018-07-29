#pragma once
#include "stdafx.h"
#include "slabcc_math.hpp"

using namespace std;

//returns a matrix of the vectors for the repeated images which are inside the radius in the lattice
mat generate_shells(const rowvec3 &lattice_vectors, const double &radius);
//reciprocal sum of the madelung
double madelung_reciprocal_sum(const mat &shells, const rowvec3 &reciprocal_lattice_vec, const double &G);
//real sum of the madelung
double madelung_real_sum(const mat &shells, const rowvec3 &lattice_vectors, const double &G);
//returns madelung constant for a point charge in jellium
double jellium_madelung_constant(const mat &shells, const rowvec3 &lattice_vectors, const double &G);