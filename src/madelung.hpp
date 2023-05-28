// Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include "slabcc_math.hpp"

// returns a matrix of the vectors for the repeated images which are inside the
// radius in the lattice
arma::mat generate_shells(const arma::rowvec3 &lattice_vectors,
                          const double &radius);

// reciprocal sum of the madelung
double madelung_reciprocal_sum(const arma::mat &shells,
                               const arma::rowvec3 &reciprocal_lattice_vec,
                               const double &G);
// real sum of the madelung
double madelung_real_sum(const arma::mat &shells,
                         const arma::rowvec3 &lattice_vectors, const double &G);

// returns madelung constant for a point charge in jellium
double jellium_madelung_constant(const arma::mat &shells,
                                 const arma::rowvec3 &lattice_vectors,
                                 const double &G);