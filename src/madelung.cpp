// Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "madelung.hpp"

arma::mat generate_shells(const arma::rowvec3 &lattice_vectors,
                          const double &radius) {
  // we will always have symmetry in our orthogonal cell (orthorombic)
  // so we check the distance of the 1st octant and use mirror planes to
  // generate full list

  // vectors with positive coefficients inside a box with smaller side than
  // radius
  arma::umat positive_box_vectors;

  // index of the box to search for distances
  arma::rowvec3 max_indexes = arma::floor(radius / lattice_vectors);
  // max_indexes.print("max_index");

  for (arma::uword i = 0; i <= max_indexes(0); ++i) {
    for (arma::uword j = 0; j <= max_indexes(1); ++j) {
      for (arma::uword k = 0; k <= max_indexes(2); ++k) {
        positive_box_vectors.insert_rows(positive_box_vectors.n_rows,
                                         arma::urowvec{i, j, k});
      }
    }
  }

  // n=0 term is omitted
  positive_box_vectors.shed_row(0);
  // positive_box_vectors.print("positive_box_vectors");

  // vectors which end up inside the sphere with requested radius
  arma::umat positive_inside_vectors;
  for (arma::uword i = 0; i < positive_box_vectors.n_rows; ++i) {
    if (arma::norm(positive_box_vectors.row(i) % lattice_vectors) < radius) {
      positive_inside_vectors.insert_rows(positive_inside_vectors.n_rows,
                                          positive_box_vectors.row(i));
    }
  }

  // vectors in all octants with smaller distance than radius
  arma::mat all_inside_vectors =
      arma::conv_to<arma::mat>::from(positive_inside_vectors);
  for (arma::uword i = 0; i < 3; ++i) {
    for (arma::uword j = 0; j < all_inside_vectors.n_rows; ++j) {
      if (all_inside_vectors(j, i) > 0) {
        all_inside_vectors.insert_rows(all_inside_vectors.n_rows,
                                       all_inside_vectors.row(j));
        all_inside_vectors(all_inside_vectors.n_rows - 1, i) *= -1;
      }
    }
  }

  return all_inside_vectors;
}

double madelung_real_sum(const arma::mat &shells,
                         const arma::rowvec3 &lattice_vectors,
                         const double &G) {
  double sum = 0;
  for (arma::uword i = 0; i < shells.n_rows; ++i) {
    const double distance = arma::norm(shells.row(i) % lattice_vectors);
    sum += std::erfc(G * distance) / distance;
  }
  return sum;
}

double madelung_reciprocal_sum(const arma::mat &shells,
                               const arma::rowvec3 &reciprocal_lattice_vec,
                               const double &G) {
  double sum = 0;
  for (arma::uword i = 0; i < shells.n_rows; ++i) {
    const double G2 =
        arma::accu(arma::square(shells.row(i) % reciprocal_lattice_vec));
    sum += std::exp(-G2 / (4 * std::pow(G, 2))) / (G2 / (4 * std::pow(G, 2)));
  }
  return sum / std::pow(G, 2);
}

double jellium_madelung_constant(const arma::mat &shells,
                                 const arma::rowvec3 &lattice_vectors,
                                 const double &G) {
  const double real_sum = madelung_real_sum(shells, lattice_vectors, G);
  const double reciprocal_sum =
      madelung_reciprocal_sum(shells, 2 * PI / lattice_vectors, G);
  const double inverse_unit_vol = PI / arma::prod(lattice_vectors);
  double madelung_constant =
      -(inverse_unit_vol * reciprocal_sum + real_sum - 2 * G / std::sqrt(PI) -
        inverse_unit_vol / std::pow(G, 2));
  madelung_constant *= std::pow(arma::prod(lattice_vectors), 1.0 / 3);
  return madelung_constant;
}
