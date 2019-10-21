#pragma once

#include <armadillo>

const double ang_to_bohr = 1e-10 / arma::datum::a_0; // 1.88972612546
const double Hartree_to_eV = arma::datum::R_inf * arma::datum::h *
                             arma::datum::c_0 / arma::datum::eV *
                             2; // 27.2113860193
const double PI = arma::datum::pi;