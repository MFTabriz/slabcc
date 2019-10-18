// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include <armadillo>
#include <iomanip>
#include <iostream>

// single-line output for vec
template <typename T>
std::ostream &operator<<(std::ostream &o, const arma::Row<T> &vec) {
  for (arma::uword elem = 0; elem < vec.n_elem; ++elem) {
    o << vec(elem);
    if (elem != vec.n_elem - 1) {
      o << " ";
    }
  }
  return o;
}

template <typename T>
std::ostream &operator<<(std::ostream &o, const arma::subview<T> &vec) {
  for (arma::uword elem = 0; elem < vec.n_elem; ++elem) {
    o << vec(elem);
    if (elem != vec.n_elem - 1) {
      o << " ";
    }
  }

  return o;
}

// single-line output for mat
template <typename T>
std::ostream &operator<<(std::ostream &o, const arma::Mat<T> &mat) {
  mat.each_row([&o](const arma::Row<T> &row) { o << row << "; "; });
  return o;
}

template <typename T>
std::istream &operator>>(std::istream &o, arma::Row<T> &vec) {
  vec.for_each([&o](T &elem) { o >> elem; });
  return o;
}

template <typename T>
std::istream &operator>>(std::istream &o, arma::subview_row<T> vec) {
  vec.for_each([&o](T &elem) { o >> elem; });
  return o;
}

template <typename T>
std::istream &operator>>(std::istream &o, arma::Cube<T> &c) {
  c.for_each([&o](T &elem) { o >> elem; });
  return o;
}

template <typename T>
std::istream &operator>>(std::istream &o, arma::Mat<T> &m) {
  m.for_each([&o](T &elem) { o >> elem; });
  return o;
}

template <typename T> std::string to_string(T input) {
  std::ostringstream output;
  output << std::setprecision(15) << input;
  return output.str();
}
