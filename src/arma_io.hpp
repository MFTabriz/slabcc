// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include <armadillo>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace arma;

//single-line output for vec
template<typename T>
ostream &operator << (ostream &o, const Row<T> &vec) {
	for (uword elem = 0; elem < vec.n_elem; ++elem) {
		o << vec(elem);
		if (elem != vec.n_elem - 1) {
			o << " ";
		}
	}
	return o;
}

template<typename T>
ostream &operator << (ostream &o, const subview<T> &vec) {
	for (uword elem = 0; elem < vec.n_elem; ++elem) {
		o << vec(elem);
		if (elem != vec.n_elem - 1) {
			o  << " ";
		}
	}

	return o;
}

//single-line output for mat
template<typename T>
ostream &operator << (ostream &o, const Mat<T> &mat) {
	mat.each_row([&o](const Row<T> &row) { o << row << "; "; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, Row<T> &vec) {
	vec.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, subview_row<T> vec) {
	vec.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, Cube<T> &c) {
	c.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template<typename T>
istream &operator >> (istream &o, Mat<T> &m) {
	m.for_each([&o](T &elem) { o >> elem; });
	return o;
}

template <typename T>
string to_string(T input) {
	ostringstream output;
	output << setprecision(15);
	output << input;
	return output.str();
}


