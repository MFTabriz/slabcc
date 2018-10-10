// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include <armadillo>

using namespace std;
using namespace arma;

//single-line output for vec
template<typename T>
ostream &operator << (ostream &o, const Row<T> &vec) {
	stringstream ss;
	ss << setprecision(12);
	vec.for_each([&ss](const T &elem) { ss << elem << " "; });
	string st = ss.str();
	st.pop_back();
	o << st;
	return o;
}

template<typename T>
ostream &operator << (ostream &o, const subview<T> &vec) {
	stringstream ss;
	ss << setprecision(12);
	vec.for_each([&ss](const T &elem) { ss << elem << " "; });
	string st = ss.str();
	st.pop_back();
	o << st;
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
	output << setprecision(12);
	output << input;
	return output.str();
}


