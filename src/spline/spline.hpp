#pragma once

/* "THE BEER-WARE LICENSE" (Revision 42): Devin Lane wrote this file. As long as you retain
* this notice you can do whatever you want with this stuff. If we meet some day, and you
* think this stuff is worth it, you can buy me a beer in return.
https://shiftedbits.org/2011/01/30/cubic-spline-interpolation/

Code has been converted to Armadillo data type to use in slabcc by https://github.com/MFTabriz

*/
#include <vector>
#include <iostream>
#include <armadillo>


template <typename X>
class Spline {
public:

	// A spline with x and y values
	Spline(const arma::Row<X>& x, const arma::Col<X>& y) {
		if (x.n_elem != y.n_elem) {
			std::cerr << "X and Y must be the same size " << std::endl;
			return;
		}

		if (x.n_elem < 3) {
			std::cerr << "Must have at least three points for interpolation" << std::endl;
			return;
		}

		const arma::uword n = y.n_elem - 1;

		arma::Row<X> b(n), d(n), a(n), c(n + 1), l(n + 1), u(n + 1), z(n + 1), h(n + 1);

		l(0) = 1;
		u(0) = 0; z(0) = 0;
		h(0) = x(1) - x(0);
		for (arma::uword i = 1; i < n; ++i) {
			h(i) = x(i + 1) - x(i);
			l(i) = 2 * (x(i + 1) - x(i - 1)) - h(i - 1) * u(i - 1);
			u(i) = h(i) / l(i);
			a(i) = (3 / h(i)) * (y(i + 1) - y(i)) - (3 / h(i - 1)) * (y(i) - y(i - 1));
			z(i) = (a(i) - h(i - 1) * z(i - 1)) / l(i);
		}
		l(n) = 1;
		z(n) = 0; c(n) = 0;
		for (arma::sword j = n - 1; j >= 0; --j) {
			c(j) = z(j) - u(j) * c(j + 1);
			b(j) = (y(j + 1) - y(j)) / h(j) - (h(j) * (c(j + 1) + 2 * c(j))) / 3;
			d(j) = (c(j + 1) - c(j)) / (3 * h(j));
		}
		for (arma::uword i = 0; i < n; ++i) {
			mElements.push_back(Element(x(i), y(i), b(i), c(i), d(i)));
		}
	}

	//return the value of the spline function for x
	X interpolate(const X&x) const {
		if (mElements.empty()) return X();

		auto it = std::lower_bound(mElements.begin(), mElements.end(), element_type(x));
		if (it != mElements.begin()) { --it; }
		return it->eval(x);
	}

	//return the value of the spline function for a Row
	arma::Col<X> interpolate(const arma::Row<X>& xx) const {
		if (mElements.empty()) return arma::Row<X>(xx.size());

		arma::Col<X> ys = arma::zeros<arma::Col<X>>(xx.n_elem);

		for (arma::uword it = 0; it < xx.n_elem; ++it) {
			ys(it) = interpolate(xx(it));
		}
		return ys;
	}

protected:

	class Element {
	public:

		X x = 0, a = 0, b = 0, c = 0, d = 0;

		Element(X _x) noexcept : x(_x) {}
		Element(X _x, X _a, X _b, X _c, X _d) noexcept
			: x(_x), a(_a), b(_b), c(_c), d(_d) {}

		X eval(const X& xx) const noexcept {
			X xix(xx - x);
			return a + b * xix + c * (xix * xix) + d * (xix * xix * xix);
		}

		bool operator<(const Element& e) const noexcept {
			return x < e.x;
		}
	};
	typedef Element element_type;
	std::vector<element_type> mElements;

};
