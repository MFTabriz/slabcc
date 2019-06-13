// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "arma_io.hpp"

void write_mat2file(const mat& input, const string& output_file) {
	ofstream out_file;
	out_file.open(output_file);
	out_file << fixed << showpos << setprecision(15);
	input.each_row([&out_file](const rowvec & row) {
		row.for_each([&out_file](const double& val) { out_file << val << " "; });
		out_file << '\n';
		});
	out_file.close();
}

