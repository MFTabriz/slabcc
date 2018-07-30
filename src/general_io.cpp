// Copyright (c) 2018, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "general_io.hpp"


void write_vec2file(const vector<double>& input, const string& output_file) {
	ofstream out_file;
	out_file.open(output_file);
	out_file << fixed << showpos << setprecision(10);
	if (is_active(verbosity::more_digits)) {
		out_file << setprecision(15);
	}
	for (const auto& i : input) {
		out_file << i << endl;
	}
	out_file.close();

}

vector<double> read_file2vec(const string& input_file) {
	vector<double> vector;
	string temp;
	ifstream in_file(input_file);

	while (getline(in_file, temp)) {
		vector.push_back(stod(temp));
	}
	in_file.close();

	return vector;
}

unsigned int xyz2int(const string& s) {
	const char dir_char = tolower(s.at(0));
	const unordered_map<char, unsigned int> dir_map{
		{ '0', 0 }, { 'x', 0 }, { 'a', 0 },
		{ '1', 1 }, { 'y', 1 }, { 'b', 1 },
		{ '2', 2 }, { 'z', 2 }, { 'c', 2 },
	};
	return dir_map.at(dir_char);
}

string timing() {
	if (is_active(verbosity::timing)) {
		stringstream streamObj;
		auto time_diff_h = chrono::duration_cast<chrono::hours>(chrono::steady_clock::now() - t0);
		auto time_diff_m = chrono::duration_cast<chrono::minutes>(chrono::steady_clock::now() - t0);
		auto time_diff_s = chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now() - t0);
		streamObj << fixed << setfill('0');
		streamObj << setw(2) << time_diff_h.count() << ":" << setw(2) << time_diff_m.count() % 60 << ":" << setw(2) << time_diff_s.count() % 60 << " ";
		return streamObj.str();
	}
	return "";
}

string tolower(string in_str) noexcept {
	for (char& c : in_str) {
		c = tolower(c);
	}
	return in_str;
}

void parse_cli(int argc, char *argv[], string& input_file, string &output_file) {

	bool showHelp = false, showManual = false, showVer = false, showAttr = false;
	auto cli = clara::Help(showHelp) |
		clara::Opt(input_file, "input_file")
		["-i"]["--input"]
		("slabcc input file name") |
		clara::Opt(output_file, "input_file")
		["-o"]["--output"]
		("slabcc output file name") |
		clara::Opt(showManual)
		["-m"]["--man"]
		("show quick start guide") |
		clara::Opt(showVer)
		["-v"]["--version"]
		("show version and compilation date") |
		clara::Opt(showAttr)
		["-c"]["--copyright"]
		("show copyright information and the attributions");

	auto cli_result = cli.parse(clara::Args(argc, argv));
	if (!cli_result) {
		cerr << "Error in command line: " << cli_result.errorMessage() << endl;
		exit(1);
	}

	if (showHelp) {
		ostringstream oss;
		oss << (clara::Parser() | cli);
		cout << oss.str();
		exit(0);
	}

	if (showVer) {
		cout << "SLAB Charge Correction (slabcc)\n"
			"Version: " << version << endl;
		cout << "Compiled: " << __DATE__ << " " << __TIME__ << endl;

		exit(0);
	}

	if (showManual) {
		cout << "Quick user guide:\n"
			"To calculate the charge correction slabcc needs the following files:\n"
			" 1. Input parameters file (default: slabcc.in)\n"
			" 2. CHGCAR of the neutral system\n"
			" 3. CHGCAR of the charged system\n"
			" 4. LOCPOT of the neutral system\n"
			" 5. LOCPOT of the charged system\n\n"
			"The input parameters file for a slab should minimally include:\n"
			" - charge_position : approximate position of the localized charge\n"
			" - diel_in : dielectric constant of the slab\n"
			" - normal_direction : direction normal to the surface\n"
			" - interfaces : surface positions of the slab\n";
		exit(0);
	}

	if (showAttr) {
		cout << "Copyright (c) 2018, Bremen Center for Computational Materials Science (BCCMS), M. Farzalipour Tabriz\n"
			"The source code and all the documentations are available under The 2-Clause BSD License. For more information see LICENSE.TXT.\n\n"
			"Included libraries\n\n"
			"-Armadillo C++ Linear Algebra Library: licensed under the Apache License 2.0\n"
			"  Copyright 2008 - 2018 Conrad Sanderson\n"
			"  Copyright 2008 - 2016 National ICT Australia (NICTA)\n"
			"  Copyright 2017 - 2018 Arroyo Consortium\n"
			"  Copyright 2017 - 2018 Data61, CSIRO\n"
			"  This product includes software developed by Conrad Sanderson\n"
			"  This product includes software developed at National ICT Australia (NICTA)\n"
			"  This product includes software developed at Arroyo Consortium\n"
			"  This product includes software developed at Data61, CSIRO\n\n"
			"-inih(INI Not Invented Here): licensed under the 3-clause BSD license\n"
			"  (c) 2009, Ben Hoyt, et al.\n\n"
			"-clara: licensed under the Boost Software License 1.0\n"
			"  (c) 2017-2018 Phil Nash, Martin Horenovsky, et al.\n\n"
			"-Spline (Cubic Spline Interpolation): licensed under the Beer-Ware License 42\n"
			"  (c) 2011, Devin Lane\n\n"
			"-NLOPT: licensed under GNU Lesser General Public License (LGPL)\n"
			"  (c) 2007-2010, Massachusetts Institute of Technology\n";
		exit(0);
	}
}
