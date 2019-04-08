// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "general_io.hpp"


void write_vec2file(const vector<double>& input, const string& output_file) {
	ofstream out_file;
	auto log = spdlog::get("loggers");
	log->trace("Started writing " + output_file);

	out_file.open(output_file);
	out_file << fixed << showpos << setprecision(15);
	for (const auto& i : input) { out_file << i << '\n'; }
	out_file.close();

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

char int2xyz(const unsigned int& i) {
	const unordered_map<unsigned int, char> dir_map{
		{ 0 , 'X'}, { 1 , 'Y' }, { 2 , 'Z' },
	};
	return dir_map.at(i);
}


string tolower(string in_str) noexcept {
	for (char& c : in_str) {
		c = tolower(c);
	}
	return in_str;
}

void parse_cli(int argc, char *argv[], cli_params& cli_parameters) {

	bool showHelp = false, showManual = false, showVer = false, showAttr = false;
	auto cli = clara::Help(showHelp) |
		clara::Opt(cli_parameters.input_file, "input_file")
		["-i"]["--input"]
		("slabcc input file name") |
		clara::Opt(cli_parameters.output_file, "input_file")
		["-o"]["--output"]
		("slabcc output file name") |
		clara::Opt(cli_parameters.log_file, "log_file")
		["-l"]["--log"]
		("slabcc log file name") |
		clara::Opt(cli_parameters.diff_only)
		["-d"]["--diff"]
		("calculate the charge and the potential differences only") |
		clara::Opt(showManual)
		["-m"]["--man"]
		("show the quick start guide") |
		clara::Opt(showVer)
		["-v"]["--version"]
		("show the slabcc version and its compilation date") |
		clara::Opt(showAttr)
		["-c"]["--copyright"]
		("show the copyright information and the attributions");

	auto cli_result = cli.parse(clara::Args(argc, argv));
	if (!cli_result) {
		cerr << "Error in command line: " << cli_result.errorMessage() << '\n';
		exit(1);
	}

	if (showHelp) {
		ostringstream oss;
		oss << (clara::Parser() | cli);
		cout << oss.str();
		exit(0);
	}

	if (showVer) {
		cout << "SLAB Charge Correction (slabcc)" << '\n';
		cout << "Version: " << SLABCC_VERSION_MAJOR << "." << SLABCC_VERSION_MINOR << "." << SLABCC_VERSION_PATCH << '\n';
		cout << "Compilation: " << __DATE__ << " " << __TIME__ << '\n';
		exit(0);
	}

	if (showManual) {
		cout << "Quick start guide:\n"
			"To calculate the charge correction, slabcc needs the following files:\n"
			" 1. Input parameters file (default: slabcc.in)\n"
			" 2. CHGCAR of the neutral system\n"
			" 3. CHGCAR of the charged system\n"
			" 4. LOCPOT of the neutral system\n"
			" 5. LOCPOT of the charged system\n\n"
			"The input parameters file for a slab should minimally include:\n"
			" - charge_position : approximate position of the localized charge center\n"
			" - diel_in : dielectric constant/tensor of the slab\n"
			" - normal_direction : direction normal to the surface\n"
			" - interfaces : surface positions of the slab\n";
		exit(0);
	}

	if (showAttr) {
		cout << "Copyright (c) 2018-2019, Bremen Center for Computational Materials Science (BCCMS), M. Farzalipour Tabriz\n"
			"The source code and all the documentations are available under The 2-Clause BSD License. For more information see LICENSE.TXT.\n\n"
			"Included libraries\n\n"
			"-Armadillo C++ Linear Algebra Library: licensed under the Apache License 2.0\n"
			"  Copyright 2008-2018, Conrad Sanderson\n"
			"  Copyright 2008-2016, National ICT Australia (NICTA)\n"
			"  Copyright 2017-2018, Arroyo Consortium\n"
			"  Copyright 2017-2018, Data61, CSIRO\n"
			"  This product includes software developed by Conrad Sanderson\n"
			"  This product includes software developed at National ICT Australia (NICTA)\n"
			"  This product includes software developed at Arroyo Consortium\n"
			"  This product includes software developed at Data61, CSIRO\n\n"
			"-inih(INI Not Invented Here): licensed under the 3-clause BSD license\n"
			"  (c) 2009, Ben Hoyt, et al.\n\n"
			"-clara: licensed under the Boost Software License 1.0\n"
			"  (c) 2017-2018, Phil Nash, Martin Horenovsky, et al.\n\n"
			"-Spline (Cubic Spline Interpolation): licensed under the Beer-Ware License 42\n"
			"  (c) 2011, Devin Lane\n\n"
			"-NLOPT: licensed under GNU Lesser General Public License (LGPL)\n"
			"  (c) 2007-2010, Massachusetts Institute of Technology\n\n"
			"-spdlog: licensed under the MIT License\n"
		    "  (c) 2016, Gabi Melman\n\n"
			"-Boost.Predef: licensed under the Boost Software License 1.0\n"
			"  (c) 2005-2018 Rene Rivera\n"
			"  (c) 2015 Charly Chevalier\n"
			"  (c) 2015 Joel Falcou\n";

		exit(0);
	}
}


void initialize_logger(string log_file) {

	vector<spdlog::sink_ptr> sinks;
	sinks.push_back(make_shared<spdlog::sinks::stdout_color_sink_mt>());
	sinks.push_back(make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true));
	sinks.push_back(make_shared<spdlog::sinks::basic_file_sink_mt>("MSG", true));
	auto combined_logger = make_shared<spdlog::logger>("loggers", begin(sinks), end(sinks));
	sinks.at(2)->set_level(spdlog::level::warn);
	sinks.at(2)->set_pattern("[%l] %v");
	combined_logger->set_pattern("%v");
	combined_logger->set_level(spdlog::level::info);
	spdlog::register_logger(combined_logger);
}

void logger_update() {

	auto log = spdlog::get("loggers");
	string log_pattern = "%v";

	if (is_active(verbosity::info)) {
		log->set_level(spdlog::level::info);
	}
	if (is_active(verbosity::debug)) {
		log->set_level(spdlog::level::debug);
		log_pattern = "[%^%l%$] " + log_pattern;
	}
	if (is_active(verbosity::trace)){
		log->set_level(spdlog::level::trace);
	}

	if (is_active(verbosity::timing)) {
		log_pattern = "[%K]" + log_pattern;
	}

	log->set_pattern(log_pattern);
	log->sinks().at(2)->set_pattern("[%^%l%$] %v");
}

string to_string(const bool& b) {
	return b ? "yes" : "no";
}
