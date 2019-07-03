// Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include "stdafx.h"
#include "clara.hpp"
#include "spdlog.h"
#include "arma_io.hpp"
#include "nlopt.hpp"
#include "sinks/basic_file_sink.h"
#include "sinks/stdout_color_sinks.h"

using namespace std;

extern int verbosity_level;

struct cli_params {
	string &input_file, &output_file, &log_file;
	bool &diff_only;

	// reads the command line and sets the input_file and output_file
	void parse(int argc, char *argv[]);

};

//defines the minimum verbosity level for each type of action
enum class verbosity:int {
	info = 1,					//spdlog->info()
	write_normal_planarAvg = 1, //write the planar average of defect and model LOCPOT files
	debug = 2,					//spdlog->debug()
	write_defect_file = 2,		//write extra charge density, extra charge potential
	write_dielectric_file = 2,	//write the dielectric profile of model
	write_planarAvg_file = 3,	//write the planar average of the CHGCAR-LOCPOT files for neutral/charged/defect/model 
	trace = 4,					//spdlog->trace()
	timing = 4,					//log the time passed from the start of the program
};

//converts the first letter of the string from "a-b-c"/"x-y-z" to 0/1/2
unsigned int xyz2int(const string& s);

//converts 0/1/2 to X/Y/Z
char int2xyz(const unsigned int& i);

inline bool file_exists(const string& name) {
	return ifstream(name.c_str()).good();
}

//decides if a functionality is active with the current verbosity level
inline bool is_active(const verbosity& action) noexcept {
	return verbosity_level >= static_cast<int>(verbosity(action));
}

inline bool file_is_empty(ifstream& file) {
	return file.peek() == ifstream::traits_type::eof();
}

string tolower(string in_str) noexcept;

void initialize_loggers(const string& log_file, const string& output_file);
void update_loggers();
void finalize_loggers();

//renames the old output file if possible or changes the output file name
void prepare_output_file(string& output_file);

// bool to yes/no conversion
string to_string(const bool& b);



