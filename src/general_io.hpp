// Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#pragma once
#include "arma_io.hpp"
#include "clara.hpp"
#include "nlopt.hpp"
#include "sinks/basic_file_sink.h"
#include "sinks/stdout_color_sinks.h"
#include "slabcc_info.hpp"
#include "spdlog.h"

extern int verbosity_level;

struct cli_params {
  std::string &input_file, &output_file, &log_file;
  bool &diff_only;

  // reads the command line and sets the input_file and output_file
  void parse(int argc, char *argv[]);
};

// defines the minimum verbosity level for each type of action
enum class verbosity : int {
  info = 1, // spdlog->info()
  write_normal_planarAvg = 1,
  debug = 2, // spdlog->debug()
  write_defect_file = 2,
  write_dielectric_file = 2,
  write_planarAvg_file = 3,
  trace = 4, // spdlog->trace()
  timing = 4,
};

// converts the first letter of the string from "a-b-c"/"x-y-z" to 0/1/2
unsigned int xyz2int(const std::string &s);

// converts 0/1/2 to X/Y/Z
char int2xyz(const unsigned int &i);

inline bool file_exists(const std::string &name) {
  return std::ifstream(name.c_str()).good();
}

// decides if a functionality is active with the current verbosity level
inline bool is_active(const verbosity &action) noexcept {
  return verbosity_level >= static_cast<int>(verbosity(action));
}

inline bool file_is_empty(std::ifstream &file) {
  return file.peek() == std::ifstream::traits_type::eof();
}

std::string tolower(std::string in_str) noexcept;

void initialize_loggers(const std::string &log_file,
                        const std::string &output_file);
void update_loggers();
void finalize_loggers();

// renames the old output file if possible or choose a new output file name
void prepare_output_file(std::string &output_file);

// bool to yes/no conversion
std::string to_string(const bool &b);
