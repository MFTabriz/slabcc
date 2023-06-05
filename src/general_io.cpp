// Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz
// Copyrights licensed under the 2-Clause BSD License.
// See the accompanying LICENSE.txt file for terms.

#include "general_io.hpp"

void cli_params::parse(int argc, char *argv[]) {

  bool showHelp = false, showManual = false, showVer = false, showAttr = false;
  auto cli =
      clara::Help(showHelp) |
      clara::Opt(input_file,
                 "input_file")["-i"]["--input"]("slabcc input file name") |
      clara::Opt(output_file,
                 "output_file")["-o"]["--output"]("slabcc output file name") |
      clara::Opt(log_file, "log_file")["-l"]["--log"]("slabcc log file name") |
      clara::Opt(diff_only)["-d"]["--diff"](
          "calculate charge and the potential differences only") |
      clara::Opt(showManual)["-m"]["--man"]("show quick start guide") |
      clara::Opt(showVer)["-v"]["--version"]("show slabcc version") |
      clara::Opt(showAttr)["-c"]["--copyright"]("show copyright information");

  auto cli_result = cli.parse(clara::Args(argc, argv));
  if (!cli_result) {
    std::cerr << "Error in command line: " << cli_result.errorMessage() << '\n';
    exit(1);
  }

  if (showHelp) {
    std::ostringstream oss;
    oss << (clara::Parser() | cli);
    std::cout << oss.str();
    exit(0);
  }

  if (showVer) {
    std::cout << "SLAB Charge Correction (slabcc) (version "
              << SLABCC_VERSION_MAJOR << "." << SLABCC_VERSION_MINOR << "."
              << SLABCC_VERSION_PATCH << ")\n";
#include "gitinfo.inc"
    std::cout << git_info << "\n\n";
    std::cout << "Included libraries: " << '\n';
    std::cout << "> Armadillo " << ARMA_VERSION_MAJOR << "."
              << ARMA_VERSION_MINOR << "." << ARMA_VERSION_PATCH << '\n';
    std::cout << "> NLOPT " << nlopt::version_major() << "."
              << nlopt::version_minor() << "." << nlopt::version_bugfix()
              << '\n';
    std::cout << "> SPDLOG " << SPDLOG_VER_MAJOR << "." << SPDLOG_VER_MINOR
              << "." << SPDLOG_VER_PATCH << '\n';
    exit(0);
  }

  if (showManual) {
    std::cout
        << "Quick start guide:\n"
           "To calculate the charge correction, slabcc needs the following "
           "files:\n"
           " 1. Input parameters file (default: slabcc.in)\n"
           " 2. CHGCAR of the neutral system (default: CHGCAR.N)\n"
           " 3. CHGCAR of the charged system (default: CHGCAR.C)\n"
           " 4. LOCPOT of the neutral system (default: LOCPOT.N)\n"
           " 5. LOCPOT of the charged system (default: LOCPOT.C)\n\n"
           "The input parameters file for a slab should minimally include:\n"
           " - charge_position : approximate position of the localized charge "
           "center\n"
           " - diel_in : dielectric constant/tensor of the slab\n"
           " - normal_direction : direction normal to the surface\n"
           " - interfaces : surface positions of the slab\n\n"
           "Complete manual: https://meisam.codeberg.page/slabcc";
    exit(0);
  }

  if (showAttr) {
    std::cout
        << "Copyright (c) 2018-2023, Bremen Center for Computational Materials "
           "Science (BCCMS), M. Farzalipour Tabriz\n"
           "The source code and all the documentations are available under The "
           "2-Clause BSD License. For more information see LICENSE.TXT.\n\n"
           "Included libraries\n\n"
           "-Armadillo C++ Linear Algebra Library: licensed under the Apache "
           "License 2.0\n"
           "  Copyright 2008-2023 Conrad Sanderson\n"
           "  Copyright 2008-2016 National ICT Australia (NICTA)\n"
           "  Copyright 2017-2023 Data61 / CSIRO\n"
           "  This product includes software developed by Conrad Sanderson\n"
           "  This product includes software developed at National ICT "
           "Australia (NICTA)\n"
           "  This product includes software developed at Data61 / CSIRO\n\n"
           "-inih(INI Not Invented Here): licensed under the 3-clause BSD "
           "license\n"
           "  (c) 2009, Ben Hoyt, et al.\n\n"
           "-clara: licensed under the Boost Software License 1.0\n"
           "  (c) 2017-2018, Phil Nash, Martin Horenovsky, et al.\n\n"
           "-Spline (Cubic Spline Interpolation): licensed under the Beer-Ware "
           "License 42\n"
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

unsigned int xyz2int(const std::string &s) {
  const char dir_char = std::tolower(s.at(0));
  const std::unordered_map<char, unsigned int> dir_map{
      {'0', 0}, {'x', 0}, {'a', 0}, {'1', 1}, {'y', 1},
      {'b', 1}, {'2', 2}, {'z', 2}, {'c', 2},
  };
  return dir_map.at(dir_char);
}

char int2xyz(const unsigned int &i) {
  const std::unordered_map<unsigned int, char> dir_map{
      {0, 'X'},
      {1, 'Y'},
      {2, 'Z'},
  };
  return dir_map.at(i);
}

std::string tolower(std::string in_str) noexcept {
  for (char &c : in_str) {
    c = tolower(c);
  }
  return in_str;
}

void prepare_output_file(std::string &output_file) {
  if (file_exists(output_file)) {
    const int max_file_tests = 100000;
    int counter = 1;
    std::string backup_name = output_file + ".old" + to_string(counter);

    while (file_exists(backup_name)) {
      counter++;
      backup_name = output_file + ".old" + to_string(counter);
      if (counter > max_file_tests) {
        std::cerr << "ERROR: no suitable name was found for backing-up the old "
                     "output file: "
                  << output_file << std::endl;
        break;
      }
    }

    if (std::rename(output_file.c_str(), backup_name.c_str())) {
      // if we don't have the permission to rename the existing file
      counter = 1;
      std::string new_name = output_file + ".new" + to_string(counter);
      while (file_exists(new_name)) {
        counter++;
        new_name = output_file + ".new" + to_string(counter);
        if (counter > max_file_tests) {
          // cleanup your old files
          std::cerr << "ERROR: output file '" << output_file
                    << "' exists and no replacement was found!" << std::endl;
          std::exit(2);
        }
      }
      output_file = new_name;
    } else {
      std::cout << "The old output file '" << output_file << "' was renamed to "
                << backup_name << std::endl;
    }
  }
}

void initialize_loggers(const std::string &log_file,
                        const std::string &output_file) {

  const std::string tmp_file = "slabcc.tmp";
  std::vector<spdlog::sink_ptr> sinks;
  sinks.push_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
  sinks.push_back(
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true));
  sinks.push_back(
      std::make_shared<spdlog::sinks::basic_file_sink_mt>(tmp_file, true));
  auto combined_logger =
      std::make_shared<spdlog::logger>("loggers", begin(sinks), end(sinks));
  sinks.at(2)->set_level(spdlog::level::warn);
  sinks.at(2)->set_pattern("[%l] %v");
  combined_logger->set_pattern("%v");
  combined_logger->set_level(spdlog::level::info);
  spdlog::register_logger(combined_logger);

  auto output_logger = spdlog::basic_logger_mt("output", output_file);
  output_logger->set_pattern("%v");
  output_logger->set_level(spdlog::level::info);
}
void finalize_loggers() {
  auto log = spdlog::get("loggers");
  auto output_log = spdlog::get("output");
  const std::string tmp_file = "slabcc.tmp";
  log->flush();
  output_log->flush();

  std::ifstream messages_list(tmp_file);
  std::string tmp_line = "";
  if (!file_is_empty(messages_list)) {
    output_log->info("\n[Messages]");
    while (getline(messages_list, tmp_line)) {
      output_log->info(tmp_line);
    }
  }
  messages_list.close();
  std::remove(tmp_file.c_str());
  output_log->flush();
}
void update_loggers() {

  auto log = spdlog::get("loggers");
  std::string log_pattern = "%v";

  if (is_active(verbosity::info)) {
    log->set_level(spdlog::level::info);
  }
  if (is_active(verbosity::debug)) {
    log->set_level(spdlog::level::debug);
    log_pattern = "[%^%l%$] " + log_pattern;
  }
  if (is_active(verbosity::trace)) {
    log->set_level(spdlog::level::trace);
  }

  if (is_active(verbosity::timing)) {
    log_pattern = "[%K]" + log_pattern;
  }

  log->set_pattern(log_pattern);
  log->sinks().at(2)->set_pattern("[%^%l%$] %v");
}

std::string to_string(const bool &b) { return b ? "yes" : "no"; }
