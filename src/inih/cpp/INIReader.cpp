// Read an INI file into easy-to-access name/value pairs.

// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
// https://github.com/benhoyt/inih

// NOTE: the "sections" parameter has been removed from the interface functions

#include "../ini.h"
#include "INIReader.h"



using std::string;

INIReader::INIReader(const string& filename) {
	_error = ini_parse(filename.c_str(), ValueHandler, this);
}

int INIReader::ParseError() const noexcept {
	return _error;
}

string INIReader::Get(const string& name, const string default_value) const {

	const string key = "=" + tolower(name);
	return _values.count(key) ? _values.find(key)->second : default_value;
}

void INIReader::dump_compilation_info() const {
	auto log = spdlog::get("loggers");
	log->debug("------------compilation machine------------");
	#if defined(BOOST_ARCH_X86_64_AVAILABLE)
		log->debug("Architecture: {}", BOOST_ARCH_X86_64_NAME);
	#endif
	#if defined(BOOST_HW_SIMD_AVAILABLE)
		string simd_features = "";
		if (BOOST_HW_SIMD_X86 >= BOOST_HW_SIMD_X86_SSE3_VERSION) { simd_features = "SSE3 "; }
		if (BOOST_HW_SIMD_X86 >= BOOST_HW_SIMD_X86_SSE4_1_VERSION) { simd_features = "SSE4.1 "; }
		if (BOOST_HW_SIMD_X86 >= BOOST_HW_SIMD_X86_SSE4_2_VERSION) { simd_features = "SSE4.2 "; }
		if (BOOST_HW_SIMD_X86 >= BOOST_HW_SIMD_X86_AVX_VERSION) { simd_features += "AVX "; }
		if (BOOST_HW_SIMD_X86 >= BOOST_HW_SIMD_X86_FMA3_VERSION) { simd_features += "FMA "; }
		if (BOOST_HW_SIMD_X86 >= BOOST_HW_SIMD_X86_AVX2_VERSION) { simd_features += "AVX2 "; }
		log->debug("SIMD features: {}", simd_features);
	#endif

	#if defined(BOOST_OS_WINDOWS_AVAILABLE)
		log->debug("OS: {}", BOOST_OS_WINDOWS_NAME);
	#endif
	#if defined(BOOST_OS_LINUX_AVAILABLE)
		log->debug("OS: {}", BOOST_OS_LINUX_NAME);
	#endif

	#if defined(BOOST_COMP_INTEL_AVAILABLE) 
		log->debug("Available compiler: {} version {}", BOOST_COMP_INTEL_NAME, __INTEL_COMPILER);
	#endif
	#if defined(BOOST_COMP_GNUC_AVAILABLE) 
		log->debug("Available compiler: {} version {}", BOOST_COMP_GNUC_NAME, __GNUC__);
	#endif
	#if defined(BOOST_COMP_CLANG_AVAILABLE) 
		log->debug("Available compiler: {} version {}", BOOST_COMP_CLANG_NAME, __clang__);
	#endif
	#if defined(BOOST_COMP_MSVC_AVAILABLE)
		log->debug("Available compiler: {} version {}", BOOST_COMP_MSVC_NAME, _MSC_VER);
	#endif

		log->debug("Compilation date: {} {}", __DATE__ , __TIME__);

}

void INIReader::dump_env_info() const {
	auto log = spdlog::get("loggers");
	log->debug("-----------enviroment variables------------");
	const vector<string> env_variables = { "OMP_DYNAMIC", "OMP_SCHEDULE","OMP_NUM_THREADS", "MKL_NUM_THREADS", "KMP_AFFINITY", "OMP_PROC_BIND", "OMP_PLACES", "GOMP_CPU_AFFINITY" };
	const vector<string> slurm_vars = { "SLURM_JOB_ID", "SLURM_SUBMIT_DIR", "SLURM_NTASKS", "SLURM_JOB_NODELIST" };
	const vector<string> pbs_vars = { "PBS_JOBID", "PBS_O_WORKDIR", "PBS_NP", "PBS_NODEFILE" };

	if (const char* env_p = getenv(slurm_vars.at(0).c_str())) {
		for (const auto &var : slurm_vars) {
			if (const char* env_p = getenv(var.c_str())) {
				log->debug(">> {}={}", var, env_p);
			}
		}
	}

	if (const char* env_p = getenv(pbs_vars.at(0).c_str())) {
		for (const auto &var : pbs_vars) {
			if (const char* env_p = getenv(var.c_str())) {
				log->debug(">> {}={}", var, env_p);
			}
		}
	}

	for (const auto &var : env_variables) {
		if (const char* env_p = getenv(var.c_str())) {
			log->debug(">> {}={}", var, env_p);
		}
	}

	auto start_date_time = chrono::system_clock::now();
	auto tt = chrono::system_clock::to_time_t(start_date_time);
	auto timeinfo = localtime(&tt);
	char time_string_buffer[80];
	strftime(time_string_buffer, 80, "System clock: %F %T", timeinfo);
	log->debug(time_string_buffer);
}

void INIReader::dump_parsed() const {

	auto log = spdlog::get("loggers");
	auto output_log = spdlog::get("output");
	sort(_parsed.begin(), _parsed.end(), [](const auto& lhs, const auto& rhs) {
		return tolower(rhs.at(0)) > tolower(lhs.at(0));
	});
	update_loggers();
	dump_compilation_info();
	dump_env_info();

	log->info("-------------slabcc parameters-------------");
	output_log->info("# Parameters read from the file or their default values:");
	for (const auto &i : _parsed) {
		log->info(i.at(0) + " = " + i.at(1));
		output_log->info(i.at(0) + " = " + i.at(1));
	}
	log->info("-----------------------------------------");
	log->flush();

	for (const auto &i : _error_msgs) {
		log->error(i);
	}
	
	// add the deprectated parameters here!
	const vector<string> deprecated_params{"optimize_charge"};

	for (const auto &val : _values) {
		string param_in_file = val.first;
		// get rid of the "=" at the start
		param_in_file.erase(param_in_file.begin());
		bool has_parsed = false;
		for (auto const &param_parsed : _parsed) {
			if (tolower(param_parsed.at(0)) == param_in_file) {
				has_parsed = true;
			}
		}
		if (!has_parsed) {
			if (find(deprecated_params.begin(), deprecated_params.end(), param_in_file) != deprecated_params.end()) {
				log->warn("The parameter \"{}\" in deprecated! Please refer to this version of the slabcc's manual for a complete list of the supported parameters.", param_in_file);
			}
			else {
				log->warn("Unrecognized parameter in the input file: " + param_in_file);
			}
		}
	}
}

string INIReader::GetStr(const string& name, const string& default_value) const
{
	string valstr = Get(name, default_value);
	valstr.erase(std::remove(valstr.begin(), valstr.end(), '"'), valstr.end());
	replace(valstr, "\n", " ");
	_parsed.push_back({ name, valstr });
	return valstr;
}

arma::rowvec INIReader::GetVec(const string& name, const arma::rowvec default_value) const
{
	arma::rowvec result;
	string error = "";
	string valstr = Get(name);
	replace(valstr, "\n", " ");
	replace(valstr, ";", " ");
	try {
		result = arma::rowvec(valstr);
	}
	catch (const exception &e) {
		error = string(e.what());
	}
	if (error != "") {
		_error_msgs.push_back("Error in parsing \"" + name + "\": " + error);
	}

	const arma::rowvec out = result.is_empty() ? default_value : result;
	_parsed.push_back({ name, to_string(out) });
	return out;
}

arma::mat INIReader::GetMat(const string& name, const arma::mat default_value) const
{
	auto log = spdlog::get("loggers");
	string error = "";
	string valstr = Get(name);
	replace(valstr, "\n", ";");
	replace(valstr, ";;", ";");
	arma::mat result;
	try {
		 result = arma::mat(valstr);
	}
	catch (const exception &e) {
		error = string(e.what());
	}
	if (error != "") {
		_error_msgs.push_back("Error in parsing \"" + name + "\": " + error);
	}
	const arma::mat out = result.is_empty() ? default_value : result;
	_parsed.push_back({ name, to_string(out) });
	return out;
}



long INIReader::GetInteger(const string& name, const long default_value) const
{
	const string valstr = Get(name);
	const char* value = valstr.c_str();
	char* end;
	const long n = strtol(value, &end, 0);
	const long result = end > value ? n : default_value;
	_parsed.push_back({ name, std::to_string(result) });
	return result;
}

double INIReader::GetReal(const string& name, const double default_value) const
{
	const string valstr = Get(name);
	const char* value = valstr.c_str();
	char* end;
	const double n = strtod(value, &end);
	const double result = end > value ? n : default_value;
	_parsed.push_back({ name, ::to_string(result) });
	return result;
}

bool INIReader::GetBoolean(const string& name, const bool default_value) const {
	const string valstr = tolower(Get(name));
	bool result = default_value;

	if (valstr == ".true." || valstr == "true" || valstr == "yes" || valstr == "on" || valstr == "1")
		result = true;
	else if (valstr == ".false." || valstr == "false" || valstr == "no" || valstr == "off" || valstr == "0")
		result = false;
	else if (valstr != "")
		_error_msgs.push_back("Error in parsing \"" + name + "\"");

	_parsed.push_back({ name, to_string(result) });
	return result;
}

int INIReader::ValueHandler(void* user, const char* section, const char* name,
	const char* value)
{
	INIReader* reader = static_cast<INIReader*>(user);
	const string key = "=" + tolower(string(name));
	if (reader->_values[key].size() > 0)
		reader->_values[key] += "\n";
	reader->_values[key] += value;
	return 1;
}

void INIReader::replace(std::string& str, const std::string& from, const std::string& to) const {
	
	size_t start_pos = str.find(from);
	while (start_pos != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos = str.find(from);
	}
}