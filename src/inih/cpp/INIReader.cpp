// Read an INI file into easy-to-access name/value pairs.

// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
// https://github.com/benhoyt/inih

// Modifications: sections parameter is removed from the interface fucntions

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
	// Use _values.find() here instead of _values.at() to support pre C++11 compilers
	return _values.count(key) ? _values.find(key)->second : default_value;
}

void INIReader::dump_parsed(std::ofstream& out_file) const {

	std::sort(_parsed.begin(), _parsed.end(), [](const auto& lhs, const auto& rhs) {
		return tolower(rhs.at(0)) > tolower(lhs.at(0));
	});

	if (is_active(verbosity::steps)) {
		cout << timing() << "-------------slabcc parameters-------------" << endl;
		for (const auto &i : _parsed) {
			std::cout << timing() << i.at(0) << " = " << i.at(1) << std::endl;
		}
		std::cout << timing() << "-----------------------------------------" << std::endl;
	}

	out_file << "# Parameters read from the file or their default values:" << std::endl;
	for (const auto &i : _parsed) {
		out_file << i.at(0) << " = " << i.at(1) << std::endl;
	}

	out_file.flush();

	for (const auto &val : _values) {
		string param_file = val.first;
		// get rid of the "=" at the start
		param_file.erase(param_file.begin());
		bool has_parsed = false;
		for (auto const &param_parsed : _parsed) {
			if (tolower(param_parsed.at(0)) == param_file) {
				has_parsed = true;
			}
		}
		if (!has_parsed) {
			cout << endl << timing() << ">> WARNING <<: Unrecognized parameter in the input file: " << param_file << endl << endl;
		}
	}
}

string INIReader::GetStr(const string& name, const string& default_value) const
{
	string valstr = Get(name, default_value);
	valstr.erase(std::remove(valstr.begin(), valstr.end(), '"'), valstr.end());
	_parsed.push_back({ name, valstr });
	return valstr;
}

arma::rowvec INIReader::GetVec(const string& name, const arma::rowvec default_value) const
{
	string valstr = Get(name);
	std::replace(valstr.begin(), valstr.end(), '\n', ' ');
	std::replace(valstr.begin(), valstr.end(), ';', ' ');
	const arma::rowvec result(valstr);
	const arma::rowvec out = result.is_empty() ? default_value : result;
	_parsed.push_back({ name, to_string(out) });
	return out;
}

arma::mat INIReader::GetMat(const string& name, const arma::mat default_value) const
{
	string valstr = Get(name);
	std::replace(valstr.begin(), valstr.end(), '\n', ';');
	const arma::mat result(valstr);
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
	_parsed.push_back({ name, to_string(result) });
	return result;
}

bool INIReader::GetBoolean(const string& name, const bool default_value) const {
	// Convert to lower case to make string comparisons case-insensitive
	const string valstr = tolower(Get(name));
	bool result = default_value;

	if (valstr == "true" || valstr == "yes" || valstr == "on" || valstr == "1")
		result = true;
	else if (valstr == "false" || valstr == "no" || valstr == "off" || valstr == "0")
		result = false;

	_parsed.push_back({ name, std::to_string(result) });
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

std::string INIReader::to_string(const double& d) {
	std::ostringstream output;
	output << std::setprecision(12) << d;
	return output.str();
}
