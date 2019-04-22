// Read an INI file into easy-to-access name/value pairs.

// inih and INIReader are released under the New BSD license (see LICENSE.txt).
// Go to the project home page for more info:
//
// https://github.com/benhoyt/inih

#ifndef __INIREADER_H__
#define __INIREADER_H__

#include <boost/predef.h>
#include <cctype>
#include <cstdlib>
#include <armadillo>
#include "general_io.hpp"
#include "stdafx.h"

// Read an INI file into easy-to-access name/value pairs. 
class INIReader
{
public:
	// Construct INIReader and parse given filename. See ini.h for more info
	// about the parsing.
	INIReader(const std::string& filename);

	// Return the result of ini_parse(), i.e., 0 on success, line number of
	// first error on parse error, or -1 on file open error.
	int ParseError() const noexcept;

	std::string Get(const std::string& name, const std::string default_value = "") const;

	// Get a string value from INI file, returning default_value if not found.
	std::string GetStr(const std::string& name, const std::string& default_value) const;

	// Get a row vector value from INI file as a set of doubles seperated by space,
	// returning default_value if not found.
	arma::rowvec GetVec(const std::string& name, const arma::rowvec default_value) const;

	// Get a matrix value from INI file as a set of doubles seperated by space for 
	// each column and seperated by ; or \n for each row, returning default_value if not found.
	arma::mat GetMat(const std::string& name, const arma::mat default_value) const;

	// Get an integer (long) value from INI file, returning default_value if
	// not found or not a valid integer
	long GetInteger(const std::string& name, const long default_value) const;

	// Get a real (floating point double) value from INI file, returning
	// default_value if not found or not a valid floating point value
	// according to strtod().
	double GetReal(const std::string& name, const double default_value) const;

	// Get a boolean value from INI file, returning default_value if not found or if
	// not a valid true/false value. Valid true values are "true", "yes", "on", "1",
	// and valid false values are "false", "no", "off", "0" (not case sensitive).
	bool GetBoolean(const std::string& name, const bool default_value) const;

	//writes the parsed variables to the output file and also log
	void dump_all() const;
	void dump_compilation_info() const;
	void dump_env_info() const;

	//replace all occurrences of the "from" string with "to" string inside "str" recursively
	void replace(std::string& str, const std::string from, const std::string to) const;

protected:
	mutable std::vector<std::vector<std::string>> _parsed;
private:
	int _error;
	std::map<std::string, std::string> _values;
	static int ValueHandler(void* user, const char* section, const char* name,
		const char* value);

};



#endif  // __INIREADER_H__
