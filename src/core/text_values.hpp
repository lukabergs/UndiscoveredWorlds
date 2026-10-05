#pragma once

// Text forms of values: numbers with separators, values read from saved text, direction names.

#include <string>
#include <vector>

#include "core/classes.hpp"

void toggle(bool &val);
std::string formatnumber(int val);
bool stob(std::string const& instring);
short stos(std::string const& instring);
unsigned short stous(std::string const& instring);
char stoc(std::string const& instring);
unsigned char stouc(std::string const& instring);
std::string getdirstring(int dir);
