#pragma once

// Tunable values of the simulation modules, by name. Each module's values are
// variables declared in its generated <module>_parameters.hpp (scripts/parameters.py,
// from configs/parameters/<module>.json); an Entry names one of them so that a run
// can set it from text and record it.

#include <string>
#include <string_view>
#include <variant>

namespace parameters
{
struct Entry
{
    const char* name;   // <module>.<group>.<key>
    std::variant<bool*, short*, int*, float*, double*> value;
};

// Sets the value from text ("true"/"false"/"on"/"off" for a bool, a number otherwise).
bool settext(const Entry& entry, std::string_view text, std::string* error = nullptr);

// The value as JSON text: the shortest number that reads back to the same value.
std::string text(const Entry& entry);
}
