#pragma once

#include "social_generation_internal.hpp"

// Authored CSV definitions. Callers choose a root directory; absent or invalid
// files retain the existing profile defaults and empty optional tables.
namespace socialgen
{
bool loadsocialprofilecsv(const std::string& filepath, SocialProfile& profile);
void loadsocialtables(const std::string& rootpath, SocialTables& tables);
}
