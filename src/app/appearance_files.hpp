#pragma once

// Saving and loading the map appearance settings (.uws files) of the desktop application.

#include <string>
#include <vector>

#include "core/classes.hpp"

class planet;
class region;

void savesettings(planet& world, std::string filename);
bool loadsettings(planet& world, std::string filename);
