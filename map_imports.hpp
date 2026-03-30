#pragma once

#include <string>
#include <vector>

#include "planet.hpp"

bool importlandheightmap(planet& world, const std::string& filepathname);
bool importseadepthmap(planet& world, const std::string& filepathname);
bool importmountainmap(planet& world, const std::string& filepathname, std::vector<std::vector<bool>>& okmountains);
bool importvolcanomap(planet& world, const std::string& filepathname);
