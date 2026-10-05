#pragma once

// The tunable values of every simulation module, as a run sets and records them.
// Defaults are compiled in from configs/parameters/<module>.json.

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "core/parameters.hpp"

// Every module's entries, in module order.
const std::vector<parameters::Entry>& recipeparameters();

// "name=value", the form of `world_pipeline run --set`.
bool setrecipeparameter(std::string_view assignment, std::string* error = nullptr);

// A JSON file of values: objects nested as the names are dotted
// ({"terrain": {"oceanridges": {"grid_size": 16}}}), or the "parameters" object of a bundle manifest.
bool loadrecipeparameters(const std::filesystem::path& path, std::string* error = nullptr);

// Conditions between values that the code relies on; checked once the overrides are in.
bool validaterecipeparameters(std::string* error = nullptr);

// Every value as one JSON object keyed by full name, for the bundle manifest.
std::string recipeparametersjson();
