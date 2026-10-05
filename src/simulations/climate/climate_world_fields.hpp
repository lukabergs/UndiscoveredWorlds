#pragma once

// Catalogue access to climateworld: the climate-owned fields and, for a climate
// world that owns its terrain, the bound terrain and planet inputs.

#include <string>
#include <vector>

#include "core/world/field_encoding.hpp"

class climateworld;

// False when the field is not available from the climate world.
bool readclimateworldfield(const climateworld& world, worldmodel::FieldId id, worldmodel::FieldData& data, bool includeinputs);
// Loads a climate field stored on the world grid. False for the solver's grids,
// levels and diagnostics, which only a solve in the same run provides.
bool loadclimateworldfield(climateworld& world, const worldmodel::FieldData& data);
bool readclimateworldtable(const climateworld& world, worldmodel::TableId id, std::string& text, std::string& extension);
// Diagnostics of the last solve (climatediagnostics); absent until a solve records them.
bool isclimatediagnosticfield(worldmodel::FieldId id);
// Fields available from a climate world, in catalogue order.
std::vector<worldmodel::FieldId> climateworldfields(bool includeinputs);
