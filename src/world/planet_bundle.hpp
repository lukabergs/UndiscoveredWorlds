#pragma once

#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "io/run_bundle.hpp"
#include "core/world/world_fields.hpp"

class planet;

// Last writer of each field or table, by step label; supplied by the scheduler.
struct BundleProvenance
{
    std::map<worldmodel::FieldId, std::string> fields;
    std::map<worldmodel::TableId, std::string> tables;
};

// A bundle a run starts from (--from-bundle). The world loads the fields it
// stores; tables are text and the climate solver's grids, levels and diagnostics
// belong to a solve, so those stay in the source and are carried into the
// bundles of the run unchanged unless one of its steps writes them.
struct PlanetBundleSource
{
    std::shared_ptr<const runbundle::Reader> reader;
    std::vector<worldmodel::FieldId> loadedfields;
    std::vector<worldmodel::TableId> loadedtables;      // equal to the cleared world's: nothing to load
    std::vector<worldmodel::FieldId> carriedfields;
    std::vector<worldmodel::TableId> carriedtables;
    BundleProvenance provenance;                        // the source's last writers
};

// Sizes the cleared world to the bundle's grid and loads its fields.
PlanetBundleSource loadplanetbundle(planet& world, const std::filesystem::path& directory);

// Writes every catalogue field and entity table of the world as a run bundle.
// With a source, its carried fields and tables are written from it except the
// ones in the written lists (what the steps of this run wrote).
void writeplanetbundle(const planet& world, const std::filesystem::path& directory,
    const runbundle::Metadata& metadata, const BundleProvenance& provenance = {},
    const PlanetBundleSource* source = nullptr, const std::vector<worldmodel::FieldId>& writtenfields = {},
    const std::vector<worldmodel::TableId>& writtentables = {});
