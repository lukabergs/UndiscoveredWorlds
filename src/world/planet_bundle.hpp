#pragma once

#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "io/run_bundle.hpp"
#include "core/world/world_fields.hpp"
#include "world/planet_export.hpp"

class planet;

// What a cube world's bundle takes from outside the world for its lat-lon export
// (sphere-mesh.md sections 6, 7 and 10; ignored on a lat-lon world). Supplied by the
// pipeline (export_routing.hpp, bundleexport).
struct PlanetBundleExport
{
    // Recomputes the `recomputed` fields on the export (the routing's directions and
    // accumulated fluxes); their exports are the recomputed values, marked recomputed.
    // Without it direction fields are exported as 0.
    LatlonDirectionRecompute recompute;
    std::vector<worldmodel::FieldId> recomputed;
    // The climate adapter's lat-lon shadow (of the export's grid) and the fields it owns.
    // A shadow-owned field on the shadow's grid is exported from it; one off it (a zonal
    // profile, the solver's grids, the diagnostics) is stored from it in its lat-lon
    // shape, as on a lat-lon world. The climate tables (atmosphere levels, weather
    // state) are the shadow's too.
    const planet* climateshadow = nullptr;
    std::function<bool(worldmodel::FieldId)> shadowowned;
};

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
// ones in the written lists (what the steps of this run wrote). A cube world's
// lat-lon export takes the recomputed fields and the shadow's from `exports`.
void writeplanetbundle(const planet& world, const std::filesystem::path& directory,
    const runbundle::Metadata& metadata, const BundleProvenance& provenance = {},
    const PlanetBundleSource* source = nullptr, const std::vector<worldmodel::FieldId>& writtenfields = {},
    const std::vector<worldmodel::TableId>& writtentables = {}, const PlanetBundleExport& exports = {});
