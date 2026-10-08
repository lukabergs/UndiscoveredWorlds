#pragma once

// Regional bundles (docs/architecture/erosion-hydrology-plan.md, Phase E, step E2d, and
// section 6, the user's decisions of 2026-10-08, evening): the fine fields of a regional
// pass written as a run bundle (schema uw-run-bundle/v1, runbundle::Writer) whose grid is
// the window's fine grid (finewidth x fineheight, cell-centred, north row first). Fields
// take the world catalogue's ids, element types, layouts and units where the meaning is
// the same; two have new ids: the sea mask geosphere.sea (no catalogue grid holds one) and
// hydrosphere.channel_direction (the pinned channels only, where a world's river_direction
// covers every land cell). The manifest adds three blocks: window, source and
// diagnostics; every step record lists the fields written. Planet-free.

#include <filesystem>
#include <string>

#include "io/run_bundle.hpp"
#include "simulations/landscape/regional_pass.hpp"

namespace regionalbundle
{
struct Metadata
{
    // Seed, recipe, command, steps, properties and parameters; the grid is set from the
    // window and the window, source and diagnostics blocks are appended to its blocks.
    runbundle::Metadata run;
    std::string sourcebundle;          // the world bundle, as given on the command line
    std::string sourcemanifestsha256;  // SHA-256 of its manifest.json, lower-case hex
    long sourceseed = 0;
    int sealevel = 0;                  // the world's geosphere.sea_level, stored as that scalar
    std::string writtenby = "region";  // written_by of every field
};

// Throws std::invalid_argument when a field's size is not the kept window's fine grid,
// and as runbundle::Writer (an existing non-empty directory is refused).
void writeregionalbundle(const std::filesystem::path& directory, const landscape::RegionalFields& fields,
                         const Metadata& metadata);

// The window block of a manifest: scale, world grid, the kept box in global cells, the
// absolute fine origin of the kept grid on the kW x kH grid, the box's edges in degrees
// (west > east when it crosses the date line) and the halo.
std::string windowjson(const landscape::RegionalWindow& window);

// SHA-256 of a file's bytes, lower-case hex; throws std::runtime_error when unreadable.
std::string filesha256(const std::filesystem::path& path);
}
