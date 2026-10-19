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
#include <optional>
#include <string>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "io/run_bundle.hpp"
#include "simulations/landscape/refinement_levels.hpp"
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

    // The refinement hierarchy (docs/architecture/multiscale-plan.md, section 3, "Fields
    // and units"; step S0.4): appended to the window block as level, factor,
    // metres_per_pixel and resampled, each only when set, so a bundle that sets none is
    // the one written before them. For a refinement (refinement_levels.hpp) the window
    // block's halo is the target's halo in global cells (level 1's); the deeper levels'
    // halos are larger (refinement_grid.hpp, levelhalo) and are not written.
    std::optional<int> level;              // the level the fields are on (k = factor^level)
    std::optional<int> factor;             // refinement.factor
    std::optional<double> metresperpixel;  // a resampled output's requested metres per pixel
    std::optional<bool> resampled;         // true when the output is not a level grid
    // A resampled output's pixel grid (refinement_levels.hpp, resampleregionalfields; step
    // S2.4): the fields are pixelcolumns x pixelrows instead of the window's fine grid, and
    // the window block gains pixel_columns, pixel_rows and pixel_cells (the source level's
    // cells per pixel edge, n). The pixels start at the kept box's north-west corner; a
    // partial pixel at the east and south edges is dropped, so they may cover less than the
    // box. All three or none.
    std::optional<int> pixelcolumns;
    std::optional<int> pixelrows;
    std::optional<double> pixelcells;
    // A refinement's per-level diagnostics (RefinementResult::levels, step S2.3): written as
    // the diagnostics block's "levels" array only when not empty.
    std::vector<landscape::LevelDiagnostics> levels;
};

// Throws std::invalid_argument when a field's size is not the kept window's fine grid (or
// the pixel grid when Metadata sets one), when the pixel grid is set in part, and as
// runbundle::Writer (an existing non-empty directory is refused).
void writeregionalbundle(const std::filesystem::path& directory, const landscape::RegionalFields& fields,
                         const Metadata& metadata);

// A cube window's regional bundle (docs/architecture/sphere-mesh.md section 7, grid kind
// cube-window; T8 W2-W1c): the fields on the kept level cells of `kept` (a
// mesh::cubewindowgrid; every field holds kept.cellcount() values in its local order, ascending
// level keys), shape [n]; the cell table under grid/ (cell_ids: the u8 level keys; centre,
// latitude_deg, longitude_deg, area_m2, neighbours as local ids with -1 outside the kept cells,
// reverse); and every grid field's lat-lon export of the window's box (fields.window's kept box
// at its scale k: finewidth x fineheight cells of the level's lat-lon grid, kW x kH) by the field's
// remap rule through mesh::makeremap on the window itself, a direction field as 0 (not
// recomputed); an export cell that overlaps no kept cell is NaN in a float field and 0 in an
// integer one (no extrapolation; the catalogue has no missing value), marked by the export-only
// mask geosphere.export_covered (u1, rule "coverage").
// The window block is the lat-lon one's (fields.window on the export grid W = 4 N).
// Throws std::invalid_argument when kept is not a cube window of N_L = basefacecells x scale,
// fields.window.worldwidth is not 4 basefacecells, a pixel grid is set, or a field's size is not
// kept.cellcount(); and as runbundle::Writer.
void writecubewindowbundle(const std::filesystem::path& directory, const landscape::RegionalFields& fields,
                           const mesh::CellGrid& kept, int basefacecells, int level, const Metadata& metadata);

// The ascending level keys of a cube-window bundle (grid/cell_ids.npy). Throws
// std::runtime_error for another grid kind or a table that does not hold cell_count keys.
std::vector<mesh::LevelKey> cubewindowkeys(const runbundle::Reader& reader);

// The window block of a manifest: scale, world grid, the kept box in global cells, the
// absolute fine origin of the kept grid on the kW x kH grid, the box's edges in degrees
// (west > east when it crosses the date line) and the halo.
std::string windowjson(const landscape::RegionalWindow& window);

// SHA-256 of a file's bytes, lower-case hex; throws std::runtime_error when unreadable.
std::string filesha256(const std::filesystem::path& path);
}
