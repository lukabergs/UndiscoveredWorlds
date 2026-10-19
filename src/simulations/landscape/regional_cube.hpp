#pragma once

// The regional pass and the refinement loop on cube windows (T8 W2-W2f; docs/architecture/sphere-mesh.md
// section 8, "Window grids", "Halo by rings" and "Cube rules of the regional units"): the steps of
// regionallevelsteps (regional_pass.hpp) in LevelStep order and finishregionallevel, each through its
// unit's cube form, on window grids of the key chain (refinement_grid.hpp, keyhalochain and
// keywindowgrid) instead of RegionalWindow boxes. The lat-lon pass and loop keep their own drivers
// (regional_pass.cpp, refinement_levels.cpp), so lat-lon windows are not changed by this unit.
// Every value of a kept cell is a pure function of the seed, its key, its centre and the global fields
// within the halo, so two cube windows give the same values on the cells they share (the overlap
// property; tests Landscape.RegionalCube*). Not ported: the tributaries (refused while
// regional.tributary_enabled is 1) and the resample to metres per pixel. Planet-free.

#include <cstddef>
#include <memory>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// A cube window's result: the kept cells' fields in ascending key order (fields.window is not set: the
// caller gives the lat-lon export box), those cells as a cube window grid, the per-level record.
struct CubeRegionResult
{
    RegionalFields fields;
    std::shared_ptr<const mesh::CellGrid> kept;
    // The refinement's levels 1 .. L (index level - 1; LevelDiagnostics::window is lat-lon and stays
    // empty but for level, factor and halo; domaincells holds the cube sizes); empty in the single pass.
    std::vector<LevelDiagnostics> levels;
    std::vector<std::size_t> keptcells, domaincells;  // per level 1 .. L (the single pass: one entry)
    std::vector<KeyLevelWindow> chain;                // the key chain, levels 0 .. L (empty in the single pass)
};

// The refinement loop on a cube world: inputs.grid the world's full cube level 0, kept0 the level-0 keys
// of the target's kept box (cellsinpixelbox), halo the target's halo in global cells (level 1's rings of
// level-0 cells), level >= 1. The chain is keyhalochain(*inputs.grid, factor, level, kept0, halo,
// prolongationreach, solvemargin(relief.catchmentreach)); level l's domain is keywindowgrid(chain[l]) and
// its kept cells the descendants of chain[l].keptglobal; level 1 projects onto the level-0 cells of
// chain[0], level l >= 2 onto level l - 1's kept cells. The order, the reaches and every rule are
// runrefinement's (refinement_levels.hpp) with the units' cube forms; the parent's masks (S2.8) are its
// kept cells' sea and lake surface at a cell's parent and the parent's 8 neighbours at the parent level
// (mesh::cubekeyneighbour). The solve pins the channels, the sea and the domain's outer ring (the cells
// with a table slot missing that the full level has). Throws std::invalid_argument for a grid that is not
// a full cube, an empty kept box, the halo below refinementhalo, the tributaries on, and as runrefinement.
CubeRegionResult runcuberefinement(const RegionalGlobalInputs& inputs, const std::vector<mesh::LevelKey>& kept0, int halo,
                                   int level, const RegionalParameters& parameters, const RefinementParameters& refinement);

// The single pass on a cube world at k fine cells per level-0 cell edge (any k >= 1): the domain is kept0
// dilated by `halo` rings of level-0 cells, its cells their descendants at N k; runregionalpass's steps
// (no smoothing, no projection, the unbounded solve). Throws std::invalid_argument as runregionalpass
// (halo below regionalhalo) and for the tributaries on.
CubeRegionResult runcuberegionalpass(const RegionalGlobalInputs& inputs, const std::vector<mesh::LevelKey>& kept0, int scale,
                                     int halo, const RegionalParameters& parameters);
}
