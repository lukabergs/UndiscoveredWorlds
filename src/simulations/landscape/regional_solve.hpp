#pragma once

// The solve step of a level (docs/architecture/multiscale-plan.md, §3 "Phase S3 re-planned";
// the hook written in S3.0, regional_pass.hpp): the analytical relief (analytical_relief.hpp)
// from the pinned channels, the sea and the domain's outer ring, as the regional pass ran it
// since Phase E (erosion plan E3) with R9's uplift term. S3.2 changes it: the solve against its
// parent (F* = Z_l + prolong(lift_(l-1)) on unpinned cells) and a channel cell's elevation as its
// cell mean b + min(h, theta c_l) (refinement.channel_cell_slope theta; 0: the bed, as here), the
// bed kept in LevelState::bed_m. Planet-free.

namespace landscape
{
struct LevelContext;  // regional_pass.hpp
struct LevelState;

// LevelStep::solve. On the window's domain: the solve's start is state.surface_m (the start
// surface with the land channel beds pinned) with every sea cell at the sea level; pinned are
// the channel cells, the sea cells and the domain's outer ring; the runoff is the global runoff
// upsampled (mm/yr to m/yr, at least 0); the key the cell's absolute fine index (row x the
// planet's fine width + column); with the incision active, the uplift term the incision rate
// upsampled (at least 0). Writes state.elevation_m (state.surface_m on pinned cells, the solve
// elsewhere), state.pinned, state.solved_m (a copy of that elevation), the solve's diagnostics
// (reliefconverged, changedreceivers; meanlowering_m and maxlowering_m over the kept window's
// unpinned cells) and state.solveddischarge_m3_yr. Throws as solveanalyticalrelief.
void solvelevel(const LevelContext& context, LevelState& state);
}
