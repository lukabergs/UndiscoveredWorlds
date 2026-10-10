#pragma once

// The solve step of a level (docs/architecture/multiscale-plan.md, §3 "Phase S3 re-planned";
// the hook written in S3.0, regional_pass.hpp): the analytical relief (analytical_relief.hpp)
// from the pinned channels, the sea and the domain's outer ring, as the regional pass ran it
// since Phase E (erosion plan E3) with R9's uplift term. S3.2 adds, in the refinement loop only,
// the solve against its parent (F* = Z_l + prolong(lift_(l-1)) on unpinned cells,
// refinement.parent_lift) and a channel cell's elevation as its cell mean b + min(h, theta c_l)
// (refinement.channel_cell_slope theta; 0: the bed), the bed kept in LevelState::bed_m.
// Planet-free.

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
// elsewhere), state.pinned, the solve's diagnostics (reliefconverged, changedreceivers;
// meanlowering_m and maxlowering_m over the kept window's unpinned cells, the solve's own
// lowering) and state.solveddischarge_m3_yr. Then, in the refinement loop (context.refinement set):
// - at level >= 2 with refinement.parent_lift 1 or 2 (S3.2, the solve against its parent): the
//   parent's lift prolonged as the start was (prolong, refinement_levels.hpp, with the limit when
//   refinement.prolongation_limit is 1), the lift being the parent's LevelState::lift_m (1) or its
//   final elevation less its ParentFields::solved_m (2), is taken off the solve's start on the
//   unpinned cells (the pinned beds keep theirs: de-lifting them measured worse at S3.2's review), so
//   the land beyond the knickpoints' reach is not lifted twice, and added back to their result, F* = Z_l
//   + lift (pinned cells keep state.surface_m); a negative lift takes no cell below the lowest of its solve,
//   the next float above the sea level and its parent cell's elevation. So the level's field keeps
//   what the parent's smoothing and projection (and, with 2, its forms, lakes and its own added
//   lift) put on its solve, and the projection only removes about Z_(l-1) - mean(Z_l);
//   meanlowering_m and maxlowering_m are then the start surface less F*;
// - with refinement.channel_cell_slope theta > 0: every land channel cell (not sea) takes its cell
//   mean b + min(h, theta c_l), b its bed (state.bed_m, unchanged), h the mean height above b of its
//   non-channel, non-sea 8-neighbours in the domain after the lift (at least 0; 0 without one), c_l
//   the level's north-south cell size; so b <= elevation <= b + theta c_l on land channel cells.
// state.solved_m is the result without the parent's lift (the channel cells' means kept): the next
// level's ParentFields::solved_m. Outside the loop, and with both keys 0, the solve as before bit for
// bit (state.solved_m a copy of the elevation). Every rule a pure function of the level's fields and
// the parent's kept box. Throws as solveanalyticalrelief and as prolong.
void solvelevel(const LevelContext& context, LevelState& state);
}
