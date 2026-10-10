#pragma once

// Tributaries per level (docs/architecture/multiscale-plan.md, §3 "Kernels" and "Phase S3
// re-planned", step S3.3; the hook written in S3.0, inert: no tributaries). Contract for S3.3
// (Génevaux et al. 2013's growth from a network, Hack-law lengths, junction angles, Strahler
// order; §6 "Learned in S3" for the numbers): tributaries are grown from the level's pinned
// channels (state.channels after rasterisechannels), each hashed by the absolute key of its
// junction cell (row x the level's width + column, 64 bits) and the seed, so two windows draw it
// alike; a tributary's length lies in [2 level cells, one parent cell) and the tributary stays
// inside its junction's parent cell and that cell's 8 neighbours (reach one parent cell: the
// level's halo covers it, and S4.1's tiles rely on it); it is pinned like a channel, its bed
// following the start surface down to the receiving channel's bed at the junction and its
// discharge a share of the junction's by Hack's law; the density approaches the target
// regional.tributary_spacing_km (until S3.6's valley-spacing field) without exceeding it.
// Parameters: RegionalParameters::tributaries (regional_types.hpp; enabled 0: off). Planet-free.

namespace landscape
{
struct LevelContext;  // regional_pass.hpp
struct LevelState;

// Called by channelslevel after the rasterisation and before the beds are pinned into
// state.surface_m. Adds tributary cells to state.channels (channel, direction, bed_m,
// discharge_m3s, flows, globalcell) and sets state.tributaries and state.tributarycells.
// Inert in S3.0: with tributaries.enabled 0 (and until S3.3 whatever its value) nothing is
// changed and both counts are 0.
void tributarieslevel(const LevelContext& context, LevelState& state);
}
