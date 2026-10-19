#pragma once

// The glacial fields of the regional pass on the fine grid and the ice catchment of the
// global network (docs/architecture/erosion-hydrology-plan.md, Phase E5, step E5.0, the
// contracts; section 6, 'Decided for Phase E5'). Shared by the three glacial forms, each in
// its own translation unit: the troughs (regional_troughs.hpp, E5.1), the cirques
// (regional_cirques.hpp, E5.2) and the ice above the ELA (regional_ice.hpp, E5.3). Every
// value is a pure function of the global fields, the seed and absolute position. Planet-free.
//
// The ELA. cryosphere.equilibrium_line_altitude is the erosion run's ELA, computed from
// the coarse climate cooled by landscape.glacial.cooling_c (glacial.hpp): the ELA of the
// glacial past that shaped the troughs and cirques. The present climate's ELA lies
// cooling_c / lapse rate higher (RegionalGlobalInputs::presentelaoffset_m, 881 m on
// erosion run 9: 4 C over its 4.54 C/km): the one that decides the present-day ice. Both are the global field
// upsampled at the unwarped fine centres (a smooth climate quantity, as the temperatures;
// no lapse correction: an altitude, not a temperature).
//
// Nothing here is read at the defaults: with every glacial feature off (glacialformson
// false) regionalglacialfields returns an inactive struct with empty vectors without
// touching a field, so the Phase R2 results stay bit-identical by construction.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
// Domain fields (domainfinewidth() x domainfineheight(), row-major) when `active`; empty
// otherwise.
struct FineGlacial
{
    bool active = false;              // the ELA field is present (W x H) and a glacial feature is on
    std::vector<float> runela_m;      // the run's ELA upsampled (m above datum)
    std::vector<float> presentela_m;  // runela_m + presentelaoffset_m, per cell
    std::vector<float> glaciation;    // cryosphere.glaciation upsampled, clamped to 0..1 (0 everywhere when the
                                      // global field is empty)
};

// True when a glacial form or the ELA ice is on: troughs.width_km > 0, cirques.radius_km > 0
// or ice.thicknessfraction > 0.
bool glacialformson(const RegionalParameters& parameters);

// The fine fields above; inactive (empty) unless glacialformson and inputs.equilibriumline_m
// holds width x height values.
FineGlacial regionalglacialfields(const RegionalGlobalInputs& inputs, const RegionalWindow& window,
                                  const RegionalParameters& parameters);
// The cube form (T8 W2-W2d; lat-lon windows take the form above): one value per cell of a cube
// window refining inputs.grid, in its local order; the fields gathered by G (upsampling.hpp).
FineGlacial regionalglacialfields(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window,
                                  const RegionalParameters& parameters);

// The ice catchment per global cell (W x H, km2): the sum over the cell and every cell
// upstream of it on the global network of glaciation x the cell's area (SphereMetric's
// row area), so a cell's value is the glaciated area that drained through it during the
// run; 0 on sea cells and wherever the glaciation field is empty. Receivers as
// regional_valleys.cpp reads them (a direction beyond a pole: none); a cell on a receiver
// loop (bad input) takes its own term. A pure function of the global fields, so every
// window agrees. Empty when the glaciation field is not W x H.
std::vector<double> glacialcatchment_km2(const RegionalGlobalInputs& inputs);

struct LevelContext;  // regional_pass.hpp (the hook interface, multiscale plan S3.0)
struct LevelState;

// LevelStep::glacialforms: state.glacial = regionalglacialfields, then with details.troughs
// carvetroughs and with details.cirques carvecirques (regional_troughs.hpp, regional_cirques.hpp)
// on state.elevation_m with the level's band (formband, regional_valleys.hpp: each form at the
// level where its half-width or radius in level cells lies in the band, level 1 also above it),
// each followed by channelcellmeans with channelcellcap_m (multiscale plan S3.5). The loop's
// diagnostics count the kept cells each changed (LevelDiagnostics::troughcells, cirquecells).
void glacialformslevel(const LevelContext& context, LevelState& state);
}
