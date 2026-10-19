#pragma once

// Present-day ice above the equilibrium line on the fine grid (docs/architecture/
// erosion-hydrology-plan.md, Phase E5: the E5.0 contract, which step E5.3 implements;
// section 6, 'Decided for Phase E5'). Called by regional_pass.cpp after the lakes and the
// present-day ice of regional_water.hpp, whose ice it adds to. Planet-free.
//
// Contract (E5.3). No-op while ice.thicknessfraction <= 0 or the glacial fields are
// inactive. On the domain after the lakes: a fine cell that is not sea whose surface z
// (the final elevation) satisfies z >= E + ice.edgenoise_m x n, with E the present
// climate's ELA at the cell (glacial.presentela_m: the run's ELA raised by
// presentelaoffset_m, cooling_c over the lapse rate) and n a RegionalNoise of its own salt
// (persistence coast.noisepersistence, standard deviation about 1, so the edge wanders
// by about edgenoise_m and at most a few times it; no noise with edgenoise_m 0) so
// the ice edge does not follow the interpolated ELA's contour exactly, carries ice of
// thickness max(the thickness it has, max(1 m, min(ice.thicknessmax_m,
// ice.thicknessfraction x (z - E)))): the accumulation areas of the present climate that
// the coarse grid missed (a ridge above the ELA inside a global cell without ice). Ice
// that regional_water.hpp placed (the global glaciers upsampled, tongues below the ELA
// included) is never removed or thinned; a lake cell above the ELA is iced like land (a
// frozen lake). Every value is a pure function of the fields, the seed and absolute
// position; the rule reads one cell, so the halo rule holds trivially. At k 1 the noise
// is inactive (RegionalNoise, as the coast noise) and the threshold is the plain ELA.
// E5.3 keeps the overlap tests bit-identical with the ELA ice on and adds its variants,
// updates the regional bundle's description of cryosphere.ice_thickness
// (src/io/regional_bundle.cpp) and measured the iced share of the land against the
// present-day ice of the global field on erosion run 9's window B at k 16 (the close's
// registered runs show A and C): the global ice is 2.6 m thick per metre above the
// present ELA there (ice-sheet thicknesses at the 3,000 m cap), unusable as a fraction,
// so thicknessfraction's default 0.3 is chosen for a visible thickness (median 62 m on
// B) under the cap.

#include <cstdint>
#include <vector>

#include "simulations/landscape/regional_glacial.hpp"
#include "simulations/landscape/regional_grid.hpp"
#include "simulations/landscape/regional_types.hpp"

namespace landscape
{
void iceabovetheela(const RegionalGlobalInputs& inputs, const RegionalWindow& window, const RegionalParameters& parameters,
                    const FineGlacial& glacial, const std::vector<uint8_t>& sea, const std::vector<float>& elevation_m,
                    std::vector<float>& icethickness_m);
// The cube form (T8 W2-W2d; lat-lon windows take the form above): one value per cell of a cube
// window refining inputs.grid, in its local order; the edge noise at the cell's centre.
void iceabovetheela(const RegionalGlobalInputs& inputs, const mesh::CellGrid& window, const RegionalParameters& parameters,
                    const FineGlacial& glacial, const std::vector<uint8_t>& sea, const std::vector<float>& elevation_m,
                    std::vector<float>& icethickness_m);
}
