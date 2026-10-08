#pragma once

// Upsampling of global fields onto a window's fine grid (docs/architecture/erosion-hydrology-plan.md,
// Phase E, step E2a, and section 6, "Faithful upsampling"). Every value is a pure
// function of the global field and of the fine cell's absolute position, so two
// windows compute the same value at a shared fine cell. Planet-free.

#include <vector>

#include "simulations/landscape/regional_grid.hpp"

namespace landscape
{
// Catmull-Rom bicubic interpolation of a W x H row-major global field at a point
// in global cell units, where cell (i, j)'s centre is (i + 0.5, j + 0.5): the
// column index wraps, the row index is clamped to 0 .. H - 1. Reproduces a field
// that is constant or linear in both indices exactly (away from the clamped rows).
double bicubic(const std::vector<float>& field, int width, int height, double gx, double gy);

// The field interpolated at every fine cell centre of the window's domain
// (domainfinewidth() x domainfineheight(), row-major, north row first).
std::vector<float> upsample(const std::vector<float>& field, int width, int height, const RegionalWindow& window);

// Block-mean correction, `passes` times: for every global cell of the domain the
// residual (mean of its k x k fine values minus the global value) is computed, the
// residuals are interpolated by `bicubic` at the fine centres and subtracted. A
// residual needs the cell's whole block, and the interpolation reaches two global
// cells, so a pass changes the fine values within two global cells of the domain's
// edge differently from a wider domain: the reach is 2 x passes global cells plus
// the stencil, which the halo must exceed (section 6). Rows beyond the poles do not
// exist; the residual of a domain-edge cell is computed from its block, which lies
// inside the domain by construction.
void correctblockmeans(const std::vector<float>& field, int width, int height, const RegionalWindow& window, int passes,
                       std::vector<float>& fine);
}
