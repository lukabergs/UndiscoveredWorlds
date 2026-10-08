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

// A domain warp (erosion plan Phase R, R0 contract; R2 fills it): the displacement of
// every domain fine cell's sampling position, in global cell units (dx eastwards, dy
// southwards, as gx and gy grow), a pure function of the fine cell's absolute
// position and the seed (regional_surface.hpp, regionalwarp). Empty vectors: no warp.
struct FineWarp
{
    std::vector<double> dx;
    std::vector<double> dy;
    bool empty() const { return dx.empty(); }
};

// The field interpolated by `bicubic` at every domain fine cell centre displaced by
// the warp; an empty warp gives `upsample` exactly. Throws std::invalid_argument when
// a non-empty warp does not hold one displacement per domain fine cell.
std::vector<float> upsample(const std::vector<float>& field, int width, int height, const RegionalWindow& window,
                            const FineWarp& warp);

// Block-mean correction, `passes` times: for every global cell of the domain the
// residual (mean of its k x k fine values minus the global value) is computed, the
// residuals are interpolated by `bicubic` at the fine centres and subtracted. The
// block means approach the global values without reaching them: with +-10 m of
// per-fine-cell noise on a smooth field at k = 4 the largest interior error fell
// from 8.5 m to 0.44 m after two passes (E2a, Landscape.UpsamplingBlockMeans*).
// The stencil reaches two global cells beyond the domain's edge, where no residual
// exists: those cells take the nearest domain cell's residual (a clamp, as rows are
// clamped at the poles). So a pass changes the fine values within two global cells
// of a non-polar domain edge differently from a wider domain, and the reach is
// 2 x passes global cells (upsample itself reads the global field and adds none),
// which the halo must exceed (section 6). The residual of a domain-edge cell is
// computed from its block, which lies inside the domain by construction.
void correctblockmeans(const std::vector<float>& field, int width, int height, const RegionalWindow& window, int passes,
                       std::vector<float>& fine);
}
