// Elevation smoothing, relief noise, and depression filling.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the cell grid
// (docs/architecture/sphere-mesh.md sections 3, 11 and 12) since T8 W2-T1.

#include "simulations/terrain/terrain_internal.hpp"
#include "simulations/terrain/terrain_elevation.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "core/parallel_rows.hpp"
#include "simulations/landscape/spherical_noise.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;

namespace
{
bool checkdepressiontile(planet& world, std::vector<float>& filledmap, mesh::CellId c, int e, bool somethingdone, const std::vector<int>& noise);

namespace sn = tuning::terrain::sphere_noise;

// Salts of the sphere fields that replace fractals drawn from a stage's stream.
constexpr int smoothingsalt = 0x5301;
constexpr int canyonbasesalt = 0x5302;
constexpr int canyonmasksalt = 0x5303;

landscape::FractalNoiseField makefield(const planet& world, std::uint64_t seed, double wavelength_km, double persistence)
{
    landscape::FractalNoise noise;
    noise.seed = seed;
    noise.radius_km = world.grid().radius_m() / 1000.0;
    noise.wavelength_km = wavelength_km;
    noise.persistence = persistence;
    // Octaves from the first wavelength down to terrain.sphere_noise.shortest_km.
    const double shortest = max(1.0, static_cast<double>(sn::shortestKm));
    int octaves = 1;
    while (wavelength_km / pow(2.0, octaves) >= shortest && octaves < 24)
        octaves++;
    noise.octaves = octaves;
    return landscape::FractalNoiseField(noise);
}

// The legacy share of the maximum at a point: clamp(mean + deviation x z, 0, 1).
double legacyshare(const landscape::FractalNoiseField& field, const terrain::SphereFractal& fractal, const mesh::Vec3& p)
{
    const double z = field.sum(p.x, p.y, p.z) / field.normalisation();
    return clamp(fractal.mean + fractal.deviation * z, 0.0, 1.0);
}

float legacyvalue(const landscape::FractalNoiseField& field, const terrain::SphereFractal& fractal, const mesh::Vec3& p,
                  float maximum)
{
    return max(1.0f, static_cast<float>(static_cast<double>(maximum) * legacyshare(field, fractal, p)));
}

terrain::SphereFractal warpfractal()
{
    return { sn::warpWavelengthKm, sn::warpPersistence, sn::warpMean, sn::warpDeviation };
}

terrain::SphereFractal warpvaryfractal()
{
    return { sn::warpVaryWavelengthKm, sn::warpVaryPersistence, sn::warpVaryMean, sn::warpVaryDeviation };
}

template <typename Fn>
void forcells(const planet& world, Fn&& fn)
{
    parallelforrows(0, static_cast<int>(world.cellcount()) - 1, [&](int first, int last)
    {
        for (mesh::CellId c = static_cast<mesh::CellId>(first); c <= static_cast<mesh::CellId>(last); c++)
            fn(c);
    }, 256);
}

// Averages a float field over each land cell's disc of radius legacy pixels (in km), in the
// legacy column-major order and in place, as planet::smoothoverland did over its square box.
void smoothoverland(planet& world, int amount, float (planet::*get)(mesh::CellId) const,
                    void (planet::*set)(mesh::CellId, float))
{
    const double radius_km = (static_cast<double>(amount) + 0.5) * terrain::legacycellkm(world);
    const float maximum = static_cast<float>(world.maxelevation());
    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (world.sea(c) != 0)
            continue;

        const vector<mesh::CellId> cells = terrain::cellswithinkm(world, c, radius_km);
        float ave = 0.0f;
        for (const mesh::CellId n : cells)
            ave += (world.*get)(n);
        ave /= static_cast<float>(cells.size());

        if (ave > 0.0f && ave < maximum)
            (world.*set)(c, ave);
    }
}
}

terrain::SphereFractal terrain::terrainfractal()
{
    return { sn::terrainFractalWavelengthKm, sn::terrainFractalPersistence, sn::terrainFractalMean, sn::terrainFractalDeviation };
}

terrain::SphereFractal terrain::smoothingfractal()
{
    return { sn::smoothingWavelengthKm, sn::smoothingPersistence, sn::smoothingMean, sn::smoothingDeviation };
}

terrain::SphereFractal terrain::canyonbasefractal()
{
    return { sn::canyonBaseWavelengthKm, sn::canyonBasePersistence, sn::canyonBaseMean, sn::canyonBaseDeviation };
}

terrain::SphereFractal terrain::canyonmaskfractal()
{
    return { sn::canyonMaskWavelengthKm, sn::canyonMaskPersistence, sn::canyonMaskMean, sn::canyonMaskDeviation };
}

terrain::SphereFractal terrain::archipelagofractal()
{
    return { sn::archipelagoWavelengthKm, sn::archipelagoPersistence, sn::archipelagoMean, sn::archipelagoDeviation };
}

terrain::SphereFractal terrain::craterdensityfractal()
{
    return { sn::craterDensityWavelengthKm, sn::craterDensityPersistence, sn::craterDensityMean, sn::craterDensityDeviation };
}

double terrain::legacycellkm(const planet& world)
{
    const double pi = 3.14159265358979323846;
    return 2.0 * pi * world.grid().radius_m() / 1000.0 / static_cast<double>(max(1, sn::legacyColumns));
}

double terrain::legacycellcount()
{
    const double columns = static_cast<double>(max(1, sn::legacyColumns));
    return columns * columns / 2.0;
}

vector<float> terrain::spherefractal(const planet& world, std::uint64_t seed, const SphereFractal& fractal, float maximum)
{
    const landscape::FractalNoiseField field = makefield(world, seed, fractal.wavelength_km, fractal.persistence);
    const mesh::CellGrid& grid = world.grid();
    vector<float> values(world.cellcount(), 0.0f);
    forcells(world, [&](mesh::CellId c) { values[c] = legacyvalue(field, fractal, grid.centre(c), maximum); });
    return values;
}

vector<float> terrain::spherefractalwarped(const planet& world, std::uint64_t seed, const SphereFractal& fractal,
                                           float maximum, int warpfactor)
{
    const landscape::FractalNoiseField field = makefield(world, seed, fractal.wavelength_km, fractal.persistence);
    const SphereFractal warp = warpfractal();
    const SphereFractal vary = warpvaryfractal();
    const landscape::FractalNoiseField east = makefield(world, deterministiccombine(seed, 1), warp.wavelength_km, warp.persistence);
    const landscape::FractalNoiseField north = makefield(world, deterministiccombine(seed, 2), warp.wavelength_km, warp.persistence);
    const landscape::FractalNoiseField varying = makefield(world, deterministiccombine(seed, 3), vary.wavelength_km, vary.persistence);
    const mesh::CellGrid& grid = world.grid();
    const double pixel_km = legacycellkm(world);
    vector<float> values(world.cellcount(), 0.0f);
    forcells(world, [&](mesh::CellId c)
    {
        const mesh::Vec3 p = grid.centre(c);
        const double reach = static_cast<double>(warpfactor) * legacyshare(varying, vary, p) * pixel_km;
        const double east_km = (legacyshare(east, warp, p) - warp.mean) * reach;
        const double north_km = (legacyshare(north, warp, p) - warp.mean) * reach;
        const mesh::Vec3 q = mesh::fromtangent(mesh::tangentframe(p, grid.radius_m()), east_km, north_km);
        values[c] = legacyvalue(field, fractal, q, maximum);
    });
    return values;
}

vector<mesh::CellId> terrain::cellswithinkm(const planet& world, mesh::CellId c, double radius_km)
{
    const mesh::CellGrid& grid = world.grid();
    vector<mesh::CellId> cells = mesh::cellswithin(grid, grid.centre(c), max(0.0, radius_km) * 1000.0);
    if (cells.empty())
        cells.push_back(c);
    return cells;
}

double terrain::distancekm(const planet& world, mesh::CellId a, mesh::CellId b)
{
    const mesh::CellGrid& grid = world.grid();
    return mesh::greatcircle_m(grid.centre(a), grid.centre(b), grid.radius_m()) / 1000.0;
}

mesh::CellId terrain::movedcell(const planet& world, mesh::CellId c, double east_km, double north_km)
{
    if (east_km == 0.0 && north_km == 0.0)
        return c;
    const mesh::CellGrid& grid = world.grid();
    const mesh::CellId moved = grid.locate(mesh::fromtangent(mesh::tangentframe(grid.centre(c), grid.radius_m()), east_km, north_km));
    return moved == mesh::nocell ? c : moved;
}

double terrain::legacyuniform()
{
    const double high = static_cast<double>(random(0, 32767));
    const double low = static_cast<double>(random(0, 32767));
    return (high * 32768.0 + low) / 1073741824.0;
}

mesh::CellId terrain::randomcellbyarea(const planet& world)
{
    return mesh::randomcellbyarea(world.grid(), [] { return legacyuniform(); });
}

vector<float> terrain::spherefractalwrapped(const planet& world, std::uint64_t seed, const SphereFractal& fractal, float maximum)
{
    const landscape::FractalNoiseField field = makefield(world, seed, fractal.wavelength_km, fractal.persistence);
    const mesh::CellGrid& grid = world.grid();
    vector<float> values(world.cellcount(), 0.0f);
    forcells(world, [&](mesh::CellId c)
    {
        const mesh::Vec3 p = grid.centre(c);
        const double share = fractal.mean + fractal.deviation * field.sum(p.x, p.y, p.z) / field.normalisation();
        values[c] = static_cast<float>(static_cast<double>(maximum) * (share - floor(share)));
    });
    return values;
}

// This function smoothes, but without turning any land to sea or vice versa: each sea cell is
// averaged over the cells within a radius (km, tuning::terrain::ocean_floor) chosen between the
// two radii by a sphere noise across its own range; a radius that reaches no other cell centre
// leaves the cell. Land is never written.

void smoothland(planet& world, int amount)
{
    (void)amount;
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();

    const vector<float> fractal = terrain::spherefractal(world, deterministiccontextseed(world.seed(), smoothingsalt),
                                                         terrain::smoothingfractal(), maxelev);

    float fractalmin = fractal.empty() ? 0.0f : fractal[0];
    float fractalmax = fractalmin;

    for (const float value : fractal)
    {
        fractalmin = min(fractalmin, value);
        fractalmax = max(fractalmax, value);
    }

    const float minimumradiuskm = tuning::terrain::ocean_floor::smoothingMinimumRadiusKm;
    const float maximumradiuskm = tuning::terrain::ocean_floor::smoothingMaximumRadiusKm;
    const float fractalrange = max(1.0f, fractalmax - fractalmin);

    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (world.sea(c) == 0) // Land is never written.
            continue;

        const float position = (fractal[c] - fractalmin) / fractalrange;
        const float radiuskm = minimumradiuskm + position * (maximumradiuskm - minimumradiuskm);
        const vector<mesh::CellId> cells = terrain::cellswithinkm(world, c, radiuskm);

        if (cells.size() < 2) // The radius reaches no other cell: the cell keeps its value.
            continue;

        float ave = 0.0f;

        for (const mesh::CellId n : cells)
            ave = ave + world.nomfloat(n);

        ave = ave / static_cast<float>(cells.size());

        if (ave > 0 && ave < maxelev && wholemetres(ave) <= sealevel) // The new value must stay sea.
            world.setnomfloat(c, ave);
    }
}

// This function smoothes only the sea.

void smoothonlysea(planet& world, int amount)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    float maxelev = (float)world.maxelevation();

    for (int i = 0; i <= width; i++)
    {
        for (int j = 0; j <= height; j++)
        {
            float crount = 0.0f;
            float ave = 0.0f;
            int origland = 0;
            int thisamount = amount;

            if (world.basebelowsea(i, j))
            {
                for (int k = i - thisamount; k <= i + thisamount; k++)
                {
                    int kk = k;

                    if (kk<0 || kk>width)
                        kk = wrap(kk, width);

                    for (int l = j - thisamount; l <= j + thisamount; l++)
                    {
                        if (l > 0 && l < height && world.basebelowsea(kk, l))
                        {
                            ave = ave + world.nomfloat(kk, l);
                            crount++;
                        }
                    }
                }

                if (crount > 0.0f)
                {
                    ave = ave / crount;

                    if (ave > 0.0f && ave < maxelev)
                        world.setnomfloat(i, j, ave);
                }
            }
        }
    }
}

// This function creates areas of raised elevation around where canyons might form later.

void createextraelev(planet& world)
{
    float gravity = world.gravity();

    float maxextra = 1000.0f; // Maximum amount that the extraelev array can have.

    float maxelev = (float)world.maxelevation();

    float div = maxelev / (maxextra * 2.0f);

    // First create a fractal for this map (the sphere noise of the legacy grain-16 fractal, warped).

    int warpfactor = random(20, 80);
    const vector<float> base = terrain::spherefractalwarped(world, deterministiccontextseed(world.seed(), canyonbasesalt),
                                                            terrain::canyonbasefractal(), maxelev, warpfactor);

    forcells(world, [&](mesh::CellId c) { world.setextraelev(c, base[c]); });

    warpfactor = random(20, 80);
    const vector<float> mask = terrain::spherefractalwarped(world, deterministiccontextseed(world.seed(), canyonmasksalt),
                                                            terrain::canyonmaskfractal(), maxelev, warpfactor);

    // Now alter it a bit.

    forcells(world, [&](mesh::CellId c)
    {
        if (world.sea(c) == 1)
            world.setextraelev(c, 0);
        else
        {
            float e = world.extraelevfloat(c);
            e = e / div;
            e = e - maxextra;

            if (e < 0.0f)
                e = 0.0f;

            if (e > maxextra)
                e = maxextra;

            e = e * gravity;

            world.setextraelev(c, e);
        }
    });

    // That gives us a *lot* of extra elevation across the map. We want to apply it a bit more judiciously.
    // So we use a second fractal to mask the first, so that fewer areas get extra elevation.

    float multfactor = (float)random(1, 10);
    multfactor = multfactor / 10.0f + 1.0f; // Different worlds will have different amounts of extra elevation.

    forcells(world, [&](mesh::CellId c)
    {
        if (world.sea(c) == 0)
        {
            float e = mask[c];
            e = e / div;
            e = e - maxextra;

            e = e * multfactor; // The higher this is, the more widespread the extra elevation will be.

            if (e < 0.0f)
                e = 0.0f;

            if (e < world.extraelevfloat(c))
                world.setextraelev(c, e);
        }
    });

    // Now smooth it over discs of 2 and 5 legacy pixels.

    smoothoverland(world, 2, &planet::extraelevfloat, &planet::setextraelev);
    smoothoverland(world, 5, &planet::extraelevfloat, &planet::setextraelev);
}

namespace
{
// The cells of a grid as panels of columns x rows, id = panel x columns x rows + row x columns
// + column: one panel on lat-lon (W x W/2), the six faces on the cube (N x N).
struct ScanPanels
{
    int panels = 0;
    int columns = 0;
    int rows = 0;
};

ScanPanels scanpanels(const mesh::CellGrid& grid)
{
    if (grid.kind() == mesh::GridKind::cube)
        return { 6, grid.facecells(), grid.facecells() };
    return { 1, grid.columns(), grid.rows() };
}

// One of the eight directional scans of the fill: column by column (columnmajor) or row by row,
// columns and rows ascending or descending; the panels in turn, descending when both are.
// On lat-lon these are the legacy (x, y) scans; skipfirstcolumn keeps the seventh scan's
// "x > 0" (column 0 is left to the other seven).
vector<mesh::CellId> scanorder(const ScanPanels& panels, bool columnmajor, bool columnsdescending, bool rowsdescending,
                               bool skipfirstcolumn = false)
{
    vector<mesh::CellId> cells;
    cells.reserve(static_cast<size_t>(panels.panels) * panels.columns * panels.rows);
    const int firstcolumn = skipfirstcolumn ? 1 : 0;
    for (int p = 0; p < panels.panels; p++)
    {
        const int panel = columnsdescending && rowsdescending ? panels.panels - 1 - p : p;
        const mesh::CellId base = static_cast<mesh::CellId>(panel) * panels.columns * panels.rows;
        const int outercount = columnmajor ? panels.columns : panels.rows;
        const int innercount = columnmajor ? panels.rows : panels.columns;
        const bool outerdescending = columnmajor ? columnsdescending : rowsdescending;
        const bool innerdescending = columnmajor ? rowsdescending : columnsdescending;
        for (int o = 0; o < outercount; o++)
        {
            const int outer = outerdescending ? outercount - 1 - o : o;
            for (int n = 0; n < innercount; n++)
            {
                const int inner = innerdescending ? innercount - 1 - n : n;
                const int column = columnmajor ? outer : inner;
                const int row = columnmajor ? inner : outer;
                if (column < firstcolumn)
                    continue;
                cells.push_back(base + static_cast<mesh::CellId>(row * panels.columns + column));
            }
        }
    }
    return cells;
}

// A land cell with a sea cell among its four cardinal neighbours (planet::outline on lat-lon:
// a cell missing a cardinal neighbour, a polar row there, is not).
bool coastalcell(const planet& world, mesh::CellId c)
{
    if (world.sea(c) != 0)
        return false;
    const mesh::CellGrid& grid = world.grid();
    bool coastal = false;
    for (int slot = 1; slot <= mesh::slotcount; slot += 2)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);
        if (neighbour < 0)
            return false;
        if (world.sea(static_cast<mesh::CellId>(neighbour)) == 1)
            coastal = true;
    }
    return coastal;
}
}

// Fill land depressions using Planchon-Darboux directional scans. The noisy
// minimum gradient prevents flat areas from blocking later river routing.

void depressionfill(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = world.cellcount();
    int maxelev = world.maxelevation();

    int e = 2; // This is the extra amount added to ensure that everywhere slopes.

    vector<int> noise(count, 0); // This will contain noise that allows us to vary the value of e from tile to tile.

    for (const mesh::CellId c : columnmajorcells(grid))
        noise[c] = random(0, 5);

    vector<float> filledmap(count, 0.0f); // This will be the new version of the map.

    for (mesh::CellId c = 0; c < count; c++) // First, fill the new map to a huge height.
    {
        if (world.basebelowsea(c))
            filledmap[c] = world.nomfloat(c); // Sea tiles start with the normal heights.
        else
        {
            if (coastalcell(world, c)) // If this is a coastal tile
                filledmap[c] = world.nomfloat(c); // Coastal tiles start with the normal heights.
            else
                filledmap[c] = static_cast<float>(maxelev * 2); // Other tiles start very high
        }
    }

    // Eight scans, because there are eight possible ways to scan the map and we must rotate between them.
    const ScanPanels panels = scanpanels(grid);
    const vector<mesh::CellId> scans[8] =
    {
        scanorder(panels, true, false, false),
        scanorder(panels, false, true, true),
        scanorder(panels, true, true, false),
        scanorder(panels, false, true, false),
        scanorder(panels, true, false, true),
        scanorder(panels, false, false, true),
        scanorder(panels, true, true, true, true),
        scanorder(panels, false, false, false),
    };

    // Now we're ready to start!

    bool somethingdone = 0;

    do
    {
        for (int scan = 0; scan < 8; scan++)
        {
            if (scan > 0)
            {
                if (somethingdone == 0)
                    break;
                somethingdone = 0;
            }

            for (const mesh::CellId c : scans[scan])
            {
                if (filledmap[c] > world.nomfloat(c))
                    somethingdone = checkdepressiontile(world, filledmap, c, e, somethingdone, noise);
            }
        }
    } while (somethingdone == 1);

    for (mesh::CellId c = 0; c < count; c++)
        world.setnomfloat(c, filledmap[c]);
}

namespace
{
// This does the actual checking/filling of each tile: its eight neighbours from a random
// one, a cardinal step one metre steeper than a diagonal.

bool checkdepressiontile(planet& world, vector<float>& filledmap, mesh::CellId c, int e, bool somethingdone, const vector<int>& noise)
{
    const mesh::CellGrid& grid = world.grid();

    int start = random(0, 7);

    for (int n = start; n <= start + 7; n++)
    {
        int nn = wrap(n, 7);

        const std::int32_t neighbour = grid.neighbour(c, nn + 1);

        if (neighbour >= 0)
        {
            const mesh::CellId k = static_cast<mesh::CellId>(neighbour);
            int ee = 0;

            if ((nn + 1) % 2 == 1)
                ee = e + 1 + noise[c];

            else
                ee = e + noise[c];

            if (world.nomfloat(c) >= filledmap[k] + ee)
            {
                filledmap[c] = world.nomfloat(c);
                somethingdone = 1;

            }
            else
            {
                if (filledmap[c] > filledmap[k] + ee)
                {
                    filledmap[c] = filledmap[k] + ee;
                    somethingdone = 1;
                }
            }
        }
    }

    return (somethingdone);
}
}

// This function adds a bit of noise to the land.

void addlandnoise(planet& world)
{
    int sealevel = world.sealevel();
    int maxelev = world.maxelevation();

    int maxadjust = 10;
    int adjustchance = 5; // The lower this is, the more there will be.

    for (const mesh::CellId c : columnmajorcells(world.grid()))
    {
        if (!world.basebelowsea(c) && random(1, adjustchance) == 1)
        {
            int adjust = randomsign(random(1, maxadjust));

            float newamount = world.nomfloat(c) + static_cast<float>(adjust);

            if (newamount <= static_cast<float>(sealevel))
                newamount = static_cast<float>(sealevel);

            if (newamount > static_cast<float>(maxelev))
                newamount = static_cast<float>(maxelev);

            world.setnomfloat(c, newamount);
        }
    }
}
