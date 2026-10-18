// Sea connectivity, coastal elevation, and island cleanup.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the cell grid since T8 W2-T2
// (docs/architecture/sphere-mesh.md sections 3 and 4): loops over cell ids in the legacy order,
// neighbours from the table, boxes as walks (terrain_coasts.hpp).

#include "simulations/terrain/terrain_coasts.hpp"

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "core/parallel_rows.hpp"
#include "simulations/landscape/spherical_noise.hpp"
#include "simulations/terrain/coast_queries.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;
using namespace terrain::coastgrid;
using mesh::CellId;

// The cell-grid helpers (terrain_coasts.hpp).

vector<CellId> terrain::coastgrid::columnorder(const mesh::CellGrid& grid)
{
    vector<CellId> cells;
    cells.reserve(grid.cellcount());
    if (grid.kind() == mesh::GridKind::latlon && !grid.issubgrid())
    {
        const CellId columns = static_cast<CellId>(grid.columns());
        const CellId rows = static_cast<CellId>(grid.rows());
        for (CellId x = 0; x < columns; x++)
            for (CellId y = 0; y < rows; y++)
                cells.push_back(y * columns + x);
    }
    else if (grid.kind() == mesh::GridKind::cube && !grid.issubgrid())
    {
        const CellId n = static_cast<CellId>(grid.facecells());
        for (CellId face = 0; face < 6; face++)
            for (CellId i = 0; i < n; i++)
                for (CellId j = 0; j < n; j++)
                    cells.push_back(face * n * n + j * n + i);
    }
    else
    {
        for (CellId c = 0; c < grid.cellcount(); c++)
            cells.push_back(c);
    }
    return cells;
}

std::int32_t terrain::coastgrid::offset(const mesh::CellGrid& grid, CellId c, int dx, int dy, bool clamped)
{
    int turn = 0;
    std::int32_t at = static_cast<std::int32_t>(c);
    for (int n = 0; n < abs(dx) && at >= 0; n++)
        at = step(grid, static_cast<CellId>(at), dx < 0 ? 7 : 3, turn);
    for (int n = 0; n < abs(dy) && at >= 0; n++)
    {
        int nextturn = turn;
        const std::int32_t next = step(grid, static_cast<CellId>(at), dy < 0 ? 1 : 5, nextturn);
        if (next < 0 && clamped)
            break;
        at = next;
        turn = nextturn;
    }
    return at;
}

std::int32_t terrain::coastgrid::rawneighbour(const mesh::CellGrid& grid, CellId c, int slot)
{
    if (grid.kind() == mesh::GridKind::latlon)
    {
        const int x = grid.column(c);
        if (x == 0 && (slot == 6 || slot == 7 || slot == 8))
            return -1;
        if (x == grid.columns() - 1 && (slot == 2 || slot == 3 || slot == 4))
            return -1;
    }
    return grid.neighbour(c, slot);
}

namespace
{
// planet::seawrap on a cell: base elevation in whole metres at or below sea level, no lake.
bool seawrapcell(const planet& world, CellId c)
{
    return world.nom(c) <= world.sealevel() && world.lakesurface(c) == 0;
}

// One of the four cardinal neighbours is of the other kind (sea for land, land for sea); false
// where a cardinal neighbour is missing.
bool borders(const planet& world, CellId c, bool seaneighbour)
{
    const mesh::CellGrid& grid = world.grid();
    bool found = false;
    for (int slot = 1; slot <= mesh::slotcount; slot += 2)
    {
        const std::int32_t n = grid.neighbour(c, slot);
        if (n < 0)
            return false;
        if (seawrapcell(world, static_cast<CellId>(n)) == seaneighbour)
            found = true;
    }
    return found;
}
}

bool terrain::coastgrid::outline(const planet& world, CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    if (grid.kind() == mesh::GridKind::latlon)
        return world.outline(grid.column(c), grid.row(c));
    return world.sea(c) == 0 && borders(world, c, true);
}

bool terrain::coastgrid::coast(const planet& world, CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    if (grid.kind() == mesh::GridKind::latlon)
        return world.coast(grid.column(c), grid.row(c));
    return world.sea(c) == 1 && borders(world, c, false);
}

void terrain::coastgrid::legacyslot(const mesh::CellGrid& grid, CellId c, int& x, int& y)
{
    if (grid.kind() == mesh::GridKind::latlon)
    {
        x = grid.column(c);
        y = grid.row(c);
        return;
    }
    x = static_cast<int>(c % ARRAYWIDTH);
    y = static_cast<int>(c / ARRAYWIDTH);
}

double terrain::coastgrid::legacypixelkm(const planet& world)
{
    const double columns = static_cast<double>(std::max(1, tuning::terrain::sphere_noise_coast::legacyColumns));
    return 2.0 * 3.14159265358979323846 * world.grid().radius_m() / 1000.0 / columns;
}

vector<CellId> terrain::coastgrid::cellswithinkm(const planet& world, const mesh::Vec3& p, double radius_km)
{
    const mesh::CellGrid& grid = world.grid();
    vector<CellId> cells = mesh::cellswithin(grid, p, std::max(0.0, radius_km) * 1000.0);
    if (cells.empty())
    {
        const CellId c = grid.locate(p);
        if (c != mesh::nocell)
            cells.push_back(c);
    }
    return cells;
}

double terrain::coastgrid::distancekm(const planet& world, const mesh::Vec3& p, const mesh::Vec3& q)
{
    return mesh::greatcircle_m(p, q, world.grid().radius_m()) / 1000.0;
}

mesh::Vec3 terrain::coastgrid::movedpoint(const planet& world, CellId c, double east_km, double north_km)
{
    const mesh::CellGrid& grid = world.grid();
    if (east_km == 0.0 && north_km == 0.0)
        return grid.centre(c);
    return mesh::fromtangent(mesh::tangentframe(grid.centre(c), grid.radius_m()), east_km, north_km);
}

double terrain::coastgrid::uniformdraw()
{
    const double high = static_cast<double>(random(0, 32767));
    const double low = static_cast<double>(random(0, 32767));
    return (high * 32768.0 + low) / 1073741824.0;
}

mesh::Vec3 terrain::coastgrid::randompoint()
{
    const double u = uniformdraw();
    const double v = uniformdraw();
    return mesh::uniformpoint(u, v);
}

double terrain::coastgrid::cellchance(const planet& world, CellId c, double perpixel)
{
    if (perpixel >= 1.0)
        return 1.0;
    if (perpixel <= 0.0)
        return 0.0;
    const double pixel_m = legacypixelkm(world) * 1000.0;
    const double pixels = world.grid().area_m2(c) / (pixel_m * pixel_m);
    return 1.0 - std::pow(1.0 - perpixel, pixels);
}

vector<CellId> terrain::coastgrid::stamptemplate(const planet& world, const mesh::Vec3& p, const boolshapetemplate& shape, bool mirroreast, bool mirrorsouth)
{
    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = legacypixelkm(world);
    const mesh::TangentFrame frame = mesh::tangentframe(p, grid.radius_m());
    const int imwidth = shape.xsize();
    const int imheight = shape.ysize();
    vector<CellId> cells;
    vector<CellId> seen;
    for (int i = 0; i <= imwidth; i++)
    {
        const int pi = mirroreast ? imwidth - i : i;
        for (int j = 0; j <= imheight; j++)
        {
            const int pj = mirrorsouth ? imheight - j : j;
            if (!shape.point(pi, pj))
                continue;
            const double east = static_cast<double>(i - imwidth / 2) * pixelkm;
            const double south = static_cast<double>(j - imheight / 2) * pixelkm;
            const CellId c = grid.locate(mesh::fromtangent(frame, east, -south));
            if (c == mesh::nocell || find(seen.begin(), seen.end(), c) != seen.end())
                continue;
            seen.push_back(c);
            cells.push_back(c);
        }
    }
    return cells;
}

bool terrain::coastgrid::getoceanridge(const planet& world, CellId c, int slot)
{
    return slot >= 1 && slot <= 8 && (world.oceanridges(c) & (1 << (slot - 1))) != 0;
}

namespace
{
// Removes the bit of slot at c and the reverse bit at its neighbour; get/set read and write the code.
template <typename Get, typename Set, typename Emptied>
void deleteridgecode(const mesh::CellGrid& grid, CellId c, int slot, Get&& get, Set&& set, Emptied&& emptied)
{
    if (slot < 1 || slot > 8)
        return;
    const int code = 1 << (slot - 1);
    int current = get(c);
    if ((current & code) == 0) // If there isn't a ridge going this way
        return;
    current -= code;
    set(c, current);
    if (current == 0)
        emptied(c);

    // Now remove one going the other way.
    const std::int32_t n = grid.neighbour(c, slot);
    if (n < 0)
        return;
    const int back = 1 << (grid.reverse(c, slot) - 1);
    int other = get(static_cast<CellId>(n));
    if ((other & back) == 0)
        return;
    other -= back;
    set(static_cast<CellId>(n), other);
    if (other == 0)
        emptied(static_cast<CellId>(n));
}
}

void terrain::coastgrid::deleteridge(planet& world, CellId c, int slot)
{
    deleteridgecode(world.grid(), c, slot,
        [&](CellId cell) { return world.mountainridge(cell); },
        [&](CellId cell, int value) { world.setmountainridge(cell, value); },
        [&](CellId cell) { world.setmountainheight(cell, 0); });
}

void terrain::coastgrid::deleteoceanridge(planet& world, CellId c, int slot)
{
    deleteridgecode(world.grid(), c, slot,
        [&](CellId cell) { return world.oceanridges(cell); },
        [&](CellId cell, int value) { world.setoceanridges(cell, value); },
        [](CellId) {});
}

namespace sn = tuning::terrain::sphere_noise_coast;

terrain::coastgrid::LegacyNoise terrain::coastgrid::grain8noise()
{
    return { sn::grain8WavelengthKm, sn::grain8Persistence, sn::grain8Mean, sn::grain8Deviation };
}

terrain::coastgrid::LegacyNoise terrain::coastgrid::grain16noise()
{
    return { sn::grain16WavelengthKm, sn::grain16Persistence, sn::grain16Mean, sn::grain16Deviation };
}

terrain::coastgrid::LegacyNoise terrain::coastgrid::grain4noise()
{
    return { sn::grain4WavelengthKm, sn::grain4Persistence, sn::grain4Mean, sn::grain4Deviation };
}

terrain::coastgrid::LegacyNoise terrain::coastgrid::finenoise()
{
    return { sn::fineWavelengthKm, sn::finePersistence, sn::fineMean, sn::fineDeviation };
}

namespace
{
landscape::FractalNoise noisesettings(const planet& world, int salt, const terrain::coastgrid::LegacyNoise& legacy)
{
    landscape::FractalNoise noise;
    noise.seed = deterministiccontextseed(world.seed(), salt);
    noise.radius_km = world.grid().radius_m() / 1000.0;
    noise.wavelength_km = std::max(1.0, legacy.wavelength_km);
    noise.persistence = legacy.persistence;
    // Octaves from the first wavelength down to sphere_noise_coast.shortest_km.
    const double shortest = std::max(1.0, static_cast<double>(sn::shortestKm));
    int octaves = 1;
    while (noise.wavelength_km / std::ldexp(1.0, octaves) >= shortest && octaves < 24)
        octaves++;
    noise.octaves = octaves;
    return noise;
}
}

terrain::coastgrid::SphereNoise::SphereNoise(const planet& world, int salt, const LegacyNoise& noise)
    : legacy_(noise), field_(noisesettings(world, salt, noise))
{
    normalisation_ = field_.normalisation() > 0.0 ? field_.normalisation() : 1.0;
}

int terrain::coastgrid::SphereNoise::value(const mesh::Vec3& p, int maximum) const
{
    const double share = std::clamp(legacy_.mean + legacy_.deviation * z(p), 0.0, 1.0);
    return std::max(1, static_cast<int>(share * static_cast<double>(maximum)));
}

void terrain::coastgrid::totangent(const mesh::TangentFrame& frame, const mesh::Vec3& p, double& east_km, double& north_km)
{
    const double along = mesh::dot(p, frame.origin);
    const double scale = along > 1e-9 ? 1.0 / along : 1e9;
    east_km = mesh::dot(p, frame.east) * scale * frame.radius_m / 1000.0;
    north_km = mesh::dot(p, frame.north) * scale * frame.radius_m / 1000.0;
}

vector<double> terrain::coastgrid::spherenoise(const planet& world, int salt, const LegacyNoise& legacy)
{
    const mesh::CellGrid& grid = world.grid();
    const SphereNoise noise(world, salt, legacy);

    vector<double> values(grid.cellcount(), 0.0);
    parallelforrows(0, static_cast<int>(grid.cellcount()) - 1, [&](int first, int last)
    {
        for (int c = first; c <= last; c++)
            values[static_cast<size_t>(c)] = noise.z(grid.centre(static_cast<CellId>(c)));
    });
    return values;
}

vector<int> terrain::coastgrid::spherefractal(const planet& world, int salt, const LegacyNoise& legacy, int maximum)
{
    const vector<double> z = spherenoise(world, salt, legacy);
    vector<int> values(z.size(), 1);
    for (size_t c = 0; c < z.size(); c++)
    {
        const double share = std::clamp(legacy.mean + legacy.deviation * z[c], 0.0, 1.0);
        values[c] = std::max(1, static_cast<int>(share * static_cast<double>(maximum)));
    }
    return values;
}

namespace
{
// The legacy loops' range in rows and columns: on lat-lon the cell lies at least margin cells
// from the first and last column and row; on the cube every cell does.
bool interior(const mesh::CellGrid& grid, CellId c, int margin)
{
    if (grid.kind() != mesh::GridKind::latlon)
        return true;
    const int x = grid.column(c);
    const int y = grid.row(c);
    return x >= margin && x < grid.columns() - margin && y >= margin && y < grid.rows() - margin;
}

// Rows 1 to height - 1 of the legacy loops: the cell has its north and south neighbours.
bool hasnorthandsouth(const mesh::CellGrid& grid, CellId c)
{
    return grid.neighbour(c, 1) >= 0 && grid.neighbour(c, 5) >= 0;
}

// planet::sea at a neighbour: out of the grid (-1) reads as land, as sea(x, y) does.
bool seaat(const planet& world, std::int32_t n)
{
    return n >= 0 && world.sea(static_cast<CellId>(n));
}

// planet::surfacebelowsea at a neighbour: out of the grid reads as sea (map(x, y) is 0 there).
bool surfacebelowseaat(const planet& world, std::int32_t n)
{
    return n < 0 || world.surfacebelowsea(static_cast<CellId>(n));
}

// One step of a legacy walk that did not wrap x: on lat-lon rawneighbour, on the cube a straight step.
std::int32_t rawstep(const mesh::CellGrid& grid, CellId c, int slot, int& turn)
{
    if (grid.kind() == mesh::GridKind::latlon)
        return rawneighbour(grid, c, slot);
    return step(grid, c, slot, turn);
}

// The legacy cell-size in km: the planet's radius x pi / rows (rows: the lat-lon export's on the cube).
float cellsizekm(const planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    const int rows = grid.kind() == mesh::GridKind::latlon ? grid.rows() : 2 * grid.facecells();
    return world.radius() / 1000.0f * 3.14159265358979f / static_cast<float>(rows);
}

// The sizes of the 4-connected components of sea (sea: base at or below sea level) on cells
// with a south neighbour (the legacy flood's "y < height"), started from interior cells
// (legacy: x and y from 1 to the last but one) in descending legacy order; 0 elsewhere.
vector<int> seaareas(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    vector<int> area(grid.cellcount(), 0);
    vector<char> checked(grid.cellcount(), 0);
    vector<CellId> order = columnorder(grid);
    vector<CellId> component;

    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const CellId start = *it;
        if (!interior(grid, start, 1) || !world.basebelowsea(start) || checked[start] || area[start] != 0)
            continue;

        component.clear();
        queue<CellId> q;
        q.push(start);
        checked[start] = 1;
        while (!q.empty())
        {
            const CellId c = q.front();
            q.pop();
            if (!world.basebelowsea(c))
                continue;
            component.push_back(c);
            for (int slot = 1; slot <= mesh::slotcount; slot += 2)
            {
                const std::int32_t n = grid.neighbour(c, slot);
                if (n < 0 || grid.neighbour(static_cast<CellId>(n), 5) < 0)
                    continue;
                if (world.basebelowsea(static_cast<CellId>(n)) && !checked[n])
                {
                    checked[n] = 1;
                    q.push(static_cast<CellId>(n));
                }
            }
        }
        for (const CellId c : component)
            area[c] = static_cast<int>(component.size());
    }
    return area;
}
}

// This function removes any small seas from the map. level is level to fill them in to.

void removesmallseas(planet& world, int minseasize, int level)
{
    const mesh::CellGrid& grid = world.grid();

    // The legacy last column (from when the grid had a column more) is the mean of its neighbours.
    if (grid.kind() == mesh::GridKind::latlon)
    {
        for (int j = 0; j < grid.rows(); j++)
        {
            const CellId c = grid.cellof(grid.columns() - 1, j);
            const float amount = (world.nomfloat(static_cast<CellId>(grid.neighbour(c, 7))) + world.nomfloat(static_cast<CellId>(grid.neighbour(c, 3)))) / 2.0f;
            world.setnomfloat(c, amount);
        }
    }

    const vector<int> area = seaareas(world);

    // Now turn all the points whose area is too small into land.

    for (const CellId c : columnorder(grid))
    {
        if (world.basebelowsea(c) && area[c] < minseasize)
            world.setnomfloat(c, static_cast<float>(level));
    }
}

// This function removes straight edges on the coastlines.

namespace
{
// The cells of a straight walk from c in slot while keep(cell, turn) holds and stopafter(next
// cell, turn) does not: the visited cells, the first that fails last (-1 when the walk left the
// grid). Their count is the legacy do-while loops' advance (they stepped once more after the
// failing cell).
template <typename Keep, typename Stop>
vector<std::int32_t> straightrun(const mesh::CellGrid& grid, CellId c, int slot, bool raw, Keep&& keep, Stop&& stopafter)
{
    vector<std::int32_t> cells;
    int turn = 0;
    std::int32_t at = static_cast<std::int32_t>(c);
    while (true)
    {
        cells.push_back(at);
        const bool ok = at >= 0 && keep(static_cast<CellId>(at), turn);
        if (!ok)
            break;
        at = raw ? rawstep(grid, static_cast<CellId>(at), slot, turn) : step(grid, static_cast<CellId>(at), slot, turn);
        if (stopafter(at, turn))
            break;
    }
    return cells;
}
}

void removestraights(planet& world)
{
    const mesh::CellGrid& grid = world.grid();

    int maxlength = 4; // Any straights longer than this will be disrupted.

    int minsize = 2;
    int maxsize = 4; // Range of sizes for the chunks to be cut out of the coasts.

    auto never = [](std::int32_t, int) { return false; };
    // The legacy "y > height - 1": the next cell is the last row (no south neighbour) or beyond.
    auto lastrow = [&](std::int32_t at, int turn) { return at < 0 || grid.neighbour(static_cast<CellId>(at), turned(5, turn)) < 0; };

    auto disruptrun = [&](const vector<std::int32_t>& run, CellId c)
    {
        const int length = static_cast<int>(run.size()); // x - i of the legacy loop: one step past the last kept cell
        if (length > maxlength) // If this line is too long
        {
            const std::int32_t middle = run[static_cast<size_t>(length / 2)]; // Move to the middle of the line
            int size = random(minsize, maxsize);
            if (middle >= 0)
                terrain::coastgrid::disruptseacoastline(world, static_cast<CellId>(middle), world.mapfloat(c), 0, size);
        }
    };

    for (const CellId c : columnorder(grid)) // First do vertical/horizontal lines.
    {
        if (!world.surfacebelowsea(c))
            continue;

        if (grid.neighbour(c, 1) >= 0 && !world.surfacebelowsea(static_cast<CellId>(grid.neighbour(c, 1)))) // Sea to the south, land to the north
        {
            const vector<std::int32_t> run = straightrun(grid, c, 3, true, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(1, turn)));
            }, never);
            disruptrun(run, c);
        }

        if (grid.neighbour(c, 5) >= 0 && !world.surfacebelowsea(static_cast<CellId>(grid.neighbour(c, 5)))) // Sea to the north, land to the south
        {
            const vector<std::int32_t> run = straightrun(grid, c, 3, true, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(5, turn)));
            }, never);
            disruptrun(run, c);
        }

        if (!surfacebelowseaat(world, rawneighbour(grid, c, 3))) // Sea to the west, land to the east
        {
            const vector<std::int32_t> run = straightrun(grid, c, 5, false, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(3, turn)));
            }, lastrow);
            disruptrun(run, c);
        }

        if (!surfacebelowseaat(world, rawneighbour(grid, c, 7))) // Sea to the east, land to the west
        {
            const vector<std::int32_t> run = straightrun(grid, c, 5, false, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, rawneighbour(grid, at, turned(7, turn)));
            }, lastrow);
            disruptrun(run, c);
        }
    }

    auto nb = [&](CellId c, int slot) { return grid.neighbour(c, slot); };
    auto setfrom = [&](std::int32_t target, std::int32_t source)
    {
        if (target >= 0 && source >= 0)
            world.setnomfloat(static_cast<CellId>(target), world.nomfloat(static_cast<CellId>(source)));
    };

    for (const CellId c : columnorder(grid)) // Now do diagonals.
    {
        if (!hasnorthandsouth(grid, c) || world.surfacebelowsea(c))
            continue;

        if (terrain::coastgrid::northwestlandonly(world, c) == 1)
        {
            if (terrain::coastgrid::northwestlandonly(world, nb(c, 2)) == 1 && terrain::coastgrid::northwestlandonly(world, nb(c, 6)) == 1)
            {
                int choose = random(1, 3);

                if (choose == 1)
                    setfrom(static_cast<std::int32_t>(c), nb(c, 4));

                if (choose == 2)
                    setfrom(nb(c, 3), static_cast<std::int32_t>(c));

                if (choose == 3)
                    setfrom(nb(c, 5), static_cast<std::int32_t>(c));
            }
        }

        if (terrain::coastgrid::northeastlandonly(world, c) == 1)
        {
            if (terrain::coastgrid::northeastlandonly(world, nb(c, 4)) == 1 && terrain::coastgrid::northeastlandonly(world, nb(c, 8)) == 1)
            {
                int choose = random(1, 3);

                if (choose == 1)
                    setfrom(static_cast<std::int32_t>(c), nb(c, 6));

                if (choose == 2)
                    setfrom(nb(c, 7), static_cast<std::int32_t>(c));

                if (choose == 3)
                    setfrom(nb(c, 5), static_cast<std::int32_t>(c));
            }
        }

        if (terrain::coastgrid::southwestlandonly(world, c) == 1)
        {
            if (terrain::coastgrid::southwestlandonly(world, nb(c, 4)) == 1 && terrain::coastgrid::southwestlandonly(world, nb(c, 8)) == 1)
            {
                int choose = random(1, 3);

                if (choose == 1)
                    setfrom(static_cast<std::int32_t>(c), nb(c, 2));

                if (choose == 2)
                    setfrom(nb(c, 3), static_cast<std::int32_t>(c));

                if (choose == 3)
                    setfrom(nb(c, 1), static_cast<std::int32_t>(c));
            }
        }

        if (terrain::coastgrid::southeastlandonly(world, c) == 1)
        {
            if (terrain::coastgrid::southeastlandonly(world, nb(c, 2)) == 1 && terrain::coastgrid::southeastlandonly(world, nb(c, 6)) == 1)
            {
                int choose = random(1, 3);

                if (choose == 1)
                    setfrom(static_cast<std::int32_t>(c), nb(c, 8));

                if (choose == 2)
                    setfrom(nb(c, 7), static_cast<std::int32_t>(c));

                if (choose == 3)
                    setfrom(nb(c, 1), static_cast<std::int32_t>(c));
            }
        }
    }
}

// This function cuts a slice out of the coastline, to try to make it more interesting.

void terrain::coastgrid::disruptseacoastline(planet& world, CellId centre, float avedepth, bool raise, int size)
{
    (void)raise;
    int half = size / 2;

    box(world.grid(), centre, half, [&](CellId c, int i, int j)
    {
        if (world.nomfloat(c) > avedepth && i * i + j * j < half * half + half)
            world.setnomfloat(c, avedepth);
    });
}

void disruptseacoastline(planet& world, int centrex, int centrey, float avedepth, bool raise, int size)
{
    const mesh::CellGrid& grid = world.grid();
    const int width = world.width();
    if (centrex < 0 || centrex > width)
        centrex = wrap(centrex, width);
    const CellId centre = grid.cellof(centrex, centrey);
    if (centre != mesh::nocell)
        terrain::coastgrid::disruptseacoastline(world, centre, avedepth, raise, size);
}

// This function ensures that channels of sea are at least two pixels wide.

void widenchannels(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    const float lowered = static_cast<float>(world.sealevel() - 15);
    auto sea = [&](std::int32_t n) { return seaat(world, n); };
    auto lower = [&](std::int32_t n)
    {
        if (n >= 0)
            world.setnomfloat(static_cast<CellId>(n), lowered);
    };
    const vector<CellId> order = columnorder(grid);

    for (const CellId c : order) // First remove any one-pixel channels.
    {
        if (!interior(grid, c, 1) || !world.sea(c))
            continue;

        if (sea(grid.neighbour(c, 1)) == 0 && sea(grid.neighbour(c, 5)) == 0)
            lower(grid.neighbour(c, 1));

        if (sea(grid.neighbour(c, 7)) == 0 && sea(grid.neighbour(c, 3)) == 0)
            lower(grid.neighbour(c, 7));

        if (sea(grid.neighbour(c, 8)) == 0 && sea(grid.neighbour(c, 4)) == 0)
            lower(grid.neighbour(c, 8));

        if (sea(grid.neighbour(c, 6)) == 0 && sea(grid.neighbour(c, 2)) == 0)
            lower(grid.neighbour(c, 6));
    }

    for (const CellId c : order) // Now remove any two-pixel channels.
    {
        if (!interior(grid, c, 2) || !world.sea(c))
            continue;

        if (sea(grid.neighbour(c, 1)) == 0 && sea(grid.neighbour(c, 5)) == 1 && sea(offset(grid, c, 0, 2)) == 0)
            lower(static_cast<std::int32_t>(c));

        if (sea(grid.neighbour(c, 7)) == 0 && sea(grid.neighbour(c, 3)) == 1 && sea(offset(grid, c, 2, 0)) == 0)
            lower(grid.neighbour(c, 3));
    }
}

void loweroceans(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    int sealevel = world.sealevel();

    for (const CellId c : columnorder(grid))
    {
        if (world.sea(c) && world.nomfloat(c) > static_cast<float>(sealevel - 170))
        {
            bool nearland = 0;

            box(grid, c, 2, [&](CellId n, int, int)
            {
                if (world.sea(n) == 0)
                    nearland = 1;
                return nearland == 0;
            });

            if (nearland == 0)
                world.setnomfloat(c, static_cast<float>(sealevel - 170 + randomsign(random(0, 10))));
        }
    }
}

// This function forces land and sea at coastlines to normalise towards certain values, to improve the appearance of the coastline. Note that landheight and seadepth are relative to the sea level, not absolute.

void normalisecoasts(planet& world, int landheight, int seadepth, int severity)
{
    const mesh::CellGrid& grid = world.grid();
    const size_t count = grid.cellcount();
    int sealevel = world.sealevel();

    float depth, elev;

    const vector<CellId> order = columnorder(grid);

    for (const CellId c : order)
    {
        if (world.volcanofloat(c) == 0 && world.mountainisland(c) == 0)
        {
            if (coast(world, c))
            {
                depth = sealevel - world.nomfloat(c);
                float newdepth = (depth + static_cast<float>(seadepth * severity)) / static_cast<float>(severity + 1);

                world.setnomfloat(c, sealevel - newdepth);
            }

            if (outline(world, c))
            {
                elev = world.nomfloat(c) - sealevel;
                float newelev = (elev + static_cast<float>(landheight * severity)) / static_cast<float>(severity + 1);

                world.setnomfloat(c, sealevel + newelev);
            }
        }
    }

    // Now we need to smooth the seabeds around them.

    // The sea cells within two cells of a coastal sea cell are not shaded, as before the reach
    // of the smoothing was set in km (a presentation mark, independent of the smoothing).
    const int shadereach = 2;

    for (const CellId c : order)
    {
        if (coast(world, c) == 0)
            continue;

        box(grid, c, shadereach, [&](CellId n, int, int)
        {
            if (world.sea(n) == 1 && coast(world, n) == 0)
                world.setnoshade(n, 1);
        });
    }

    // The reach of the smoothing and its square are set in km (tuning::terrain::ocean_floor);
    // in cells they are the km over the cell size, and under half a cell there is none.
    const float cellkm = cellsizekm(world);
    const int reach = cellkm > 0.0f ? static_cast<int>(lroundf(tuning::terrain::ocean_floor::coastSmoothingRadiusKm / cellkm)) : 0;

    if (reach < 1)
        return;

    vector<char> done(count, 0);

    int amount = reach;
    int amount2 = reach;

    auto smootharound = [&](CellId c)
    {
        box(grid, c, amount, [&](CellId kk, int, int)
        {
            if (done[kk] == 0 && world.sea(kk) == 1 && coast(world, kk) == 0)
            {
                float total = 0.0f;
                int crount = 0;

                box(grid, kk, amount2, [&](CellId mm, int, int)
                {
                    if (world.sea(mm) == 1)
                    {
                        crount++;
                        total = total + world.nomfloat(mm);
                    }
                });

                if (total != 0)
                {
                    float newnom = total / static_cast<float>(crount);
                    world.setnomfloat(kk, newnom);
                    done[kk] = 1;

                    world.setnoshade(kk, 1);
                }
            }
        });
    };

    for (const CellId c : order)
    {
        if (coast(world, c) == 1)
            smootharound(c);
    }

    // Do it again, going the other way.

    fill(done.begin(), done.end(), 0);

    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        if (coast(world, *it) == 1)
            smootharound(*it);
    }
}

namespace
{
// terrain::detail::makemountainisland on a cube cell (its (x, y) form serves lat-lon): the
// same draws and stamps, the neighbour by slot and the ridge stored at both ends.
void makemountainislandcell(planet& world, CellId c, float peakheight)
{
    const mesh::CellGrid& grid = world.grid();
    int sealevel = world.sealevel();

    int iterations = random(1, 3);

    for (int n = 1; n <= iterations; n++)
    {
        int landheight = sealevel + random(5, 50);

        int dx = 0;
        int dy = 0;

        do
        {
            dx = -2 + random(1, 3);
            dy = -2 + random(1, 3);

        } while (dx == 0 && dy == 0);

        static const int slots[3][3] = { { 8, 1, 2 }, { 7, 0, 3 }, { 6, 5, 4 } }; // [dy + 1][dx + 1]
        const int slot = slots[dy + 1][dx + 1];
        const std::int32_t next = grid.neighbour(c, slot);

        if (next < 0)
            return;

        const CellId cc = static_cast<CellId>(next);

        if (world.mountainheightfloat(c) < peakheight)
            world.setmountainheight(c, peakheight);

        if (world.mountainheightfloat(cc) < peakheight)
            world.setmountainheight(cc, peakheight);

        box(grid, c, 1, [&](CellId m, int, int) { world.setmountainisland(m, 1); }); // Mark it all out as mountain islands!
        box(grid, cc, 1, [&](CellId m, int, int) { world.setmountainisland(m, 1); });

        if (world.nomfloat(c) < landheight)
            world.setnomfloat(c, static_cast<float>(landheight));

        if (world.nomfloat(cc) < landheight)
            world.setnomfloat(cc, static_cast<float>(landheight));

        // If there isn't a ridge between them, make one.

        const int code = 1 << (slot - 1);

        if ((world.mountainridge(c) & code) == 0)
        {
            world.setmountainridge(c, world.mountainridge(c) + code);

            const int back = 1 << (grid.reverse(c, slot) - 1);

            if ((world.mountainridge(cc) & back) == 0)
                world.setmountainridge(cc, world.mountainridge(cc) + back);
        }

        c = cc;

        peakheight = peakheight + randomsign(random(static_cast<int>(peakheight) / 4, static_cast<int>(peakheight) / 2));
    }
}
}

// This notes down all the very small islands on the map.

void checkislands(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    int sealevel = world.sealevel();

    for (const CellId c : columnorder(grid))
    {
        if (!hasnorthandsouth(grid, c))
            continue;

        if (world.sea(c) == 0)
        {
            int crount = -1; // Because there will definitely be one positive, i.e. the current cell.

            box(grid, c, 1, [&](CellId n, int, int)
            {
                if (world.sea(n) == 0)
                    crount++;
            });

            if (crount == 0) // This is a one-tile island!
            {
                world.setisland(c, 1);

                box(grid, c, 1, [&](CellId n, int, int)
                {
                    if (world.sea(n) == 1)
                        world.setnoshade(n, 1);
                });
            }
        }
        else
        {
            float volcano = world.volcanofloat(c);

            if (volcano != 0)
            {
                bool extinct = 0;

                if (volcano < 0)
                {
                    extinct = 1;
                    volcano = 0 - volcano;
                }

                const float totalheight = world.nomfloat(c) + volcano;

                if (wholemetres(totalheight) > sealevel)
                {
                    volcano = totalheight - static_cast<float>(sealevel);

                    if (grid.kind() == mesh::GridKind::latlon)
                        terrain::detail::makemountainisland(world, grid.column(c), grid.row(c), volcano);
                    else
                        makemountainislandcell(world, c, volcano);

                    if (extinct == 1)
                        volcano = 0 - volcano;

                    world.setvolcano(c, volcano);
                }
            }
        }
    }
}

// This extends the noshade areas by one cell in each direction.

void extendnoshade(planet& world)
{
    const mesh::CellGrid& grid = world.grid();

    vector<char> extra(grid.cellcount(), 0); // This is to mark all cells touching the no-shade ones as themselves no-shade.

    for (const CellId c : columnorder(grid))
    {
        if (world.noshade(c) == 1)
        {
            box(grid, c, 1, [&](CellId n, int, int)
            {
                if (world.sea(n) == 1)
                    extra[n] = 1;
            });
        }
    }

    for (CellId c = 0; c < grid.cellcount(); c++)
    {
        if (extra[c] == 1)
            world.setnoshade(c, 1);
    }
}
