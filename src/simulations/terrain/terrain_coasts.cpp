// Sea connectivity, coastal elevation, and island cleanup.
// Algorithms by Jonathan Hill, extracted from globalterrain.cpp; on the cell grid since T8 W2-T2
// (docs/architecture/sphere-mesh.md sections 3, 4 and 11): loops over cell ids in the legacy
// order, neighbours from the table; neighbour tests in cells, distances in legacy pixels of km.

#include "simulations/terrain/terrain_coasts.hpp"

#include "simulations/terrain/terrain_internal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <queue>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "core/grid_shapes.hpp"
#include "simulations/terrain/coast_queries.hpp"
#include "simulations/terrain/terrain_elevation.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "simulations/terrain/terrain_volcanoes.hpp"
#include "world/functions.hpp"
#include "world/planet.hpp"

using namespace std;
using namespace terrain::coastgrid;
using mesh::CellId;

// The cell-grid helpers (terrain_coasts.hpp).

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
    return world.sea(c) == 0 && borders(world, c, true);
}

bool terrain::coastgrid::coast(const planet& world, CellId c)
{
    return world.sea(c) == 1 && borders(world, c, false);
}

bool terrain::coastgrid::seaat(const planet& world, std::int32_t n)
{
    return n >= 0 && world.sea(static_cast<CellId>(n));
}

terrain::SphereFractal terrain::coastgrid::finefractal()
{
    namespace sc = tuning::terrain::sphere_noise_coast;
    return { sc::fineWavelengthKm, sc::finePersistence, sc::fineMean, sc::fineDeviation };
}

terrain::SphereFractal terrain::coastgrid::trenchfractal()
{
    namespace sn = tuning::terrain::sphere_noise;
    return { sn::warpVaryWavelengthKm, sn::warpVaryPersistence, sn::warpVaryMean, sn::warpVaryDeviation };
}

terrain::SphereFractal terrain::coastgrid::ridgeanglefractal()
{
    namespace sn = tuning::terrain::sphere_noise;
    return { sn::warpVaryWavelengthKm, sn::warpVaryPersistence, 0.5, tuning::terrain::sphere_noise_coast::ridgeAngleDeviation };
}

double terrain::coastgrid::cellchance(const planet& world, CellId c, double perpixel)
{
    if (perpixel >= 1.0)
        return 1.0;
    if (perpixel <= 0.0)
        return 0.0;
    const double pixel_m = terrain::legacycellkm(world) * 1000.0;
    const double pixels = world.grid().area_m2(c) / (pixel_m * pixel_m);
    return 1.0 - std::pow(1.0 - perpixel, pixels);
}

vector<CellId> terrain::coastgrid::stamptemplate(const planet& world, const mesh::Vec3& p, const boolshapetemplate& shape, bool mirroreast, bool mirrorsouth)
{
    const mesh::CellGrid& grid = world.grid();
    const double pixelkm = terrain::legacycellkm(world);
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

void terrain::coastgrid::totangent(const mesh::TangentFrame& frame, const mesh::Vec3& p, double& east_km, double& north_km)
{
    const double along = mesh::dot(p, frame.origin);
    const double scale = along > 1e-9 ? 1.0 / along : 1e9;
    east_km = mesh::dot(p, frame.east) * scale * frame.radius_m / 1000.0;
    north_km = mesh::dot(p, frame.north) * scale * frame.radius_m / 1000.0;
}

namespace
{
// Rows 1 to height - 1 of the legacy loops: the cell has its north and south neighbours.
bool hasnorthandsouth(const mesh::CellGrid& grid, CellId c)
{
    return grid.neighbour(c, 1) >= 0 && grid.neighbour(c, 5) >= 0;
}

// planet::surfacebelowsea at a neighbour: out of the grid reads as sea (map(x, y) is 0 there).
bool surfacebelowseaat(const planet& world, std::int32_t n)
{
    return n < 0 || world.surfacebelowsea(static_cast<CellId>(n));
}

// The sizes (cell counts) of the 4-connected components of sea (base at or below sea level); 0 on land.
vector<int> seaareas(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    vector<int> area(grid.cellcount(), 0);
    vector<char> checked(grid.cellcount(), 0);
    vector<CellId> order = columnmajorcells(grid);
    vector<CellId> component;

    for (auto it = order.rbegin(); it != order.rend(); ++it)
    {
        const CellId start = *it;
        if (!world.basebelowsea(start) || checked[start] || area[start] != 0)
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
                if (n < 0)
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

    const vector<int> area = seaareas(world);

    // Now turn all the points whose area is too small into land.

    for (const CellId c : columnmajorcells(grid))
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
vector<std::int32_t> straightrun(const mesh::CellGrid& grid, CellId c, int slot, Keep&& keep, Stop&& stopafter)
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
        at = step(grid, static_cast<CellId>(at), slot, turn);
        // A run that comes back round to its start (a whole parallel of coast) ends there.
        if (stopafter(at, turn) || at == static_cast<std::int32_t>(c) || cells.size() >= grid.cellcount())
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

    for (const CellId c : columnmajorcells(grid)) // First do vertical/horizontal lines.
    {
        if (!world.surfacebelowsea(c))
            continue;

        if (grid.neighbour(c, 1) >= 0 && !world.surfacebelowsea(static_cast<CellId>(grid.neighbour(c, 1)))) // Sea to the south, land to the north
        {
            const vector<std::int32_t> run = straightrun(grid, c, 3, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(1, turn)));
            }, never);
            disruptrun(run, c);
        }

        if (grid.neighbour(c, 5) >= 0 && !world.surfacebelowsea(static_cast<CellId>(grid.neighbour(c, 5)))) // Sea to the north, land to the south
        {
            const vector<std::int32_t> run = straightrun(grid, c, 3, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(5, turn)));
            }, never);
            disruptrun(run, c);
        }

        if (!surfacebelowseaat(world, grid.neighbour(c, 3))) // Sea to the west, land to the east
        {
            const vector<std::int32_t> run = straightrun(grid, c, 5, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(3, turn)));
            }, never);
            disruptrun(run, c);
        }

        if (!surfacebelowseaat(world, grid.neighbour(c, 7))) // Sea to the east, land to the west
        {
            const vector<std::int32_t> run = straightrun(grid, c, 5, [&](CellId at, int turn)
            {
                return world.surfacebelowsea(at) && !surfacebelowseaat(world, grid.neighbour(at, turned(7, turn)));
            }, never);
            disruptrun(run, c);
        }
    }

    auto nb = [&](CellId c, int slot) { return grid.neighbour(c, slot); };
    auto setfrom = [&](std::int32_t target, std::int32_t source)
    {
        if (target >= 0 && source >= 0)
            world.setnomfloat(static_cast<CellId>(target), world.nomfloat(static_cast<CellId>(source)));
    };

    for (const CellId c : columnmajorcells(grid)) // Now do diagonals.
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
    const vector<CellId> order = columnmajorcells(grid);

    auto land = [&](std::int32_t n) { return n >= 0 && !world.sea(static_cast<CellId>(n)); }; // a cell, and not sea

    for (const CellId c : order) // First remove any one-pixel channels.
    {
        if (!world.sea(c))
            continue;

        if (land(grid.neighbour(c, 1)) && land(grid.neighbour(c, 5)))
            lower(grid.neighbour(c, 1));

        if (land(grid.neighbour(c, 7)) && land(grid.neighbour(c, 3)))
            lower(grid.neighbour(c, 7));

        if (land(grid.neighbour(c, 8)) && land(grid.neighbour(c, 4)))
            lower(grid.neighbour(c, 8));

        if (land(grid.neighbour(c, 6)) && land(grid.neighbour(c, 2)))
            lower(grid.neighbour(c, 6));
    }

    for (const CellId c : order) // Now remove any two-pixel channels.
    {
        if (!world.sea(c))
            continue;

        if (land(grid.neighbour(c, 1)) && sea(grid.neighbour(c, 5)) == 1 && land(offset(grid, c, 0, 2)))
            lower(static_cast<std::int32_t>(c));

        if (land(grid.neighbour(c, 7)) && sea(grid.neighbour(c, 3)) == 1 && land(offset(grid, c, 2, 0)))
            lower(grid.neighbour(c, 3));
    }
}

void loweroceans(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    int sealevel = world.sealevel();

    for (const CellId c : columnmajorcells(grid))
    {
        if (world.sea(c) && world.nomfloat(c) > static_cast<float>(sealevel - 170))
        {
            bool nearland = 0;

            for (const CellId n : terrain::cellswithinkm(world, c, 2.5 * terrain::legacycellkm(world)))
            {
                if (world.sea(n) == 0)
                {
                    nearland = 1;
                    break;
                }
            }

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

    const vector<CellId> order = columnmajorcells(grid);

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

    // The sea cells within the reach (km, tuning::terrain::ocean_floor) of a coastal sea cell are
    // each averaged over the sea cells within the reach of themselves.
    const double reach_km = tuning::terrain::ocean_floor::coastSmoothingRadiusKm;

    vector<char> done(count, 0);

    auto smootharound = [&](CellId c)
    {
        for (const CellId kk : terrain::cellswithinkm(world, c, reach_km))
        {
            if (done[kk] == 0 && world.sea(kk) == 1 && coast(world, kk) == 0)
            {
                float total = 0.0f;
                int crount = 0;

                for (const CellId mm : terrain::cellswithinkm(world, kk, reach_km))
                {
                    if (world.sea(mm) == 1)
                    {
                        crount++;
                        total = total + world.nomfloat(mm);
                    }
                }

                if (total != 0)
                {
                    float newnom = total / static_cast<float>(crount);
                    world.setnomfloat(kk, newnom);
                    done[kk] = 1;

                    world.setnoshade(kk, 1);
                }
            }
        }
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

// This notes down all the very small islands on the map.

void checkislands(planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    int sealevel = world.sealevel();

    for (const CellId c : columnmajorcells(grid))
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

                    terrain::makemountainisland(world, c, volcano);

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

    for (const CellId c : columnmajorcells(grid))
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
