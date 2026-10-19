#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "core/classes.hpp"
#include "core/deterministic_random.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"
#include "simulations/hydrology/flow_queries.hpp"
#include "simulations/hydrology/hydrology_cells.hpp"
#include "simulations/hydrology/wetlands.hpp"
#include "simulations/terrain/terrain_elevation.hpp"

using namespace std;

// On the world's cell grid: neighbourhoods from the neighbour table, the patches placed in km
// through the cell's tangent frame, the drainage noise on the sphere.

namespace
{
// Salt of the wetlands' drainage noise (terrain::spherefractal).
constexpr std::uint64_t drainagesalt = 0x7303;

// The cells within two steps of the neighbour table (the 5 x 5 box of the lat-lon grid).
vector<mesh::CellId> twosteps(const mesh::CellGrid& grid, mesh::CellId c)
{
    vector<mesh::CellId> cells = hydrologycells::square(grid, c);
    const size_t first = cells.size();
    for (size_t index = 0; index < first; index++)
        for (const mesh::CellId cell : hydrologycells::square(grid, cells[index]))
            cells.push_back(cell);
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());
    return cells;
}

// This actually pastes a patch of wetlands onto the map: the template's pixels are legacy
// pixels (terrain::legacycellkm, a cell of the 2048-column world) laid out east and north of
// the centre cell in its tangent frame, top row north.

void pastewetland(planet& world, mesh::CellId centre, int shapenumber, boolshapetemplate smalllake[])
{
    const mesh::CellGrid& grid = world.grid();
    const double pixel_km = terrain::legacycellkm(world);

    int imheight = smalllake[shapenumber].ysize() - 1;
    int imwidth = smalllake[shapenumber].xsize() - 1;

    // The legacy left-right and top-bottom reversals: drawn as before, they changed only the
    // order in which the pixels were pasted, not the patch.
    random(0, 1);
    random(0, 1);

    const mesh::TangentFrame frame = mesh::tangentframe(grid.centre(centre), grid.radius_m());

    for (int i = 0; i <= imwidth; i++)
    {
        for (int j = 0; j <= imheight; j++)
        {
            if (smalllake[shapenumber].point(i, j) == 1)
            {
                const double east_km = static_cast<double>(i - imwidth / 2) * pixel_km;
                const double north_km = -static_cast<double>(j - imheight / 2) * pixel_km;
                const mesh::CellId c = grid.locate(mesh::fromtangent(frame, east_km, north_km));

                if (c == mesh::nocell)
                    continue;

                if (world.sea(c) == 0 && world.mountainheightfloat(c) == 0) // Don't try to put wetlands on top of sea or mountains...
                    world.setspecial(c, 130);
            }
        }
    }
}
}

// This lays down wetlands.

void createwetlands(planet& world, boolshapetemplate smalllake[])
{
    const mesh::CellGrid& grid = world.grid();
    const vector<mesh::CellId> order = hydrologycells::visitorder(grid);
    int maxelev = world.maxelevation();

    int minflow = 100; // Average flow must be at least this to have wetlands.
    int maxflatness = 5; // Must be no less flat than this.
    float factor = (float)maxelev / 5.0f;

    int shapenumber;

    // First, the drainage map: the legacy grain-8 fractal (valuemod 0.2 to 3-6) on the sphere,
    // as terrain::terrainfractal stands for it.

    const vector<float> drainage = terrain::spherefractal(world, deterministiccontextseed(world.seed(), drainagesalt),
                                                          terrain::terrainfractal(), static_cast<float>(maxelev));

    // Now go through the map and place wetlands.

    for (const mesh::CellId c : order)
    {
        if (world.sea(c) == 0 && world.special(c) == 0 && world.lakesurfacefloat(c) == 0 && world.climate(c) != 31 && world.climate(c) != 5)
        {
            if (world.deltadir(c) != 0)
            {
                shapenumber = random(0, 5);
                pastewetland(world, c, shapenumber, smalllake);
            }
            else
            {
                if (world.riveraveflow(c) >= minflow)
                {
                    int flatness = getflatness(world, c);

                    if (flatness <= maxflatness)
                    {
                        float drain = (float)(int)drainage[c];

                        drain = drain / factor;

                        int d = (int)drain; // This should be from 1 to 5.
                        d = d * 2;

                        int prob = flatness + d; // The flatter the land, the more likely wetlands are.

                        int flow = world.riveraveflow(c) / 1000;

                        prob = prob - flow; // The bigger the flow here, the more likely wetlands are.

                        for (const mesh::CellId cell : hydrologycells::square(grid, c)) // The more lakes/sea there are around here, the more likely wetlands are.
                        {
                            if (world.lakesurfacefloat(cell) != 0)
                                prob = prob - 3;

                            if (world.sea(cell) == 1)
                                prob = prob - 4;
                        }

                        prob = prob + static_cast<int>(
                            40.0f - world.averainfloat(c) / 10.0f);

                        if (world.climate(c) == 30 && world.maxtempfloat(c) >= 5) // Much likelier in tundra with warmish summers
                            prob = prob / 4;

                        if (prob < 1)
                            prob = 1;

                        if (random(1, prob) == 1) // Put some wetlands here.
                        {
                            shapenumber = random(0, 5);
                            pastewetland(world, c, shapenumber, smalllake);
                        }
                    }
                }
            }
        }
    }

    int deltaswampchance = 6; // The BIGGER this is, the more probable wetlands will be around delta branches.

    for (const mesh::CellId c : order) // Add wetlands to the actual delta branches
    {
        if (world.deltadir(c) > 0 && world.sea(c) == 0 && random(1, deltaswampchance) != 1)
            world.setspecial(c, 130);
    }

    // Now fill in extra wetlands if necessary to reach the sea or other wetlands: a cell between
    // a wetland and a wetland or the sea across it (the opposite slots of its own frame).

    const auto specialat = [&](std::int32_t cell) { return cell < 0 ? 0 : world.special(static_cast<mesh::CellId>(cell)); };
    const auto seaat = [&](std::int32_t cell) { return cell >= 0 && world.sea(static_cast<mesh::CellId>(cell)) == 1; };

    for (const mesh::CellId c : order)
    {
        if (world.sea(c) == 0 && world.special(c) == 0 && world.climate(c) != 31)
        {
            // The legacy order: north-south, south-north, west-east, east-west, then the diagonals.
            for (const int slot : { 1, 5, 7, 3, 8, 2, 6, 4 })
            {
                const std::int32_t here = grid.neighbour(c, slot);
                const std::int32_t across = grid.neighbour(c, (slot + 3) % 8 + 1);

                if (specialat(here) == 130)
                {
                    if (specialat(across) == 130)
                        world.setspecial(c, 130);

                    if (seaat(across))
                        world.setspecial(c, 130);
                }
            }
        }
    }

    for (const mesh::CellId c : order) // For some reason it keeps putting wetlands all over the sea, so we remove them.
    {
        if (world.sea(c) == 1)
            world.setspecial(c, 0);
    }

    // Now make the wetlands salty/brackish if necessary.

    for (const mesh::CellId c : order)
    {
        if (world.special(c) == 130)
        {
            int salt = 0;

            for (const mesh::CellId cell : hydrologycells::square(grid, c))
            {
                if (world.sea(cell) == 1)
                    salt = 1;

                if (world.special(cell) == 100 || world.special(cell) == 110)
                    salt = 1;
            }

            if (salt == 1)
                world.setspecial(c, 132);
            else
            {
                for (const mesh::CellId cell : twosteps(grid, c))
                {
                    if (world.sea(cell) == 1)
                        salt = 1;

                    if (world.special(cell) == 100 || world.special(cell) == 110)
                        salt = 1;
                }

                if (salt == 1)
                    world.setspecial(c, 131);
            }
        }
    }
}

// This removes any wetlands tiles that border large lakes.

void removeexcesswetlands(planet& world)
{
    const mesh::CellGrid& grid = world.grid();

    for (const mesh::CellId c : hydrologycells::visitorder(grid))
    {
        if (world.special(c) == 130)
        {
            bool nearlake = 0;

            for (const mesh::CellId cell : hydrologycells::square(grid, c))
            {
                if (world.truelake(cell) == 1)
                {
                    nearlake = 1;
                    break;
                }
            }

            if (nearlake == 1)
            {
                world.setlakesurface(c, 0);
                world.setspecial(c, 0);
            }
        }
    }
}
