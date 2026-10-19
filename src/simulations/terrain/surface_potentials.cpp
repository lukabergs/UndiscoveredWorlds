#include "simulations/terrain/surface_potentials.hpp"
#include "core/grid_shapes.hpp"
#include "core/mesh/cell_grid.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/resources/physical_layers_internal.hpp"
#include "world/cell_queries.hpp"

namespace physical_layers
{
using namespace detail;
// The shared per-cell queries (world/cell_queries.hpp); clampunit and clampscore are detail's.
using worldcells::annualrainfall;
using worldcells::averageriverflow;
using worldcells::isdeltacell;
using worldcells::islakecell;
using worldcells::issaltfeature;
using worldcells::iswetlandcell;

namespace
{
float averagetemperature(const planet& world, mesh::CellId c)
{
    return (world.mintempfloat(c) + world.maxtempfloat(c)) / 2.0f;
}

float localrelief(const planet& world, mesh::CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    const float centre = world.mapfloat(c);
    float maxdifference = 0.0f;

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour < 0)
            continue;

        maxdifference = std::max(maxdifference, std::abs(world.mapfloat(static_cast<mesh::CellId>(neighbour)) - centre));
    }

    return maxdifference;
}

// The neighbours whose river or distributary direction leads into the cell.
int inflowcount(const planet& world, mesh::CellId c)
{
    const mesh::CellGrid& grid = world.grid();
    int inflow = 0;

    for (int slot = 1; slot <= mesh::slotcount; slot++)
    {
        const std::int32_t neighbour = grid.neighbour(c, slot);

        if (neighbour < 0)
            continue;

        const mesh::CellId n = static_cast<mesh::CellId>(neighbour);
        const int riverdir = world.riverdir(n);

        if (riverdir != 0 && grid.neighbour(n, riverdir) == static_cast<std::int32_t>(c))
            inflow++;

        const int deltadir = world.deltadir(n);

        if (deltadir > 0 && grid.neighbour(n, deltadir) == static_cast<std::int32_t>(c))
            inflow++;
    }

    return inflow;
}

float temperaturerangefactor(float temperature)
{
    if (temperature <= -10.0f || temperature >= 40.0f)
        return 0.0f;

    if (temperature <= 18.0f)
        return clampunit((temperature + 10.0f) / 28.0f);

    return clampunit(1.0f - (temperature - 18.0f) / 22.0f);
}

float rainfertilityfactor(float rainfall)
{
    if (rainfall <= 75)
        return 0.0f;

    if (rainfall <= 700)
        return clampunit(static_cast<float>(rainfall - 75) / 625.0f);

    if (rainfall <= 1800)
        return 1.0f;

    return clampunit(1.0f - static_cast<float>(rainfall - 1800) / 1200.0f);
}

// 100 x log10(1 + rate) / log10(1 + reference), rates in m/Myr; the caller clamps to 0..100.
float ratepotential(float rate, float reference)
{
    const float denominator = std::log10(1.0f + std::max(0.0f, reference));
    if (denominator <= 0.0f)
        return 0.0f;
    return 100.0f * std::log10(1.0f + std::max(0.0f, rate)) / denominator;
}
}

void generate_surface_potentials(planet& world)
{
    const int maximumflow = std::max(1, world.maxriverflow());
    const int majorriverthreshold = std::max(250, maximumflow / 20);
    // An eroded world (erosion_duration > 0) takes the erosion and deposition
    // potentials from the erosion run's rates; a procedural world keeps the proxies.
    const planet& reader = world;
    const bool eroded = reader.erosionduration() > 0.0f;
    // Net rates over the run (B4 review): the gross erosion and deposition rates of a
    // cell are both large where sediment passes through, so the scores read the net
    // rock removal and the net accumulation instead (denudation minus aggradation and
    // the reverse: rock removed and sediment laid, not the surface's movement, which
    // the isostatic response of landscape.isostasy.* makes smaller).
    const worldgrid<float>* denudation = eroded ? &reader.denudation() : nullptr;
    const worldgrid<float>* aggradation = eroded ? &reader.aggradation() : nullptr;
    const float duration = eroded ? reader.erosionduration() : 1.0f;
    const float erosionreference = tuning::landscape::potentials::erosionReferenceRateMMyr;
    const float depositionreference = tuning::landscape::potentials::depositionReferenceRateMMyr;

    parallelforrows(0, static_cast<int>(world.cellcount()) - 1, [&](int first, int last)
    {
        for (mesh::CellId c = static_cast<mesh::CellId>(first); c <= static_cast<mesh::CellId>(last); c++)
        {
            if (world.sea(c) != 0)
            {
                world.seterosionpotential(c, 0);
                world.setdepositionpotential(c, 0);
                world.setfloodplainfertility(c, 0);
                continue;
            }

            const float relief = localrelief(world, c);
            const float reliefnorm = clampunit(relief / 2200.0f);
            const float lownorm = 1.0f - clampunit(relief / 1400.0f);
            const float rainfallnorm = clampunit(static_cast<float>(annualrainfall(world, c)) / 1000.0f);
            const float flownorm = std::log1p(static_cast<float>(averageriverflow(world, c))) / std::log1p(static_cast<float>(maximumflow));
            const float inflownorm = clampunit(static_cast<float>(inflowcount(world, c)) / 4.0f);
            // Discs of one and two legacy pixels (cellqueries::anywithinpixels), as the resources read.
            const bool majorrivernear = cellqueries::anywithinpixels(world, c, 2, [&](mesh::CellId n)
            {
                return averageriverflow(world, n) >= majorriverthreshold || isdeltacell(world, n);
            });
            const bool wetlandnear = cellqueries::anywithinpixels(world, c, 1, [&](mesh::CellId n)
            {
                return iswetlandcell(world, n);
            });
            const bool lakemargin = islakecell(world, c) || cellqueries::anywithinpixels(world, c, 1, [&](mesh::CellId n)
            {
                return islakecell(world, n) || issaltfeature(world, n);
            });
            const bool terminalbasin = world.basinclass(c) == BasinClass::endorheic && (issaltfeature(world, c) || lakemargin);
            const float rivercorridor = std::max(flownorm, majorrivernear ? 0.65f : 0.0f);

            float erosionscore = 100.0f * (0.34f * reliefnorm + 0.30f * flownorm + 0.24f * rainfallnorm + 0.12f * inflownorm);

            if (world.mountainheightfloat(c) > 0)
                erosionscore += 8.0f;

            float depositionscore = 100.0f * (0.32f * lownorm + 0.27f * rivercorridor + 0.15f * rainfallnorm
                + 0.12f * (wetlandnear ? 1.0f : 0.0f)
                + 0.09f * (lakemargin ? 1.0f : 0.0f)
                + 0.05f * (terminalbasin ? 1.0f : 0.0f));

            if (isdeltacell(world, c))
                depositionscore += 20.0f;

            if (reliefnorm > 0.7f)
                depositionscore *= 0.6f;

            if (eroded)
            {
                const float net = ((*denudation)[c] - (*aggradation)[c]) / duration;
                erosionscore = ratepotential(std::max(0.0f, net), erosionreference);
                depositionscore = ratepotential(std::max(0.0f, -net), depositionreference);
            }

            const float fertilitytemperature = temperaturerangefactor(averagetemperature(world, c));
            const float fertilityrain = rainfertilityfactor(annualrainfall(world, c));
            const float riverproximity = std::max(rivercorridor, majorrivernear ? 1.0f : 0.0f);
            const float wetlandfactor = wetlandnear ? 1.0f : (lakemargin ? 0.5f : 0.0f);
            const float gentleterrain = 1.0f - clampunit(relief / 1800.0f);
            const float coldpenalty = clampunit(static_cast<float>(world.glacialtemp() - world.maxtempfloat(c) + 2) / 14.0f);
            const float aridpenalty = clampunit(static_cast<float>(200 - annualrainfall(world, c)) / 200.0f);
            const float glacialpenalty = (world.biome(c) == biomeice || world.maxtempfloat(c) <= world.glacialtemp()) ? 1.0f : 0.0f;

            float fertilityscore = 0.42f * static_cast<float>(clampscore(depositionscore)) + 22.0f * riverproximity + 14.0f * fertilityrain
                + 12.0f * fertilitytemperature + 10.0f * wetlandfactor;
            fertilityscore *= 0.55f + 0.45f * gentleterrain;
            fertilityscore -= 24.0f * coldpenalty + 20.0f * aridpenalty + 26.0f * glacialpenalty;

            world.seterosionpotential(c, clampscore(erosionscore));
            world.setdepositionpotential(c, clampscore(depositionscore));
            world.setfloodplainfertility(c, clampscore(fertilityscore));
        }
    }, 1024);
}
}
