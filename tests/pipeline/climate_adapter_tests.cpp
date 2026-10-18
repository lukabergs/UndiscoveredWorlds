// The climate adapter (docs/architecture/sphere-mesh.md section 10, T8 W2-C): on a
// cube world the climate processes run on a lat-lon shadow of the export width;
// only changed inputs are re-exported and each process's writes are imported once.
// On a lat-lon world the adapter is the identity (no shadow). Seed 1.

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "core/mesh/cell_grid.hpp"
#include "core/world/world_fields.hpp"
#include "pipeline/climate_adapter.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"
#include "world/planet.hpp"
#include "world/planet_fields.hpp"
#include "world/world_setup.hpp"

namespace
{
using S = GenerationStageId;
using worldmodel::FieldId;

// The coarse climate of the procedural recipe, then the two climate diagnosers.
const std::vector<GenerationStageId> kClimateProcesses = { S::global_temperature, S::pressure, S::winds, S::ocean_currents,
    S::sea_surface_temperatures, S::sea_ice_and_tides, S::rainfall, S::mountain_temperature_lapse, S::climates_and_biomes };

// Two continents and an ocean on a world of the given grid, seed 1, at the coarse
// climate's settings (32 columns, 2 coupling iterations).
std::unique_ptr<planet> syntheticworld(mesh::GridKind kind, int columns)
{
    auto world = std::make_unique<planet>();
    initialiseworld(*world);
    world->setseed(1);
    world->setgrid(kind, columns);
    const mesh::CellGrid& grid = world->grid();
    const mesh::Vec3 first = mesh::unitvector(0.3, 0.5);
    const mesh::Vec3 second = mesh::unitvector(-0.4, 2.8);
    const float sealevel = static_cast<float>(world->sealevel());
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        const double a = p.x * first.x + p.y * first.y + p.z * first.z;
        const double b = p.x * second.x + p.y * second.y + p.z * second.z;
        const double land = std::max(a - 0.75, b - 0.8);
        world->setnomfloat(c, land > 0.0 ? sealevel + 200.0f + static_cast<float>(land) * 8000.0f : sealevel - 3000.0f);
    }
    getlandandseatotals(*world);
    climatecontext::Context& climate = world->climatesimulation();
    climate.resolution = { 32, 32, 64, 16 };
    climate.couplingIterationLimit = 2;
    return world;
}

bool runstage(GenerationStageId stage, planet& world, GenerationScratch& scratch)
{
    WorldStageParameters parameters;
    parameters.seaLevel = world.sealevel();
    const GenerationExecutionContext context;
    std::string error;
    const bool ran = executegenerationstage(stage, world, scratch, parameters, context, &error);
    EXPECT_TRUE(ran) << getworldprocess(stage).id << ": " << error;
    return ran;
}
}

TEST(ClimateAdapter, CoarseClimateOnCubeWorld)
{
    std::unique_ptr<planet> world = syntheticworld(mesh::GridKind::cube, 64);   // N 16
    const mesh::CellGrid& grid = world->grid();
    ASSERT_EQ(grid.kind(), mesh::GridKind::cube);
    GenerationScratch scratch;
    initializegenerationscratch(scratch);

    std::vector<FieldId> imported;
    for (std::size_t index = 0; index < kClimateProcesses.size(); index++)
    {
        const GenerationStageId stage = kClimateProcesses[index];
        ASSERT_TRUE(runstage(stage, *world, scratch));
        const ClimateShadowSync& sync = scratch.climateShadow.lastsync();
        const WorldProcessDefinition& process = getworldprocess(stage);
        std::cout << "ClimateAdapter " << std::left << std::setw(34) << process.id << std::fixed << std::setprecision(4)
                  << " export " << sync.exportseconds << " s (" << sync.exported << (sync.created ? ", created" : "")
                  << ")  import " << sync.importseconds << " s (" << sync.imported << " fields, " << sync.keptcells
                  << " cells kept)\n";
        // The shadow is built once, by the first process; afterwards nothing changed on the
        // cube but climate's own imports, so nothing is sent back.
        EXPECT_EQ(sync.created, index == 0) << process.id;
        EXPECT_EQ(sync.exported, 0u) << process.id;
        // At most one import per process, of the declared writes it changed; every
        // process here but the lapse (no mountains) changes some.
        EXPECT_EQ(sync.importcalls, sync.imported > 0 ? 1u : 0u) << process.id;
        if (stage != S::mountain_temperature_lapse)
            EXPECT_GT(sync.imported, 0u) << process.id;
        EXPECT_LE(sync.imported, process.writes.size()) << process.id;
        imported.insert(imported.end(), process.writes.begin(), process.writes.end());
    }
    const planet* shadow = scratch.climateShadow.world();
    ASSERT_NE(shadow, nullptr);
    EXPECT_EQ(shadow->grid().kind(), mesh::GridKind::latlon);
    EXPECT_EQ(shadow->grid().columns(), 64);

    // Every imported float is finite on every cube cell.
    for (const FieldId id : imported)
    {
        const worldmodel::FieldSpec& spec = worldmodel::fieldspec(id);
        const PlanetFieldData data = readplanetfield(*world, id);
        if (spec.element != worldmodel::ElementType::f32 || spec.layout == worldmodel::Layout::scalar ||
            spec.layout == worldmodel::Layout::zonalprofile || data.shape.empty() || data.shape.back() != grid.cellcount())
            continue;
        std::vector<float> values(data.bytes.size() / sizeof(float));
        std::memcpy(values.data(), data.bytes.data(), data.bytes.size());
        std::size_t bad = 0;
        for (const float value : values)
            bad += std::isfinite(value) ? 0 : 1;
        EXPECT_EQ(bad, 0u) << worldmodel::qualifiedname(id);
    }

    // The annual temperature falls from the tropics to the mid-latitudes to the poles.
    double sums[3] = {}, counts[3] = {};
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
    {
        const double latitude = std::abs(grid.latitude_rad(c)) * 180.0 / 3.14159265358979323846;
        const int band = latitude < 15.0 ? 0 : (latitude > 30.0 && latitude < 45.0) ? 1 : latitude > 60.0 ? 2 : -1;
        if (band < 0)
            continue;
        const float annual = 0.5f * (world->jantempfloat(c) + world->jultempfloat(c));
        ASSERT_TRUE(std::isfinite(annual));
        sums[band] += annual;
        counts[band] += 1.0;
    }
    ASSERT_GT(counts[0] * counts[1] * counts[2], 0.0);
    const double tropics = sums[0] / counts[0], middle = sums[1] / counts[1], polar = sums[2] / counts[2];
    std::cout << "ClimateAdapter annual temperature: tropics " << tropics << ", 30-45 " << middle << ", polar " << polar << '\n';
    EXPECT_GT(tropics, middle);
    EXPECT_GT(middle, polar);

    // A changed input is re-exported, and only it: raise the land, run the temperature again.
    for (mesh::CellId c = 0; c < grid.cellcount(); c++)
        if (!world->basebelowsea(c))
            world->setnomfloat(c, world->nomfloat(c) + 300.0f);
    ASSERT_TRUE(runstage(S::global_temperature, *world, scratch));
    EXPECT_FALSE(scratch.climateShadow.lastsync().created);
    EXPECT_EQ(scratch.climateShadow.lastsync().exported, 1u);   // geosphere.base_elevation
    // Unchanged again: nothing to send.
    ASSERT_TRUE(runstage(S::pressure, *world, scratch));
    EXPECT_EQ(scratch.climateShadow.lastsync().exported, 0u);
}

TEST(ClimateAdapter, IdentityOnLatlonWorld)
{
    std::unique_ptr<planet> world = syntheticworld(mesh::GridKind::latlon, 64);
    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    ASSERT_TRUE(runstage(S::global_temperature, *world, scratch));
    EXPECT_EQ(scratch.climateShadow.world(), nullptr);
    EXPECT_EQ(scratch.climateShadow.lastsync().importcalls, 0u);
}
