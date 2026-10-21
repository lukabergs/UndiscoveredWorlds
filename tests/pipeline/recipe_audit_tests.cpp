// The access audit of runrecipe (options.auditAccess, master plan 2.6): a write a step does not
// declare is counted as a violation. The undeclared write is injected through the public
// afterstep hook, which changes hydrosphere.runoff between two runs of terrain.subgrid_statistics
// (declared writes: relief_amplitude and valley_spacing); the second step's audit then sees the
// change as its own. Cube of N 16, seed 1.

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>

#include <gtest/gtest.h>

#include "core/mesh/cell_grid.hpp"
#include "core/world/world_fields.hpp"
#include "pipeline/world_recipes.hpp"
#include "pipeline/world_stages.hpp"
#include "world/planet.hpp"
#include "world/planet_fields.hpp"
#include "world/world_parameters.hpp"

namespace
{
using S = GenerationStageId;
using worldmodel::FieldId;

constexpr int facecells = 16;
constexpr float sealevel = 31000.0f;

// Land where x > 0.1 (rising inland), sea elsewhere; uniform erodibility, runoff and erosion rate.
std::unique_ptr<planet> auditworld()
{
    auto world = std::make_unique<planet>();
    world->setgrid(mesh::GridKind::cube, 4 * facecells);
    world->setseed(1);
    world->setmaxelevation(48000);
    world->setsealevel(static_cast<int>(sealevel));
    world->setgravity(1.0f);
    const mesh::CellGrid& grid = world->grid();
    for (mesh::CellId c = 0; c < world->cellcount(); c++)
    {
        const double land = grid.centre(c).x - 0.1;
        world->setnomfloat(c, static_cast<float>(sealevel + (land > 0.0 ? 50.0 + 2500.0 * land : -3000.0)));
        world->erodibility()[c] = 1.0f;
        world->runoff()[c] = 500.0f;
        world->erosionrate()[c] = 100.0f;
    }
    return world;
}

WorldRecipe twostatistics()
{
    WorldRecipe recipe;
    recipe.name = "audit-test";
    recipe.steps = { recipestep(S::subgrid_statistics, "first"), recipestep(S::subgrid_statistics, "second") };
    return recipe;
}

// Runs the recipe audited; with inject, the afterstep hook changes the runoff after the first step.
RecipeRunResult auditedrun(bool inject)
{
    const std::unique_ptr<planet> world = auditworld();
    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    WorldStageParameters parameters;
    parameters.seaLevel = world->sealevel();
    GenerationExecutionContext context;
    RecipeRunOptions options;
    options.auditAccess = true;
    planet& target = *world;
    options.afterstep = [&target, inject](const RecipeRunResult&, const RecipeStepRecord& step, std::string&)
    {
        if (inject && step.index == 0)
            target.runoff()[mesh::CellId{ 0 }] += 1.0f;
        return true;
    };
    return runrecipe(twostatistics(), target, scratch, parameters, context, options);
}

bool hasviolation(const RecipeStepRecord& step, const std::string& text)
{
    return std::find(step.violations.begin(), step.violations.end(), text) != step.violations.end();
}
}

TEST(RecipeAudit, UndeclaredWriteIsAViolation)
{
    const RecipeRunResult result = auditedrun(true);
    ASSERT_TRUE(result.completed) << result.error;
    ASSERT_EQ(result.steps.size(), 2u);
    EXPECT_TRUE(result.steps[1].measured);
    EXPECT_TRUE(hasviolation(result.steps[1], "undeclared write " + worldmodel::qualifiedname(FieldId::runoff)));
    EXPECT_GT(result.violations, 0);
}

TEST(RecipeAudit, DeclaredWritesAloneAreClean)
{
    // The same run without the injected write: the second step's audit finds nothing undeclared.
    const RecipeRunResult result = auditedrun(false);
    ASSERT_TRUE(result.completed) << result.error;
    ASSERT_EQ(result.steps.size(), 2u);
    EXPECT_TRUE(result.steps[1].measured);
    EXPECT_FALSE(hasviolation(result.steps[1], "undeclared write " + worldmodel::qualifiedname(FieldId::runoff)));
    for (const std::string& violation : result.steps[1].violations)
        EXPECT_EQ(violation.rfind("note: ", 0), 0u) << violation;
}
