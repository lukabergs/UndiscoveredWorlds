#pragma once

// Recipes order processes; the scheduler checks them before running and can
// audit what each step actually writes (docs/architecture/world-model.md).

#include <string>
#include <string_view>
#include <vector>

#include "planet_bundle.hpp"
#include "world_processes.hpp"
#include "world_stages.hpp"

struct WorldRecipe
{
    std::string name;
    std::string description;
    std::vector<GenerationStageId> steps;
    std::vector<worldmodel::Component> frozen;      // components no step may write
    std::vector<worldmodel::FieldId> inputs;        // fields present before the first step
    std::vector<worldmodel::TableId> inputtables;
};

// procedural-world: every stage from plate tectonics.
WorldRecipe proceduralworldrecipe();
// imported-world: imported terrain, then the remaining terrain passes, climate and hydrology.
WorldRecipe importedworldrecipe();
// earth-climate: climate processes only on imported terrain; Geosphere is frozen.
WorldRecipe earthclimaterecipe();
bool findworldrecipe(std::string_view name, WorldRecipe& recipe);
std::vector<std::string> worldrecipenames();

struct RecipeIssue
{
    int step = -1;           // -1 for the recipe as a whole
    bool error = true;       // false: a warning that does not stop the run
    std::string message;
};

std::vector<RecipeIssue> validaterecipe(const WorldRecipe& recipe);

struct RecipeRunOptions
{
    bool auditAccess = false;    // hash every field around each step and check declared writes
};

struct RecipeStepRecord
{
    int index = 0;
    GenerationStageId stage{};
    std::string process;
    std::string label;
    double seconds = 0.0;
    bool measured = false;                          // writes come from the audit, not the declaration
    std::vector<worldmodel::FieldId> writes;
    std::vector<worldmodel::TableId> tablewrites;
    std::vector<std::uint64_t> writehashes;         // audit: hash of each written field after the step
    std::vector<worldmodel::FieldId> reads;         // analysis builds only
    std::vector<worldmodel::TableId> tablereads;
    bool drewrandom = false;                        // analysis builds only
    std::vector<std::string> violations;
};

struct RecipeRunResult
{
    bool completed = false;
    std::string error;
    std::vector<RecipeStepRecord> steps;
    BundleProvenance provenance;
    int violations = 0;
};

RecipeRunResult runrecipe(const WorldRecipe& recipe, planet& world, GenerationScratch& scratch,
    WorldStageParameters& parameters, const GenerationExecutionContext& context, const RecipeRunOptions& options);

// Bundle step records for a finished run.
std::vector<runbundle::StepRecord> bundlesteps(const RecipeRunResult& result);
