#pragma once

// Recipes order processes; the scheduler checks them before running and can
// audit what each step actually writes (docs/architecture/world-model.md).

#include <functional>
#include <string>
#include <string_view>
#include <vector>

#include "world/planet_bundle.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"

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
// eroded-world: procedural-world with flow routing and fluvial erosion after the climate coupling.
WorldRecipe erodedworldrecipe();
// earth-climate: climate processes only on imported terrain; Geosphere is frozen.
WorldRecipe earthclimaterecipe();
// tectonics: the Geosphere initializer and the tectonics evolver on the engine's grid.
WorldRecipe tectonicsrecipe();
// True for the recipes that start from imported terrain.
bool recipeimportsterrain(const WorldRecipe& recipe);
bool findworldrecipe(std::string_view name, WorldRecipe& recipe);
std::vector<std::string> worldrecipenames();

struct RecipeIssue
{
    int step = -1;           // -1 for the recipe as a whole
    bool error = true;       // false: a warning that does not stop the run
    std::string message;
};

std::vector<RecipeIssue> validaterecipe(const WorldRecipe& recipe);

// A contiguous part of a recipe, by step index (inclusive).
struct RecipeSteps
{
    std::size_t first = 0;
    std::size_t last = 0;
};

// One step as an index or a process id ("climate.winds").
bool findrecipestep(const WorldRecipe& recipe, std::string_view text, std::size_t& index);
// "a..b", "a..", "..b" or one step.
bool parserecipesteps(const WorldRecipe& recipe, std::string_view text, RecipeSteps& steps, std::string* reason = nullptr);

struct RecipeStepRecord;
struct RecipeRunResult;

struct RecipeRunOptions
{
    bool auditAccess = false;    // hash every field around each step and check declared writes
    bool selectSteps = false;    // run only `steps`; otherwise the whole recipe
    RecipeSteps steps;
    // World state loaded from a bundle: the fields and tables it supplied stand
    // in for the recipe's inputs and the steps before the selection.
    bool fromBundle = false;
    std::vector<worldmodel::FieldId> loadedfields;
    std::vector<worldmodel::TableId> loadedtables;
    BundleProvenance sourceprovenance;              // last writers recorded by that bundle
    // Called after each completed step (checkpoints); false stops the run.
    std::function<bool(const RecipeRunResult& sofar, const RecipeStepRecord& step, std::string& error)> afterstep;
};

// Checks that a run of part of a recipe has what its steps read: a field or
// table the full recipe would have provided must come from the bundle or from a
// step of this run, and no step may depend on state that an earlier, skipped
// step leaves outside the world fields.
std::vector<RecipeIssue> validaterecipesteps(const WorldRecipe& recipe, const RecipeRunOptions& options);

struct RecipeStepRecord
{
    int index = 0;                                  // position in the recipe
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
