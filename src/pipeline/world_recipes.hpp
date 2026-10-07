#pragma once

// Recipes order processes; the scheduler checks them before running and can
// audit what each step actually writes (docs/architecture/world-model.md).
// A step is a process invocation with a label and its own settings; a process
// may appear more than once (procedural-world runs the climate twice).

#include <functional>
#include <string>
#include <string_view>
#include <vector>

#include "world/planet_bundle.hpp"
#include "pipeline/world_processes.hpp"
#include "pipeline/world_stages.hpp"

// A setting one step changes while it runs, restored after it: a module parameter
// named as --set takes it (<module>.<group>.<key>), or a climate-context setting,
// climate.resolution (as --climate-resolution) or climate.coupling_iterations (as
// --climate-coupling-iterations). It wins over the run's value for that step only.
struct RecipeOverride
{
    std::string name;
    std::string value;
};

struct RecipeStep
{
    GenerationStageId stage{};
    // Unique within the recipe: the process id, or the id with a suffix
    // ("climate.rainfall_and_coupling@coarse") on every occurrence of a repeated process.
    std::string label;
    std::vector<RecipeOverride> overrides;
};

// A step labelled with its process id, or with "<process id>@<suffix>" when a suffix is given.
RecipeStep recipestep(GenerationStageId stage, std::string_view suffix = {}, std::vector<RecipeOverride> overrides = {});

struct WorldRecipe
{
    std::string name;
    std::string description;
    std::vector<RecipeStep> steps;
    std::vector<worldmodel::Component> frozen;      // components no step may write
    std::vector<worldmodel::FieldId> inputs;        // fields present before the first step
    std::vector<worldmodel::TableId> inputtables;
};

// True when a step of the recipe runs this process.
bool recipehasstage(const WorldRecipe& recipe, GenerationStageId stage);

// procedural-world: the terrain stages, a coarse climate (@coarse: 32 columns, 2 coupling
// iterations), the fluvial erosion, the flow routing and lake levels
// (@before-climate), the full climate at the run's settings (@full; its pressure and
// winds start from the coarse climate's sea-surface temperatures), then the flow
// routing, the rivers, the lake levels (@after-climate) and the remaining stages.
WorldRecipe proceduralworldrecipe();
// imported-world: imported terrain, then the remaining terrain passes, climate,
// the flow routing, rivers and lake levels (no erosion), and the remaining stages.
WorldRecipe importedworldrecipe();
// earth-climate: climate processes only on imported terrain; Geosphere is frozen.
WorldRecipe earthclimaterecipe();
// earth-hydrology: earth-climate followed by the flow routing and the lake levels; Geosphere is frozen.
WorldRecipe earthhydrologyrecipe();
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

// One step as an index, a label ("climate.rainfall_and_coupling@coarse") or the id of
// a process that occurs once ("climate.winds"); an id that occurs more than once is
// refused, and `reason` names its labels.
bool findrecipestep(const WorldRecipe& recipe, std::string_view text, std::size_t& index, std::string* reason = nullptr);
// The evolver of the time model named by its process id or by the label of one of
// its steps in the recipe; null, with `reason`, otherwise. A duration or checkpoint
// given for it applies to the process, so to every step of it in the recipe.
const WorldProcessDefinition* findtimedevolver(const WorldRecipe& recipe, std::string_view label,
    std::string* reason = nullptr);
// A label as a file or directory name: characters other than letters, digits and ._@- become '_'.
std::string recipesteppathname(std::string_view label);
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
    std::string label;                              // the recipe step's label
    std::vector<RecipeOverride> overrides;          // as applied: values in their canonical text
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

// Runs the recipe's steps in order; a step's overrides hold while it runs, and the
// values they replaced are back before the next step and before `afterstep`.
RecipeRunResult runrecipe(const WorldRecipe& recipe, planet& world, GenerationScratch& scratch,
    WorldStageParameters& parameters, const GenerationExecutionContext& context, const RecipeRunOptions& options);

// Bundle step records for a finished run.
std::vector<runbundle::StepRecord> bundlesteps(const RecipeRunResult& result);
