#pragma once

// Headless world generation entry points shared
// by world_pipeline and the desktop application.

#include <filesystem>
#include <string>
#include <vector>

#include "cli/climate_options.hpp"
#include "cli/tectonic_options.hpp"
#include "world/planet.hpp"
#include "simulations/society/social_generation.hpp"

struct ImportedClimateMaps;

constexpr int kEarthBenchmarkWidth = 512;

// A time given for a named evolver: the label (the process id, or the label of one
// of its steps in the recipe) is resolved once the recipe is known.
struct EvolverTimeRequest
{
    std::string label;
    std::string text;     // the time as given
    double timeMyr = 0;
};

struct CommandLineGenerationOptions
{
    bool run = false;
    bool showHelp = false;
    bool deltas = true;
    bool socialEnabled = false;
    SocialGenerationOptions::Mode socialMode = SocialGenerationOptions::Mode::static_ex_nihilo;
    bool usePrehistory = true;
    int historyYears = 1200;
    bool climateOnly = false;
    bool hasSeed = false;
    bool hasSavePath = false;
    bool hasImportLandPath = false;
    bool hasImportSeaPath = false;
    bool hasWorldWidth = false;
    bool hasWorldHeight = false;
    bool hasResolution = false;
    long seed = 0;
    double plateTectonicsDurationMyr = 0;
    int plateTectonicsPlateCount = 0;
    int worldWidth = 0;
    int worldHeight = 0;
    int resolution = 0;
    // --grid latlon|cube: the world's cell grid (docs/architecture/sphere-mesh.md);
    // the cube has N = W / 4 cells per face edge, so W must be a multiple of 4.
    mesh::GridKind grid = mesh::GridKind::latlon;
    ClimateRunOptions climate;
    std::string savePath;
    std::string bundlePath;      // Run bundle directory (--bundle); empty writes none.
    std::string recipeName;      // --recipe; empty selects the recipe from the command and imports.
    bool auditAccess = false;    // --audit-access: hash fields around each step and check declared writes.
    bool runCommand = false;     // run: the recipe is named, nothing is registered.
    std::string fromBundlePath;  // --from-bundle: the world of this bundle instead of a seed or imported terrain.
    std::string stepsText;       // --steps a..b, or --only p as p..p; empty runs the whole recipe.
    std::vector<std::string> checkpointAt;   // --checkpoint-at: steps (labels or process ids) after which a bundle is written.
    std::vector<double> checkpointTimesMyr;  // --checkpoint-at with numbers: times inside the tectonics evolver.
    double checkpointEveryMyr = 0;           // --checkpoint-every: interval inside the tectonics evolver.
    std::vector<EvolverTimeRequest> evolverDurations;         // --duration <evolver>=<Myr>
    std::vector<EvolverTimeRequest> evolverCheckpointEvery;   // --checkpoint-every <evolver>=<Myr>
    std::vector<EvolverTimeRequest> evolverCheckpointTimes;   // --checkpoint-at <evolver>@<Myr>
    tectonics::scenario::Scenario tectonicScenario = tectonics::scenario::baseline();   // with the --<scenario key> overrides
    bool tectonicScenarioGiven = false;
    bool mineralHistory = true;              // --mineral-history on|off
    std::string checkpointDirectory;         // --checkpoint-dir; default <bundle>-checkpoints.
    std::string accessReportPath;
    std::vector<std::string> commandArguments;
    std::string importLandPath;
    std::string importSeaPath;
    std::string importLakesPath;   // --import-lakes: a lake mask for hydrosphere.imported_lake (imported recipes)
};

void printcommandlineusage();
bool parsecommandlineoptions(const std::vector<std::string>& arguments, CommandLineGenerationOptions& options);
int runcommandlineworldgeneration(const CommandLineGenerationOptions& options);
SocialGenerationOptions makesocialoptions(bool enabled, SocialGenerationOptions::Mode mode, bool useprehistory, int historyyears);
struct WorldRecipe;
struct RecipeRunResult;
struct RecipeRunOptions;
struct GenerationExecutionContext;
// Runs a recipe on an imported world (imported-world or earth-climate).
RecipeRunResult runimportedworldrecipe(planet& world, const WorldRecipe& recipe, const GenerationExecutionContext& context,
    const RecipeRunOptions& runoptions);
void completeimportedworldgeneration(planet& world, bool dodeltas,
    const SocialGenerationOptions& socialoptions, boolshapetemplate smalllake[], boolshapetemplate largelake[],
    boolshapetemplate landshape[], std::vector<std::vector<bool>>& okmountains,
    const ImportedClimateMaps* importedclimate = nullptr, bool climateOnly = false);
