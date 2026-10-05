#pragma once

// World generation stages without user interface. The desktop workbench, the
// headless pipeline and the recipe runner all execute stages through
// executegenerationstage. Presentation state (render presets, previews,
// overlays) stays in the workbench.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "pipeline/generation_stage_registry.hpp"
#include "world/planet.hpp"
#include "simulations/society/social_generation.hpp"
#include "core/world/transient_state.hpp"
#include "simulations/geology/world_tectonics.hpp"

class boolshapetemplate;
struct ImportedClimateMaps;

enum class GenerationComponentPolicy
{
    keep,
    fill,
    drain
};

struct GenerationWaterComponent
{
    int id = 0;
    int seedX = 0;
    int seedY = 0;
    int cellCount = 0;
    bool touchesEdge = false;
    GenerationComponentPolicy policy = GenerationComponentPolicy::keep;
};

struct GenerationBasinComponent
{
    int id = 0;
    int seedX = 0;
    int seedY = 0;
    int cellCount = 0;
    BasinClass basinClass = BasinClass::none;
    GenerationComponentPolicy policy = GenerationComponentPolicy::keep;
};

// State carried between stages outside the world. Migration step 5 turns the
// physical parts into world fields; the component lists hold editor policies.
struct GenerationScratch
{
    std::vector<std::vector<int>> mountaindrainage;
    std::vector<int> squareroot;
    std::vector<std::vector<std::vector<float>>> saltLakeMap;
    std::vector<std::vector<int>> noLake;
    std::vector<std::vector<int>> basinSeeds;
    std::vector<std::vector<int>> inlandSeaComponentIds;
    std::vector<GenerationWaterComponent> inlandSeaComponents;
    std::vector<std::vector<int>> basinComponentIds;
    std::vector<GenerationBasinComponent> basinComponents;
    std::vector<int> heightHistogram;
    std::vector<float> heightCdf;
    std::vector<std::vector<std::uint8_t>> stageFeatureMask;   // overlay of cells changed by the last stage
    int stageFeatureCellCount = 0;
    transientstate<GeologySolverState> geologySolver;           // initializer -> tectonics, this run only
};

struct WorldCurvePoint
{
    float x = 0.0f;
    float y = 0.0f;
};

// Simulation inputs edited in the workbench; no presentation state.
struct WorldStageParameters
{
    double tectonicDurationMyr = 0;
    int tectonicPlateCount = 0;
    float tectonicDeltaTimeMyr = 0.0f;
    std::shared_ptr<const tectonics::scenario::Scenario> tectonicScenario;   // null: configs/tectonics-baseline.json
    bool tectonicWorldGrid = false;             // the engine runs on the world grid (tectonics recipe)
    bool tectonicMineralHistory = true;
    double tectonicCheckpointEveryMyr = 0;      // checkpoints inside the tectonics evolver, see evolverCheckpoint
    std::vector<double> tectonicCheckpointsMyr;
    int seaLevel = 0;                           // terraforming target; tectonics writes the result back
    bool useHypsometricRemap = false;
    std::vector<WorldCurvePoint> hypsometricCurve;  // monotone target curve, sampled by the editor
};

struct GenerationExecutionContext
{
    bool dorivers = true;
    bool dolakes = true;
    bool dodeltas = true;
    SocialGenerationOptions socialoptions;
    boolshapetemplate* landshape = nullptr;
    boolshapetemplate* chainland = nullptr;
    boolshapetemplate* smalllake = nullptr;
    boolshapetemplate* largelake = nullptr;
    std::vector<std::vector<bool>>* okmountains = nullptr;
    const ImportedClimateMaps* importedClimate = nullptr;
    bool climateOnly = false; // Recompute climate on an externally supplied, immutable terrain.
    // Called at each checkpoint of an evolver with the world as of that time (Myr).
    // An exception stops the step.
    std::function<void(GenerationStageId stage, double timeMyr)> evolverCheckpoint;
};

inline constexpr int kHeightHistogramBins = 64;

void initializegenerationscratch(GenerationScratch& scratch);
void rebuildheightdistribution(const planet& world, GenerationScratch& scratch);
bool stageusesmanualpreview(GenerationStageId stageId);
void detectinlandseas(planet& world, GenerationScratch& scratch);
void detectbasincomponents(planet& world, GenerationScratch& scratch);

// Runs one stage on the world. Returns false only when the stage cannot run.
bool executegenerationstage(GenerationStageId stageId, planet& world, GenerationScratch& scratch,
    WorldStageParameters& parameters, const GenerationExecutionContext& context, std::string* errormessage = nullptr);

// Physical layers, optional society and the optional climate benchmark record.
void finalizegeneratedworld(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context);
