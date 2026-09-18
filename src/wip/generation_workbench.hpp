#pragma once

#include <deque>
#include <string>
#include <vector>

#include "imgui.h"

#include "appearance_settings.hpp"
#include "functions.hpp"
#include "generation_stage_registry.hpp"
#include "map_imports.hpp"

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

struct GenerationScratch
{
    std::vector<std::vector<int>> mountaindrainage;
    std::vector<std::vector<bool>> shelves;
    std::vector<int> squareroot;
    std::vector<std::vector<std::vector<int>>> saltLakeMap;
    std::vector<std::vector<int>> noLake;
    std::vector<std::vector<int>> basinSeeds;
    std::vector<std::vector<int>> inlandSeaComponentIds;
    std::vector<GenerationWaterComponent> inlandSeaComponents;
    std::vector<std::vector<int>> basinComponentIds;
    std::vector<GenerationBasinComponent> basinComponents;
    std::vector<int> heightHistogram;
    std::vector<float> heightCdf;
    std::vector<std::vector<std::uint8_t>> stageFeatureMask;
    int stageFeatureCellCount = 0;
};

struct GenerationCommittedSnapshot
{
    std::size_t stageIndex = 0;
    planet world;
    GenerationScratch scratch;
};

enum class WorkbenchTerrainRenderPreset
{
    terraform,
    grayscale_heightmap,
    custom
};

enum class HypsometricCurveSource
{
    control_points,
    formula
};

struct GenerationWorkbenchUiState
{
    bool previewEnabled = true;
    bool showOceanLandContour = false;
    bool showPlateBoundaries = false;
    int tectonicSeed = 0;
    int tectonicCycleCount = 2;
    int tectonicCycleStepLimit = 600;
    int tectonicPlateCount = 10;
    bool tectonicUseSeaLevelMeters = false;
    int tectonicSeaLevelMeters = 31043;
    int tectonicAggregationOverlapAbsolute = -1;
    float tectonicAggregationOverlapRelative = 0.20f;
    float tectonicFoldingRatio = 0.08f;
    int tectonicErosionPeriod = 60;
    float tectonicErosionStrength = 1.0f;
    float tectonicLandmassRotation = 0.20f;
    float tectonicRotationStrength = 1.0f;
    float tectonicSubductionStrength = 1.0f;
    float tectonicDivergentCarveStrength = 0.015f;
    float tectonicDeltaTimeMyr = 1.0f;
    WorkbenchTerrainRenderPreset terrainRenderPreset = WorkbenchTerrainRenderPreset::terraform;
    MapGradientSettings terrainMapGradient;
    int terrainGradientSelectedStop = 0;
    bool showFeatureOverlay = true;
    float featureOverlayOpacity = 0.55f;
    int seaLevel = 0;
    bool useHypsometricRemap = false;
    HypsometricCurveSource hypsometricCurveSource = HypsometricCurveSource::control_points;
    std::vector<ImVec2> hypsometricControlPoints =
    {
        ImVec2(0.0f, 0.0f),
        ImVec2(0.22f, 0.28f),
        ImVec2(0.78f, 0.72f),
        ImVec2(1.0f, 1.0f)
    };
    std::vector<float> hypsometricControlWeights = { 1.0f, 1.0f, 1.0f, 1.0f };
    int hypsometricSelectedPoint = 1;
    std::string hypsometricFormula = "x";
    float boundaryOverlayWidth = 1.0f;
};

struct GenerationSessionState
{
    bool active = false;
    bool procedural = false;
    bool previewAvailable = false;
    bool importedClimateAvailable = false;
    std::size_t currentStageIndex = 0;
    planet previewWorld;
    GenerationScratch committedScratch;
    GenerationScratch previewScratch;
    std::deque<GenerationCommittedSnapshot> history;
    GenerationWorkbenchUiState ui;
};

struct GenerationExecutionContext
{
    bool dorivers = true;
    bool dolakes = true;
    bool dodeltas = true;
    bool appendclimateworkbook = false;
    SocialGenerationOptions socialoptions;
    boolshapetemplate* landshape = nullptr;
    boolshapetemplate* chainland = nullptr;
    boolshapetemplate* smalllake = nullptr;
    boolshapetemplate* largelake = nullptr;
    std::vector<std::vector<bool>>* okmountains = nullptr;
    const ImportedClimateMaps* importedClimate = nullptr;
};

struct GenerationWorkbenchPanelResult
{
    bool previewModeChanged = false;
    bool displayChanged = false;
    bool discardPreviewRequested = false;
    bool recomputeRequested = false;
    bool applyRequested = false;
    bool skipRequested = false;
    bool resetPreviewRequested = false;
    bool backRequested = false;
    bool abortRequested = false;
    bool controlsChanged = false;
};

void initializegenerationscratch(GenerationScratch& scratch);
void resetworkbenchsession(GenerationSessionState& session);
void initializeproceduralworkbenchsession(GenerationSessionState& session, const planet& world, int plateCycles, int plateCount);
void initializeimportedworkbenchsession(GenerationSessionState& session, const planet& world, bool importedClimateAvailable);
void clearworkbenchpreview(GenerationSessionState& session);
const planet& getworkbenchdisplayworld(const GenerationSessionState& session, const planet& committedworld);
const GenerationStageDefinition& getcurrentgenerationstage(const GenerationSessionState& session);
bool workbenchfinished(const GenerationSessionState& session);
bool canstepbackworkbenchsession(const GenerationSessionState& session);
bool stepbackworkbenchsession(GenerationSessionState& session, planet& committedworld);
mapviewenum preferredworkbenchmapview(GenerationStageId stageId);
bool previewcurrentgenerationstage(GenerationSessionState& session, const planet& committedworld, const GenerationExecutionContext& context, std::string* errormessage = nullptr);
bool applycurrentgenerationstage(GenerationSessionState& session, planet& committedworld, const GenerationExecutionContext& context, std::string* errormessage = nullptr);
void skipcurrentgenerationstage(GenerationSessionState& session, const planet& committedworld);
bool cycleinlandseacomponentpolicyat(GenerationSessionState& session, int x, int y);
bool runcurrentgenerationstagewithoutpreview(GenerationSessionState& session, planet& committedworld, const GenerationExecutionContext& context, std::string* errormessage = nullptr);
bool runremaininggenerationstages(planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui, const GenerationExecutionContext& context, std::size_t startStageIndex, std::string* errormessage = nullptr);
void finalizegeneratedworld(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context);
GenerationWorkbenchPanelResult drawgenerationworkbenchpanel(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, const planet& displayedworld, GenerationSessionState& session);
