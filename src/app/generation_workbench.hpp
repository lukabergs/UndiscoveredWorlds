#pragma once

#include <deque>
#include <string>
#include <vector>

#include "imgui.h"

#include "app/rendering/appearance_settings.hpp"
#include "world/functions.hpp"
#include "app/rendering/map_views.hpp"
#include "pipeline/generation_stage_registry.hpp"
#include "io/map_imports.hpp"
#include "pipeline/world_stages.hpp"

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
    double tectonicDurationMyr = 0;
    int tectonicPlateCount = 0;
    float tectonicDeltaTimeMyr = 0.0f;
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

void resetworkbenchsession(GenerationSessionState& session);
void initializeproceduralworkbenchsession(GenerationSessionState& session, const planet& world, double durationMyr, int plateCount);
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
// Stage parameters edited in the workbench (samples the hypsometric curve).
WorldStageParameters makestageparameters(const GenerationWorkbenchUiState& ui);
// Runs the procedural-world recipe through the recipe runner the command line uses.
bool runproceduralworldrecipe(planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui, const GenerationExecutionContext& context, std::string* errormessage = nullptr);
GenerationWorkbenchPanelResult drawgenerationworkbenchpanel(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, const planet& displayedworld, GenerationSessionState& session);
