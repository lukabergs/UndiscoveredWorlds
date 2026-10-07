#include <algorithm>
#include <array>
#include <cerrno>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <stdexcept>
#include <string>

#include "app/generation_workbench.hpp"

#include "app/rendering/map_appearance.hpp"
#include "pipeline/world_recipes.hpp"
#include "pipeline/physical_layers.hpp"
#include "app/ui_charts.hpp"

using namespace std;

namespace
{
constexpr int histogrambincount = 64;

void settooltipifhovered(const char* text)
{
    if (ImGui::IsItemHovered())
        ImGui::SetTooltip("%s", text);
}

struct ScopedWorkbenchDebugRun
{
    bool started = false;
    WorldGenerationDebugOptions options;

    ScopedWorkbenchDebugRun(long seed, double durationMyr, int platecount)
    {
        if (isworldgendebugrunactive())
            return;

        options.visualizePlateTectonicsRealtime = false;
        options.plateTectonicsDurationMyr = max(0.0, durationMyr);
        options.plateTectonicsPlateCount = max(0, platecount);
        beginworldgendebugrun(seed, &options);
        started = true;
    }

    ~ScopedWorkbenchDebugRun()
    {
        if (started)
            endworldgendebugrun();
    }
};

float clamp01(float value)
{
    return clamp(value, 0.0f, 1.0f);
}

struct HypsometricFormulaSample
{
    vector<ImVec2> points;
    string error;
};

class HypsometricFormulaParser
{
public:
    HypsometricFormulaParser(const string& source, float xvalue)
        : text(source), x(xvalue)
    {
    }

    bool parse(float& result, string& error)
    {
        position = 0;
        haserror = false;
        result = parseexpression();
        skipwhitespace();

        if (!haserror && position != text.size())
            seterror("Unexpected trailing text.");

        if (haserror)
        {
            error = errormessage;
            return false;
        }

        return true;
    }

private:
    const string& text;
    size_t position = 0;
    float x = 0.0f;
    bool haserror = false;
    string errormessage;

    void skipwhitespace()
    {
        while (position < text.size() && isspace(static_cast<unsigned char>(text[position])))
            position++;
    }

    void seterror(const char* message)
    {
        if (!haserror)
        {
            haserror = true;
            errormessage = message;
        }
    }

    bool consume(char token)
    {
        skipwhitespace();

        if (position < text.size() && text[position] == token)
        {
            position++;
            return true;
        }

        return false;
    }

    float parseexpression()
    {
        float value = parseterm();

        while (!haserror)
        {
            if (consume('+'))
                value += parseterm();
            else if (consume('-'))
                value -= parseterm();
            else
                break;
        }

        return value;
    }

    float parseterm()
    {
        float value = parsepower();

        while (!haserror)
        {
            if (consume('*'))
                value *= parsepower();
            else if (consume('/'))
            {
                const float divisor = parsepower();

                if (fabsf(divisor) <= 0.000001f)
                    seterror("Division by zero.");
                else
                    value /= divisor;
            }
            else
                break;
        }

        return value;
    }

    float parsepower()
    {
        float value = parseunary();

        if (consume('^'))
        {
            const float exponent = parsepower();
            value = powf(value, exponent);
        }

        return value;
    }

    float parseunary()
    {
        if (consume('+'))
            return parseunary();

        if (consume('-'))
            return -parseunary();

        return parseprimary();
    }

    float parseprimary()
    {
        skipwhitespace();

        if (position >= text.size())
        {
            seterror("Unexpected end of formula.");
            return 0.0f;
        }

        if (consume('('))
        {
            const float value = parseexpression();

            if (!consume(')'))
                seterror("Missing closing ')'.");

            return value;
        }

        if (isdigit(static_cast<unsigned char>(text[position])) || text[position] == '.')
            return parsenumber();

        if (isalpha(static_cast<unsigned char>(text[position])))
            return parseidentifier();

        seterror("Unexpected token.");
        return 0.0f;
    }

    float parsenumber()
    {
        skipwhitespace();
        char* end = nullptr;
        errno = 0;
        const double value = strtod(text.c_str() + position, &end);

        if (end == text.c_str() + position || errno == ERANGE)
        {
            seterror("Invalid number.");
            return 0.0f;
        }

        position = static_cast<size_t>(end - text.c_str());
        return static_cast<float>(value);
    }

    string parseidentifiername()
    {
        skipwhitespace();
        const size_t start = position;

        while (position < text.size() &&
            (isalnum(static_cast<unsigned char>(text[position])) || text[position] == '_'))
        {
            position++;
        }

        return text.substr(start, position - start);
    }

    float parseidentifier()
    {
        const string name = parseidentifiername();

        if (name == "x")
            return x;

        if (name == "pi")
            return 3.14159265358979323846f;

        if (name == "e")
            return 2.71828182845904523536f;

        if (!consume('('))
        {
            seterror("Unknown identifier.");
            return 0.0f;
        }

        vector<float> arguments;
        arguments.push_back(parseexpression());

        while (!haserror && consume(','))
            arguments.push_back(parseexpression());

        if (!consume(')'))
            seterror("Missing closing ')' after function call.");

        if (haserror)
            return 0.0f;

        return evaluatename(name, arguments);
    }

    float evaluatename(const string& name, const vector<float>& arguments)
    {
        if (name == "abs" && arguments.size() == 1)
            return fabsf(arguments[0]);
        if (name == "sqrt" && arguments.size() == 1)
            return sqrtf(max(0.0f, arguments[0]));
        if (name == "exp" && arguments.size() == 1)
            return expf(arguments[0]);
        if (name == "log" && arguments.size() == 1)
            return logf(max(arguments[0], 0.000001f));
        if (name == "log10" && arguments.size() == 1)
            return log10f(max(arguments[0], 0.000001f));
        if (name == "sin" && arguments.size() == 1)
            return sinf(arguments[0]);
        if (name == "cos" && arguments.size() == 1)
            return cosf(arguments[0]);
        if (name == "tan" && arguments.size() == 1)
            return tanf(arguments[0]);
        if (name == "floor" && arguments.size() == 1)
            return floorf(arguments[0]);
        if (name == "ceil" && arguments.size() == 1)
            return ceilf(arguments[0]);
        if (name == "pow" && arguments.size() == 2)
            return powf(arguments[0], arguments[1]);
        if (name == "min" && arguments.size() == 2)
            return min(arguments[0], arguments[1]);
        if (name == "max" && arguments.size() == 2)
            return max(arguments[0], arguments[1]);
        if (name == "clamp" && arguments.size() == 3)
            return clamp(arguments[0], arguments[1], arguments[2]);

        seterror("Unknown function or wrong argument count.");
        return 0.0f;
    }
};

bool evaluatehypsometricformula(const string& formula, float x, float& value, string& error)
{
    HypsometricFormulaParser parser(formula, x);
    return parser.parse(value, error);
}

float normalizeterrainvalue(int current, int minimum, int maximum)
{
    return clamp01(static_cast<float>(current - minimum) / static_cast<float>(max(1, maximum - minimum)));
}

float clampcurveweight(float value)
{
    return clamp(value, 0.1f, 8.0f);
}

int hypsometriccurvedegree(size_t controlpointcount)
{
    if (controlpointcount <= 1)
        return 0;

    return min(3, static_cast<int>(controlpointcount) - 1);
}

vector<float> buildhypsometricknots(size_t controlpointcount)
{
    const int degree = hypsometriccurvedegree(controlpointcount);
    const size_t knotcount = controlpointcount + static_cast<size_t>(degree) + 1;
    vector<float> knots(knotcount, 0.0f);

    if (controlpointcount == 0)
        return knots;

    const size_t interiorcount = knotcount - static_cast<size_t>((degree + 1) * 2);

    for (size_t index = 0; index < knotcount; index++)
    {
        if (index <= static_cast<size_t>(degree))
            knots[index] = 0.0f;
        else if (index >= controlpointcount)
            knots[index] = 1.0f;
        else
            knots[index] = static_cast<float>(index - degree) / static_cast<float>(interiorcount + 1);
    }

    return knots;
}

float evaluatebsplinebasis(int index, int degree, float t, const vector<float>& knots)
{
    if (degree == 0)
        return (t >= knots[index] && t < knots[index + 1]) || (t >= 1.0f && knots[index + 1] >= 1.0f) ? 1.0f : 0.0f;

    float left = 0.0f;
    float right = 0.0f;
    const float leftdenominator = knots[index + degree] - knots[index];
    const float rightdenominator = knots[index + degree + 1] - knots[index + 1];

    if (leftdenominator > 0.0001f)
        left = ((t - knots[index]) / leftdenominator) * evaluatebsplinebasis(index, degree - 1, t, knots);

    if (rightdenominator > 0.0001f)
        right = ((knots[index + degree + 1] - t) / rightdenominator) * evaluatebsplinebasis(index + 1, degree - 1, t, knots);

    return left + right;
}

ImVec2 evaluatehypsometriccurve(const vector<ImVec2>& controlpoints, const vector<float>& weights, float t)
{
    if (controlpoints.empty())
        return ImVec2(0.0f, 0.0f);

    const int degree = hypsometriccurvedegree(controlpoints.size());
    const vector<float> knots = buildhypsometricknots(controlpoints.size());
    const float clampedt = clamp01(t);
    float totalweight = 0.0f;
    ImVec2 point(0.0f, 0.0f);

    for (size_t index = 0; index < controlpoints.size(); index++)
    {
        const float basis = evaluatebsplinebasis(static_cast<int>(index), degree, clampedt, knots);
        const float weight = clampcurveweight(index < weights.size() ? weights[index] : 1.0f);
        const float scaledbasis = basis * weight;
        totalweight += scaledbasis;
        point.x += controlpoints[index].x * scaledbasis;
        point.y += controlpoints[index].y * scaledbasis;
    }

    if (totalweight <= 0.0001f)
        return ImVec2(clamp01(controlpoints.back().x), clamp01(controlpoints.back().y));

    return ImVec2(clamp01(point.x / totalweight), clamp01(point.y / totalweight));
}

HypsometricFormulaSample samplehypsometricformulacurve(const GenerationWorkbenchUiState& ui)
{
    HypsometricFormulaSample sample;
    sample.points.reserve(histogrambincount);
    float previousy = 0.0f;

    for (int index = 0; index < histogrambincount; index++)
    {
        const float x = histogrambincount > 1 ? static_cast<float>(index) / static_cast<float>(histogrambincount - 1) : 0.0f;
        float y = x;

        if (!evaluatehypsometricformula(ui.hypsometricFormula, x, y, sample.error))
        {
            sample.points.clear();
            break;
        }

        if (!std::isfinite(y))
        {
            sample.error = "Formula produced a non-finite value.";
            sample.points.clear();
            break;
        }

        y = clamp01(y);

        if (index == 0)
            y = 0.0f;
        else if (index == histogrambincount - 1)
            y = 1.0f;

        y = max(previousy, y);
        previousy = y;
        sample.points.push_back(ImVec2(x, y));
    }

    if (sample.points.empty())
    {
        sample.points.reserve(histogrambincount);

        for (int index = 0; index < histogrambincount; index++)
        {
            const float x = histogrambincount > 1 ? static_cast<float>(index) / static_cast<float>(histogrambincount - 1) : 0.0f;
            sample.points.push_back(ImVec2(x, x));
        }
    }

    return sample;
}

vector<ImVec2> samplehypsometriccurve(const GenerationWorkbenchUiState& ui)
{
    if (ui.hypsometricCurveSource == HypsometricCurveSource::formula)
        return samplehypsometricformulacurve(ui).points;

    vector<ImVec2> points;
    points.reserve(histogrambincount);

    for (int index = 0; index < histogrambincount; index++)
    {
        const float t = histogrambincount > 1 ? static_cast<float>(index) / static_cast<float>(histogrambincount - 1) : 0.0f;
        points.push_back(evaluatehypsometriccurve(ui.hypsometricControlPoints, ui.hypsometricControlWeights, t));
    }

    return points;
}

pair<int, int> terrainrendergradientrange(const planet& world)
{
    return { 1 - world.sealevel(), (world.maxelevation() - 1) - world.sealevel() };
}

const char* terrainrenderpresetlabel(WorkbenchTerrainRenderPreset preset)
{
    switch (preset)
    {
    case WorkbenchTerrainRenderPreset::terraform:
        return "Terraform";
    case WorkbenchTerrainRenderPreset::grayscale_heightmap:
        return "Grayscale heightmap";
    case WorkbenchTerrainRenderPreset::custom:
    default:
        return "Custom";
    }
}

void applyterrainrenderpreset(GenerationWorkbenchUiState& ui, const planet& world, WorkbenchTerrainRenderPreset preset)
{
    const auto [minimum, maximum] = terrainrendergradientrange(world);
    MapGradientSettings gradient;
    gradient.discrete = false;

    switch (preset)
    {
    case WorkbenchTerrainRenderPreset::grayscale_heightmap:
        gradient.stopcount = 2;
        gradient.stops[0].position = minimum;
        gradient.stops[0].colour = ImVec4(0.02f, 0.02f, 0.02f, 1.0f);
        gradient.stops[1].position = maximum;
        gradient.stops[1].colour = ImVec4(0.98f, 0.98f, 0.98f, 1.0f);
        break;

    case WorkbenchTerrainRenderPreset::terraform:
    default:
    {
        const int shoreline = minimum < 0
            ? clamp(-max(1, static_cast<int>(roundf(static_cast<float>(abs(minimum)) * 0.08f))), minimum, -1)
            : minimum;
        const int quarterland = max(0, static_cast<int>(roundf(static_cast<float>(maximum) * 0.25f)));
        const int halfland = max(quarterland, static_cast<int>(roundf(static_cast<float>(maximum) * 0.50f)));
        const int threequarterland = max(halfland, static_cast<int>(roundf(static_cast<float>(maximum) * 0.75f)));

        const int sealevelstop = clamp(0, minimum, maximum);
        gradient.stopcount = 7;
        gradient.stops[0] = { minimum, ImVec4(0.11f, 0.16f, 0.32f, 1.0f) };
        gradient.stops[1] = { shoreline, ImVec4(1.0f, 1.0f, 1.0f, 1.0f) };
        gradient.stops[2] = { sealevelstop, ImVec4(0.02f, 0.39f, 0.03f, 1.0f) };
        gradient.stops[3] = { max(1, quarterland), ImVec4(0.98f, 0.93f, 0.36f, 1.0f) };
        gradient.stops[4] = { max(gradient.stops[3].position + 1, halfland), ImVec4(1.0f, 0.75f, 0.0f, 1.0f) };
        gradient.stops[5] = { max(gradient.stops[4].position + 1, threequarterland), ImVec4(0.39f, 0.35f, 0.13f, 1.0f) };
        gradient.stops[6] = { max(gradient.stops[5].position + 1, maximum), ImVec4(0.51f, 0.0f, 0.50f, 1.0f) };
        break;
    }
    }

    normalizegradientstops(gradient);
    ui.terrainRenderPreset = preset;
    ui.terrainMapGradient = gradient;
    ui.terrainGradientSelectedStop = clamp(ui.terrainGradientSelectedStop, 0, max(0, gradient.stopcount - 1));
}



vector<float> histogramasfloats(const vector<int>& values)
{
    vector<float> result(values.size(), 0.0f);

    for (size_t index = 0; index < values.size(); index++)
        result[index] = static_cast<float>(values[index]);

    return result;
}

vector<ImVec2> cdfcurvepoints(const vector<float>& values)
{
    vector<ImVec2> points;
    if (values.empty())
        return points;

    points.reserve(values.size());

    for (size_t index = 0; index < values.size(); index++)
    {
        const float height = values.size() > 1 ? static_cast<float>(index) / static_cast<float>(values.size() - 1) : 0.0f;
        points.emplace_back(clamp01(values[index]), height);
    }

    return points;
}
}


void pushworkbenchsnapshot(GenerationSessionState& session, const planet& committedworld)
{
    session.history.emplace_back();
    GenerationCommittedSnapshot& snapshot = session.history.back();
    snapshot.stageIndex = session.currentStageIndex;
    snapshot.world = committedworld;
    snapshot.scratch = session.committedScratch;
}

void resetworkbenchsession(GenerationSessionState& session)
{
    session.active = false;
    session.procedural = false;
    session.previewAvailable = false;
    session.importedClimateAvailable = false;
    session.currentStageIndex = 0;
    session.ui = GenerationWorkbenchUiState{};
    initializegenerationscratch(session.committedScratch);
    initializegenerationscratch(session.previewScratch);
    session.history.clear();
}

void initializeproceduralworkbenchsession(GenerationSessionState& session, const planet& world, double durationMyr, int plateCount)
{
    resetworkbenchsession(session);
    session.active = true;
    session.procedural = true;
    session.previewAvailable = false;
    session.currentStageIndex = getgenerationstageindex(GenerationStageId::geosphere_initialization);
    session.previewWorld = world;
    initializegenerationscratch(session.committedScratch);
    initializegenerationscratch(session.previewScratch);
    session.ui.tectonicSeed = static_cast<int>(world.seed());
    session.ui.tectonicDurationMyr = durationMyr;
    session.ui.tectonicPlateCount = plateCount;
    session.ui.seaLevel = world.sealevel();
    applyterrainrenderpreset(session.ui, world, session.ui.terrainRenderPreset);
    rebuildheightdistribution(world, session.committedScratch);
    pushworkbenchsnapshot(session, world);
}

void initializeimportedworkbenchsession(GenerationSessionState& session, const planet& world, bool importedClimateAvailable)
{
    resetworkbenchsession(session);
    session.active = true;
    session.procedural = false;
    session.importedClimateAvailable = importedClimateAvailable;
    session.previewAvailable = false;
    session.currentStageIndex = getgenerationstageindex(GenerationStageId::mountain_bases);
    session.previewWorld = world;
    session.previewWorld.setmaxelevation(200000);
    initializegenerationscratch(session.committedScratch);
    initializegenerationscratch(session.previewScratch);
    session.ui.tectonicSeed = static_cast<int>(world.seed());
    session.ui.seaLevel = world.sealevel();
    applyterrainrenderpreset(session.ui, world, session.ui.terrainRenderPreset);
    rebuildheightdistribution(world, session.committedScratch);
    pushworkbenchsnapshot(session, world);
}

void clearworkbenchpreview(GenerationSessionState& session)
{
    session.previewAvailable = false;
}

WorldStageParameters makestageparameters(const GenerationWorkbenchUiState& ui)
{
    WorldStageParameters parameters;
    parameters.tectonicDurationMyr = ui.tectonicDurationMyr;
    parameters.tectonicPlateCount = ui.tectonicPlateCount;
    parameters.tectonicDeltaTimeMyr = ui.tectonicDeltaTimeMyr;
    parameters.seaLevel = ui.seaLevel;
    parameters.useHypsometricRemap = ui.useHypsometricRemap;
    if (ui.useHypsometricRemap)
        for (const ImVec2& point : samplehypsometriccurve(ui))
            parameters.hypsometricCurve.push_back({ point.x, point.y });
    return parameters;
}

namespace
{
// Runs a stage with the workbench's parameters; render presets stay UI state.
bool executeworkbenchstage(GenerationStageId stageid, planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui,
    const GenerationExecutionContext& context, string* errormessage)
{
    WorldStageParameters parameters = makestageparameters(ui);
    // The rivers come from a fresh routing after the fjords, as in imported-world: the workbench
    // steps over hydrology.flow_routing at its registry place, before the fjords.
    if (stageid == GenerationStageId::rivers
        && !executegenerationstage(GenerationStageId::flow_routing, world, scratch, parameters, context, errormessage))
        return false;
    const bool executed = executegenerationstage(stageid, world, scratch, parameters, context, errormessage);
    ui.seaLevel = parameters.seaLevel;
    if (executed && stageid == GenerationStageId::terraforming && ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
        applyterrainrenderpreset(ui, world, ui.terrainRenderPreset);
    return executed;
}
}

const planet& getworkbenchdisplayworld(const GenerationSessionState& session, const planet& committedworld)
{
    if (session.active && session.ui.previewEnabled && session.previewAvailable)
        return session.previewWorld;

    return committedworld;
}

const GenerationStageDefinition& getcurrentgenerationstage(const GenerationSessionState& session)
{
    const vector<GenerationStageDefinition>& stages = getgenerationstages();
    const size_t index = min(session.currentStageIndex, stages.size() - 1);
    return stages[index];
}

bool workbenchfinished(const GenerationSessionState& session)
{
    return session.currentStageIndex >= getgenerationstages().size();
}

bool canstepbackworkbenchsession(const GenerationSessionState& session)
{
    return session.history.size() > 1;
}

bool stepbackworkbenchsession(GenerationSessionState& session, planet& committedworld)
{
    if (canstepbackworkbenchsession(session) == false)
        return false;

    session.history.pop_back();
    const GenerationCommittedSnapshot& snapshot = session.history.back();
    committedworld = snapshot.world;
    session.committedScratch = snapshot.scratch;
    session.currentStageIndex = snapshot.stageIndex;
    session.previewAvailable = false;
    session.previewWorld = committedworld;
    session.ui.seaLevel = committedworld.sealevel();
    if (session.ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
        applyterrainrenderpreset(session.ui, committedworld, session.ui.terrainRenderPreset);
    return true;
}

mapviewenum preferredworkbenchmapview(GenerationStageId stageId)
{
    switch (stageId)
    {
    case GenerationStageId::geosphere_initialization:
    case GenerationStageId::plate_tectonics:
    case GenerationStageId::tectonic_trenches:
    case GenerationStageId::tectonic_volcanoes:
    case GenerationStageId::flow_routing:
    case GenerationStageId::fluvial_erosion:
        return relief;

    case GenerationStageId::terraforming:
        return elevation;

    case GenerationStageId::global_temperature:
    case GenerationStageId::sea_surface_temperatures:
    case GenerationStageId::pressure:
    case GenerationStageId::winds:
    case GenerationStageId::mountain_temperature_lapse:
        return temperature;

    case GenerationStageId::rainfall:
    case GenerationStageId::arid_features:
        return precipitation;

    case GenerationStageId::rivers:
    case GenerationStageId::lake_levels:
    case GenerationStageId::deltas_wetlands_and_roughness:
        return rivers;

    case GenerationStageId::climates_and_biomes:
        return climate;

    default:
        return relief;
    }
}


namespace
{
// The step-by-step workbench walks the registered stages with one climate pass, so it runs the
// hydrology of imported-world, not the erosion (procedural-world runs it through
// runproceduralworldrecipe): it steps over terrain.fluvial_erosion and over hydrology.flow_routing
// at its registry place (the rivers step routes first, after the fjords; erosion plan H4).
bool workbenchskipsstage(GenerationStageId stage)
{
    return stage == GenerationStageId::flow_routing || stage == GenerationStageId::fluvial_erosion;
}

void advanceworkbenchstage(GenerationSessionState& session)
{
    const vector<GenerationStageDefinition>& stages = getgenerationstages();
    if (session.currentStageIndex < stages.size())
        session.currentStageIndex++;
    while (session.currentStageIndex < stages.size() && workbenchskipsstage(stages[session.currentStageIndex].id))
        session.currentStageIndex++;

    if (workbenchfinished(session) == false && getcurrentgenerationstage(session).id == GenerationStageId::terraforming)
        session.ui.seaLevel = session.previewWorld.sealevel();

    session.previewAvailable = false;
}

void seedcommittedwatercomponents(GenerationSessionState& session)
{
    if (session.committedScratch.inlandSeaComponents.empty() && session.previewAvailable && session.previewScratch.inlandSeaComponents.empty() == false)
    {
        session.committedScratch.inlandSeaComponentIds = session.previewScratch.inlandSeaComponentIds;
        session.committedScratch.inlandSeaComponents = session.previewScratch.inlandSeaComponents;
    }
}

void cyclepolicy(GenerationComponentPolicy& policy)
{
    switch (policy)
    {
    case GenerationComponentPolicy::keep:
        policy = GenerationComponentPolicy::drain;
        break;
    case GenerationComponentPolicy::drain:
        policy = GenerationComponentPolicy::fill;
        break;
    case GenerationComponentPolicy::fill:
        policy = GenerationComponentPolicy::keep;
        break;
    }
}
}

bool previewcurrentgenerationstage(GenerationSessionState& session, const planet& committedworld, const GenerationExecutionContext& context, string* errormessage)
{
    if (session.active == false || workbenchfinished(session))
        return false;

    const GenerationStageId stageid = getcurrentgenerationstage(session).id;
    const bool tectonicstage = stageid == GenerationStageId::geosphere_initialization || stageid == GenerationStageId::plate_tectonics;
    const long debugseed = tectonicstage ? static_cast<long>(session.ui.tectonicSeed) : committedworld.seed();
    ScopedWorkbenchDebugRun debugrun(debugseed, session.ui.tectonicDurationMyr, session.ui.tectonicPlateCount);
    session.previewWorld = committedworld;
    if (tectonicstage)
        session.previewWorld.setseed(debugseed);
    session.previewScratch = session.committedScratch;
    session.previewAvailable = true;

    if (stageusesmanualpreview(stageid))
    {
        if (stageid == GenerationStageId::inland_seas && session.previewScratch.inlandSeaComponents.empty())
            detectinlandseas(session.previewWorld, session.previewScratch);

        if (isterraingenerationstage(stageid) && session.ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
            applyterrainrenderpreset(session.ui, session.previewWorld, session.ui.terrainRenderPreset);

        return true;
    }

    session.previewAvailable = executeworkbenchstage(stageid, session.previewWorld, session.previewScratch, session.ui, context, errormessage);
    if (session.previewAvailable && isterraingenerationstage(stageid) && session.ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
        applyterrainrenderpreset(session.ui, session.previewWorld, session.ui.terrainRenderPreset);
    return session.previewAvailable;
}

bool applycurrentgenerationstage(GenerationSessionState& session, planet& committedworld, const GenerationExecutionContext& context, string* errormessage)
{
    if (session.active == false || workbenchfinished(session))
        return false;

    if (session.previewAvailable == false)
    {
        if (previewcurrentgenerationstage(session, committedworld, context, errormessage) == false)
            return false;
    }

    if (stageusesmanualpreview(getcurrentgenerationstage(session).id))
    {
        if (executeworkbenchstage(getcurrentgenerationstage(session).id, session.previewWorld, session.previewScratch, session.ui, context, errormessage) == false)
            return false;
    }

    committedworld = session.previewWorld;
    session.committedScratch = session.previewScratch;
    advanceworkbenchstage(session);
    pushworkbenchsnapshot(session, committedworld);
    return true;
}

void skipcurrentgenerationstage(GenerationSessionState& session, const planet& committedworld)
{
    if (session.active == false || workbenchfinished(session))
        return;

    advanceworkbenchstage(session);
    pushworkbenchsnapshot(session, committedworld);
}

bool cycleinlandseacomponentpolicyat(GenerationSessionState& session, int x, int y)
{
    if (session.active == false || workbenchfinished(session) || getcurrentgenerationstage(session).id != GenerationStageId::inland_seas)
        return false;

    seedcommittedwatercomponents(session);

    const GenerationScratch& visible = (session.previewAvailable && session.ui.previewEnabled) ? session.previewScratch : session.committedScratch;

    if (x < 0 || y < 0 || x >= static_cast<int>(visible.inlandSeaComponentIds.size()) || y >= static_cast<int>(visible.inlandSeaComponentIds[x].size()))
        return false;

    const int componentid = visible.inlandSeaComponentIds[x][y];
    if (componentid <= 0)
        return false;

    GenerationComponentPolicy nextpolicy = GenerationComponentPolicy::keep;
    bool found = false;

    for (const GenerationWaterComponent& component : visible.inlandSeaComponents)
    {
        if (component.id != componentid)
            continue;

        nextpolicy = component.policy;
        cyclepolicy(nextpolicy);
        found = true;
        break;
    }

    if (!found)
        return false;

    auto applypolicy = [&](GenerationScratch& scratch)
    {
        for (GenerationWaterComponent& component : scratch.inlandSeaComponents)
        {
            if (component.id == componentid)
            {
                component.policy = nextpolicy;
                break;
            }
        }
    };

    applypolicy(session.committedScratch);
    applypolicy(session.previewScratch);
    return true;
}

bool runcurrentgenerationstagewithoutpreview(GenerationSessionState& session, planet& committedworld, const GenerationExecutionContext& context, string* errormessage)
{
    session.previewAvailable = false;
    return applycurrentgenerationstage(session, committedworld, context, errormessage);
}

bool runproceduralworldrecipe(planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui, const GenerationExecutionContext& context, string* errormessage)
{
    // The recipe the command line runs (world_pipeline run --recipe procedural-world): it has two
    // climate passes and repeated hydrology steps, which the registry walk of the workbench cannot
    // express. Each step still reports its stage through beginworldgenstep.
    ScopedWorkbenchDebugRun debugrun(world.seed(), ui.tectonicDurationMyr, ui.tectonicPlateCount);
    const WorldRecipe recipe = proceduralworldrecipe();
    WorldStageParameters parameters = makestageparameters(ui);
    RecipeRunOptions options;
    options.afterstep = [&](const RecipeRunResult&, const RecipeStepRecord& step, string&)
    {
        if (step.stage == GenerationStageId::terraforming)
        {
            ui.seaLevel = parameters.seaLevel;
            if (ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
                applyterrainrenderpreset(ui, world, ui.terrainRenderPreset);
        }
        return true;
    };
    const RecipeRunResult result = runrecipe(recipe, world, scratch, parameters, context, options);
    ui.seaLevel = parameters.seaLevel;
    if (!result.completed && errormessage != nullptr)
        *errormessage = result.error;
    return result.completed;
}


GenerationWorkbenchPanelResult drawgenerationworkbenchpanel(const ImGuiViewport* main_viewport, ImGuiWindowFlags window_flags, const planet& displayedworld, GenerationSessionState& session)
{
    GenerationWorkbenchPanelResult result;

    if (session.active == false || workbenchfinished(session))
        return result;

    const GenerationStageDefinition& stage = getcurrentgenerationstage(session);
    const GenerationScratch& panelscratch = (session.previewAvailable && session.ui.previewEnabled) ? session.previewScratch : session.committedScratch;

    ImGui::SetNextWindowPos(ImVec2(main_viewport->WorkPos.x + 10, main_viewport->WorkPos.y + 20), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(352, 735), ImGuiCond_FirstUseEver);
    ImGui::Begin("Generation Workbench", nullptr, window_flags);

    ImGui::Text("Stage %d / %d", static_cast<int>(session.currentStageIndex + 1), static_cast<int>(getgenerationstages().size()));
    ImGui::Separator();
    ImGui::TextUnformatted(stage.label);
    ImGui::TextWrapped("%s", stage.description);
    if (ImGui::Checkbox("Preview", &session.ui.previewEnabled))
        result.previewModeChanged = true;
    result.displayChanged |= ImGui::Checkbox("Ocean-land contour", &session.ui.showOceanLandContour);
    settooltipifhovered("Overlay the ocean-land contour on the current map preview.");
    result.displayChanged |= ImGui::Checkbox("Show plate boundaries", &session.ui.showPlateBoundaries);
    settooltipifhovered("Overlay tectonic plate boundaries on the current map preview.");

    if (isterraingenerationstage(stage.id))
    {
        const auto [gradientminimum, gradientmaximum] = terrainrendergradientrange(displayedworld);
        if (ImGui::BeginCombo("Map style", terrainrenderpresetlabel(session.ui.terrainRenderPreset)))
        {
            for (WorkbenchTerrainRenderPreset preset : { WorkbenchTerrainRenderPreset::grayscale_heightmap, WorkbenchTerrainRenderPreset::terraform })
            {
                const bool selected = preset == session.ui.terrainRenderPreset;
                if (ImGui::Selectable(terrainrenderpresetlabel(preset), selected))
                {
                    applyterrainrenderpreset(session.ui, displayedworld, preset);
                    result.displayChanged = true;
                }

                if (selected)
                    ImGui::SetItemDefaultFocus();
            }

            ImGui::EndCombo();
        }
        settooltipifhovered("Choose a terrain render preset. Editing the gradient below switches this to Custom.");

        if (drawgradienteditor("WorkbenchTerrainGradient", session.ui.terrainMapGradient, gradientminimum, gradientmaximum, " m", session.ui.terrainGradientSelectedStop))
        {
            session.ui.terrainRenderPreset = WorkbenchTerrainRenderPreset::custom;
            result.displayChanged = true;
        }

        if (panelscratch.stageFeatureCellCount > 0 && stage.id != GenerationStageId::inland_seas)
        {
            result.displayChanged |= ImGui::Checkbox("View features", &session.ui.showFeatureOverlay);
            ImGui::SameLine();
            result.displayChanged |= ImGui::SliderFloat("Opacity", &session.ui.featureOverlayOpacity, 0.0f, 1.0f, "%.2f");
        }
    }

    if (stage.id == GenerationStageId::geosphere_initialization || stage.id == GenerationStageId::plate_tectonics)
    {
        // Both tectonic stages use the same seed and engine parameters.
        result.discardPreviewRequested |= ImGui::InputInt("Seed", &session.ui.tectonicSeed);
        if (session.ui.tectonicSeed < 0)
            session.ui.tectonicSeed = 0;
        settooltipifhovered("Seed used for plate tectonics. The same seed yields the same result.");

        ImGui::SameLine();
        if (ImGui::Button("Recompute"))
            result.recomputeRequested = true;
        settooltipifhovered("Rerun this stage with the current seed and parameters.");

        ImGui::TextWrapped("Uses configs/tectonics-baseline.json and the repository's current column solver.");
        result.discardPreviewRequested |= ImGui::InputDouble("Duration (Myr; 0 = recipe)", &session.ui.tectonicDurationMyr);
        session.ui.tectonicDurationMyr = max(0.0, session.ui.tectonicDurationMyr);
        result.discardPreviewRequested |= ImGui::InputInt("Plate count (0 = recipe)", &session.ui.tectonicPlateCount);
        session.ui.tectonicPlateCount = max(0, session.ui.tectonicPlateCount);
        result.discardPreviewRequested |= ImGui::InputFloat("Timestep (Myr; 0 = recipe)", &session.ui.tectonicDeltaTimeMyr);
        session.ui.tectonicDeltaTimeMyr = max(0.0f, session.ui.tectonicDeltaTimeMyr);

        if (session.previewAvailable == false)
            ImGui::TextWrapped("Press Recompute to preview this stage with the current seed and parameters.");
    }
    else if (stage.id == GenerationStageId::terraforming)
    {
        result.controlsChanged |= ImGui::SliderInt("Sea level", &session.ui.seaLevel, 1, displayedworld.maxelevation() - 2);
        result.controlsChanged |= ImGui::Checkbox("Enable hypsometric remap", &session.ui.useHypsometricRemap);

        if (ImGui::BeginTabBar("TerraformTabs"))
        {
            if (ImGui::BeginTabItem("Histogram"))
            {
                float histogramsealevel = static_cast<float>(session.ui.seaLevel);
                if (drawheighthistogramwidget("Elevation histogram", histogramasfloats(panelscratch.heightHistogram), 1.0f, static_cast<float>(displayedworld.maxelevation() - 1), histogramsealevel, &session.ui.terrainMapGradient, " m", ImVec2(0.0f, 150.0f)))
                {
                    session.ui.seaLevel = clamp(static_cast<int>(roundf(histogramsealevel)), 1, displayedworld.maxelevation() - 2);
                    result.controlsChanged = true;
                }

                ImGui::EndTabItem();
            }

            if (ImGui::BeginTabItem("Curve"))
            {
                vector<ImVec2> curve = cdfcurvepoints(panelscratch.heightCdf);
                const HypsometricFormulaSample formulasample = samplehypsometricformulacurve(session.ui);
                float sealevel = normalizeterrainvalue(session.ui.seaLevel, 1, displayedworld.maxelevation() - 1);

                int curvesource = static_cast<int>(session.ui.hypsometricCurveSource);
                if (ImGui::RadioButton("Control points", &curvesource, static_cast<int>(HypsometricCurveSource::control_points)))
                {
                    session.ui.useHypsometricRemap = true;
                    result.controlsChanged = true;
                }
                ImGui::SameLine();
                if (ImGui::RadioButton("Formula", &curvesource, static_cast<int>(HypsometricCurveSource::formula)))
                {
                    session.ui.useHypsometricRemap = true;
                    result.controlsChanged = true;
                }
                session.ui.hypsometricCurveSource = static_cast<HypsometricCurveSource>(curvesource);

                const bool editablecurve = session.ui.hypsometricCurveSource == HypsometricCurveSource::control_points;
                const vector<ImVec2>* previewcurve = editablecurve ? nullptr : &formulasample.points;

                if (drawhypsometriccurvewidget("Hypsometric curve", curve, session.ui.hypsometricControlPoints, session.ui.hypsometricControlWeights, session.ui.hypsometricSelectedPoint, sealevel, ImVec2(0.0f, 190.0f), previewcurve, editablecurve))
                {
                    session.ui.useHypsometricRemap = true;
                    session.ui.seaLevel = clamp(1 + static_cast<int>(roundf(sealevel * static_cast<float>(displayedworld.maxelevation() - 2))), 1, displayedworld.maxelevation() - 2);
                    result.controlsChanged = true;
                }

                if (editablecurve && session.ui.hypsometricSelectedPoint > 0 && session.ui.hypsometricSelectedPoint < static_cast<int>(session.ui.hypsometricControlPoints.size()) - 1)
                {
                    float& weight = session.ui.hypsometricControlWeights[session.ui.hypsometricSelectedPoint];
                    if (ImGui::InputFloat("Selected weight", &weight, 0.1f, 0.5f, "%.3f"))
                    {
                        weight = clampcurveweight(weight);
                        session.ui.useHypsometricRemap = true;
                        result.controlsChanged = true;
                    }

                    if (session.ui.hypsometricControlPoints.size() > 4)
                    {
                        if (ImGui::Button("Delete selected point"))
                        {
                            session.ui.hypsometricControlPoints.erase(session.ui.hypsometricControlPoints.begin() + session.ui.hypsometricSelectedPoint);
                            session.ui.hypsometricControlWeights.erase(session.ui.hypsometricControlWeights.begin() + session.ui.hypsometricSelectedPoint);
                            session.ui.hypsometricSelectedPoint = clamp(session.ui.hypsometricSelectedPoint - 1, 1, static_cast<int>(session.ui.hypsometricControlPoints.size()) - 2);
                            session.ui.useHypsometricRemap = true;
                            result.controlsChanged = true;
                        }
                    }
                }
                else if (!editablecurve)
                {
                    char formulabuffer[256];
                    memset(formulabuffer, 0, sizeof(formulabuffer));
                    strncpy_s(formulabuffer, sizeof(formulabuffer), session.ui.hypsometricFormula.c_str(), _TRUNCATE);

                    if (ImGui::InputText("y(x)", formulabuffer, sizeof(formulabuffer)))
                    {
                        session.ui.hypsometricFormula = formulabuffer;
                        session.ui.useHypsometricRemap = true;
                        result.controlsChanged = true;
                    }

                    ImGui::TextUnformatted("Use x in [0, 1]. Functions: abs sqrt exp log log10 sin cos tan floor ceil pow min max clamp.");
                    if (!formulasample.error.empty())
                        ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.45f, 1.0f), "%s", formulasample.error.c_str());
                }

                if (ImGui::Button("Reset curve"))
                {
                    session.ui.hypsometricControlPoints =
                    {
                        ImVec2(0.0f, 0.0f),
                        ImVec2(0.22f, 0.28f),
                        ImVec2(0.78f, 0.72f),
                        ImVec2(1.0f, 1.0f)
                    };
                    session.ui.hypsometricControlWeights = { 1.0f, 1.0f, 1.0f, 1.0f };
                    session.ui.hypsometricSelectedPoint = 1;
                    session.ui.hypsometricCurveSource = HypsometricCurveSource::control_points;
                    session.ui.hypsometricFormula = "x";
                    result.controlsChanged = true;
                }

                ImGui::EndTabItem();
            }

            ImGui::EndTabBar();
        }
    }
    else if (stage.id == GenerationStageId::inland_seas)
    {
        seedcommittedwatercomponents(session);
        const vector<GenerationWaterComponent>& visiblecomponents = session.previewAvailable ? session.previewScratch.inlandSeaComponents : session.committedScratch.inlandSeaComponents;

        ImGui::TextUnformatted("Connected inland sea components");

        if (visiblecomponents.empty())
            ImGui::TextUnformatted("Preview this stage to detect components.");
        else
            ImGui::TextWrapped("Click inland seas on the map to cycle policies: red keeps the sea, green drains it, blue landfills it.");
    }
    else
    {
        ImGui::TextWrapped("Preview runs this stage on a scratch world cloned from the last committed state. Apply commits the result and advances.");
    }

    ImGui::Dummy(ImVec2(0.0f, 10.0f));
    if (ImGui::Button("Apply"))
        result.applyRequested = true;

    if (stage.skippable)
    {
        ImGui::SameLine();
        if (ImGui::Button("Skip"))
            result.skipRequested = true;
    }

    if (session.previewAvailable)
    {
        ImGui::SameLine();
        if (ImGui::Button("Reset preview"))
            result.resetPreviewRequested = true;
    }

    ImGui::BeginDisabled(canstepbackworkbenchsession(session) == false);
    if (ImGui::Button("Back"))
        result.backRequested = true;
    ImGui::EndDisabled();

    ImGui::SameLine();
    if (ImGui::Button("Abort"))
        result.abortRequested = true;

    ImGui::End();
    return result;
}
