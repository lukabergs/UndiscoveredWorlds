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

#include "generation_workbench.hpp"

#include "fastlem_mountains.hpp"
#include "map_appearance.hpp"
#include "physical_layers.hpp"
#include "ui_charts.hpp"

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

    ScopedWorkbenchDebugRun(long seed, int platecycles, int platecount)
    {
        if (isworldgendebugrunactive())
            return;

        options.logToProfilingWorkbook = false;
        options.visualizePlateTectonicsRealtime = false;
        options.plateTectonicsCycleCount = max(1, platecycles);
        options.plateTectonicsPlateCount = max(1, platecount);
        beginworldgendebugrun(seed, &options);
        started = true;
    }

    ~ScopedWorkbenchDebugRun()
    {
        if (started)
            endworldgendebugrun();
    }
};

vector<vector<int>> makeintgrid(int value = 0)
{
    return vector<vector<int>>(ARRAYWIDTH, vector<int>(ARRAYHEIGHT, value));
}

vector<vector<bool>> makeboolgrid(bool value = false)
{
    return vector<vector<bool>>(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, value));
}

vector<vector<std::uint8_t>> makefeaturegrid(std::uint8_t value = 0)
{
    return vector<vector<std::uint8_t>>(ARRAYWIDTH, vector<std::uint8_t>(ARRAYHEIGHT, value));
}

vector<vector<vector<int>>> makethreeintgrid()
{
    return vector<vector<vector<int>>>(ARRAYWIDTH, vector<vector<int>>(ARRAYHEIGHT, vector<int>(2, 0)));
}

void reseedterrainpass(planet& world, int salt)
{
    fast_srand(deterministicfastseed(deterministiccontextseed(world.seed(), salt)));
}

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

int wrapx(const planet& world, int x)
{
    return wrap(x, world.width());
}

PlateTectonicsSimulationOptions maketectonicsimulationoptions(const GenerationWorkbenchUiState& ui)
{
    PlateTectonicsSimulationOptions options;
    options.cycleCount = max(1, ui.tectonicCycleCount);
    options.cycleStepLimit = max(0, ui.tectonicCycleStepLimit);
    options.plateCount = max(1, ui.tectonicPlateCount);
    options.useSeaLevelMeters = ui.tectonicUseSeaLevelMeters;
    options.seaLevelMeters = clamp(ui.tectonicSeaLevelMeters, 0, 65535);
    options.aggregationOverlapAbsolute = ui.tectonicAggregationOverlapAbsolute;
    options.aggregationOverlapRelative = clamp(ui.tectonicAggregationOverlapRelative, 0.0f, 1.0f);
    options.foldingRatio = clamp(ui.tectonicFoldingRatio, 0.0f, 1.0f);
    options.erosionPeriod = max(1, ui.tectonicErosionPeriod);
    options.erosionStrength = max(0.0f, ui.tectonicErosionStrength);
    options.landmassRotation = max(0.0f, ui.tectonicLandmassRotation);
    options.rotationStrength = max(0.0f, ui.tectonicRotationStrength);
    options.subductionStrength = clamp(ui.tectonicSubductionStrength, 0.0f, 1.0f);
    options.divergentCarveStrength = max(0.0f, ui.tectonicDivergentCarveStrength);
    options.deltaTimeMyr = max(0.001f, ui.tectonicDeltaTimeMyr);
    return options;
}

int quantilethreshold(const vector<vector<int>>& fractal, int width, int height, float searatio)
{
    vector<int> values;
    values.reserve(static_cast<size_t>(width + 1) * static_cast<size_t>(height + 1));

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
            values.push_back(fractal[x][y]);
    }

    const size_t targetindex = min(values.size() - 1, static_cast<size_t>(roundf(clamp01(searatio) * static_cast<float>(values.size() - 1))));
    nth_element(values.begin(), values.begin() + targetindex, values.end());
    return values[targetindex];
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

float samplecurveaxis(const vector<ImVec2>& points, float input, bool samplexfromy)
{
    if (points.empty())
        return clamp01(input);

    const float clampedinput = clamp01(input);

    for (size_t index = 1; index < points.size(); index++)
    {
        const float lowaxis = samplexfromy ? points[index - 1].y : points[index - 1].x;
        const float highaxis = samplexfromy ? points[index].y : points[index].x;

        if (clampedinput > highaxis)
            continue;

        const float span = max(0.0001f, highaxis - lowaxis);
        const float factor = clamp01((clampedinput - lowaxis) / span);
        const float lowvalue = samplexfromy ? points[index - 1].x : points[index - 1].y;
        const float highvalue = samplexfromy ? points[index].x : points[index].y;
        return clamp01(lowvalue + (highvalue - lowvalue) * factor);
    }

    return samplexfromy ? points.back().x : points.back().y;
}

float solvebeziery(float inputx, const GenerationWorkbenchUiState& ui)
{
    return samplecurveaxis(samplehypsometriccurve(ui), inputx, false);
}

float sampleheightcdf(const vector<float>& values, float normalizedheight)
{
    if (values.empty())
        return clamp01(normalizedheight);

    const float position = clamp01(normalizedheight) * static_cast<float>(values.size() - 1);
    const size_t lowindex = static_cast<size_t>(floorf(position));
    const size_t highindex = min(lowindex + 1, values.size() - 1);
    const float factor = position - static_cast<float>(lowindex);
    return clamp01(values[lowindex] + (values[highindex] - values[lowindex]) * factor);
}

float solvetargetheightforcdf(float cumulativeprobability, const GenerationWorkbenchUiState& ui)
{
    return samplecurveaxis(samplehypsometriccurve(ui), cumulativeprobability, true);
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

void clearfeaturemask(GenerationScratch& scratch)
{
    scratch.stageFeatureMask = makefeaturegrid(0);
    scratch.stageFeatureCellCount = 0;
}

void markfeaturecell(GenerationScratch& scratch, int x, int y)
{
    if (x < 0 || x >= ARRAYWIDTH || y < 0 || y >= ARRAYHEIGHT)
        return;

    if (scratch.stageFeatureMask.empty())
        clearfeaturemask(scratch);

    if (scratch.stageFeatureMask[x][y] != 0)
        return;

    scratch.stageFeatureMask[x][y] = 1;
    scratch.stageFeatureCellCount++;
}

bool stageusesmanualpreview(GenerationStageId stageid)
{
    return stageid == GenerationStageId::inland_seas || stageid == GenerationStageId::basin_editor;
}

void rebuildheightdistribution(const planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    const int minimum = 1;
    const int maximum = world.maxelevation() - 1;

    scratch.heightHistogram.assign(histogrambincount, 0);
    scratch.heightCdf.assign(histogrambincount, 0.0f);

    int total = 0;

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            const float normalized = normalizeterrainvalue(world.nom(x, y), minimum, maximum);
            int bin = static_cast<int>(floorf(normalized * static_cast<float>(histogrambincount)));
            bin = clamp(bin, 0, histogrambincount - 1);
            scratch.heightHistogram[bin]++;
            total++;
        }
    }

    if (total <= 0)
        return;

    int running = 0;

    for (int index = 0; index < histogrambincount; index++)
    {
        running += scratch.heightHistogram[index];
        scratch.heightCdf[index] = static_cast<float>(running) / static_cast<float>(total);
    }
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

void initializegenerationscratch(GenerationScratch& scratch)
{
    scratch.mountaindrainage = makeintgrid();
    scratch.shelves = makeboolgrid(false);
    scratch.squareroot.assign((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS + 1) * 24, 0);

    for (size_t index = 1; index < scratch.squareroot.size(); index++)
        scratch.squareroot[index] = static_cast<int>(sqrt(static_cast<float>(index)));

    scratch.saltLakeMap = makethreeintgrid();
    scratch.noLake = makeintgrid();
    scratch.basinSeeds = makeintgrid();
    scratch.inlandSeaComponentIds = makeintgrid(0);
    scratch.basinComponentIds = makeintgrid(0);
    scratch.heightHistogram.assign(histogrambincount, 0);
    scratch.heightCdf.assign(histogrambincount, 0.0f);
    scratch.inlandSeaComponents.clear();
    scratch.basinComponents.clear();
    clearfeaturemask(scratch);
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

void initializeproceduralworkbenchsession(GenerationSessionState& session, const planet& world, int plateCycles, int plateCount)
{
    resetworkbenchsession(session);
    session.active = true;
    session.procedural = true;
    session.previewAvailable = false;
    session.currentStageIndex = getgenerationstageindex(GenerationStageId::plate_tectonics);
    session.previewWorld = world;
    initializegenerationscratch(session.committedScratch);
    initializegenerationscratch(session.previewScratch);
    session.ui.tectonicSeed = static_cast<int>(world.seed());
    session.ui.tectonicCycleCount = plateCycles;
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
    case GenerationStageId::plate_tectonics:
    case GenerationStageId::tectonic_trenches:
    case GenerationStageId::tectonic_volcanoes:
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

    case GenerationStageId::rivers_and_basins:
    case GenerationStageId::basin_editor:
    case GenerationStageId::lakes:
    case GenerationStageId::post_river_fastlem:
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
void rebuildshelves(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
                scratch.shelves[x][y] = false;
        }
    });

    makecontinentalshelves(world, scratch.shelves, 4);
}

void applyhypsometricremap(planet& world, GenerationScratch& scratch, const GenerationWorkbenchUiState& ui)
{
    const int width = world.width();
    const int height = world.height();
    const int minimum = 1;
    const int maximum = world.maxelevation() - 1;
    const vector<float> sourcecdf = scratch.heightCdf;
    const vector<ImVec2> targetcurve = samplehypsometriccurve(ui);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const float normalized = normalizeterrainvalue(world.nom(x, y), minimum, maximum);
                const float quantile = sampleheightcdf(sourcecdf, normalized);
                const float remapped = samplecurveaxis(targetcurve, quantile, false);
                const int newvalue = minimum + static_cast<int>(roundf(remapped * static_cast<float>(maximum - minimum)));
                world.setnom(x, y, clamp(newvalue, minimum, maximum));
            }
        }
    });

    rebuildheightdistribution(world, scratch);
}

void applyterraforming(planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui)
{
    world.setsealevel(clamp(ui.seaLevel, 1, world.maxelevation() - 2));

    if (ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
        applyterrainrenderpreset(ui, world, ui.terrainRenderPreset);

    if (ui.useHypsometricRemap)
        applyhypsometricremap(world, scratch, ui);

    rebuildshelves(world, scratch);
    rebuildheightdistribution(world, scratch);
}

vector<vector<int>> buildterrainfractal(planet& world, int salt, int maximum)
{
    vector<vector<int>> fractal = makeintgrid();
    reseedterrainpass(world, salt);
    createfractal(fractal, world.width(), world.height(), 8, 0.2f, static_cast<float>(random(3, 6)), 1, maximum, 0, 0);
    warp(fractal, world.width(), world.height(), world.maxelevation(), random(20, 80), 1);
    return fractal;
}

void metadataoceanrefinement(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    const int maxelev = world.maxelevation();
    const int sealevel = world.sealevel();
    vector<vector<int>> seafractal = makeintgrid();

    reseedterrainpass(world, 0x4110);
    createfractal(seafractal, width, height, 8, 0.2f, static_cast<float>(random(3, 6)), 1, maxelev, 0, 0);
    warp(seafractal, width, height, maxelev, random(20, 80), 1);

    const float coastalvarreduce = static_cast<float>(maxelev) / 500.0f;
    const float oceanvarreduce = static_cast<float>(maxelev) / 1000.0f;

    for (int x = 0; x <= width; x++)
    {
        for (int y = 0; y <= height; y++)
        {
            if (world.sea(x, y) == 0)
                continue;

            if (scratch.shelves[x][y])
            {
                const float var = static_cast<float>(seafractal[x][y] - maxelev / 2) / coastalvarreduce;
                int newvalue = sealevel - 200 + static_cast<int>(var);
                newvalue = min(newvalue, sealevel - 10);
                world.setnom(x, y, max(1, newvalue));
            }
            else
            {
                int xx = x + width / 2;
                if (xx > width)
                    xx -= width;

                const float var = static_cast<float>(seafractal[xx][y] - maxelev / 2) / oceanvarreduce;
                int newvalue = sealevel - 5000 + static_cast<int>(var);
                newvalue = min(newvalue, sealevel - 3000);
                world.setnom(x, y, max(1, newvalue));
            }
        }
    }
}

int wrapdistx(const planet& world, int left, int right)
{
    const int width = world.width() + 1;
    int delta = abs(left - right);
    delta = min(delta, width - delta);
    return delta;
}

void metadataactivetrenches(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    clearfeaturemask(scratch);
    reseedterrainpass(world, 0x4114);

    int maxsegmentid = 0;
    for (const TectonicBoundarySegment& segment : world.tectonicboundarysegments())
        maxsegmentid = max(maxsegmentid, segment.id);

    vector<int> segmentlookup(static_cast<size_t>(maxsegmentid + 1), -1);
    vector<int> segmentids;
    vector<vector<pair<int, int>>> segmentcells;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 0 || world.tectonicboundarytype(x, y) != BoundaryType::convergent)
                continue;

            const CrustClass crustclass = world.tectoniccrustclass(x, y);
            if (crustclass != CrustClass::oceanic && crustclass != CrustClass::transitional)
                continue;

            const int segmentid = world.tectonicboundarysegmentid(x, y);
            if (segmentid <= 0)
                continue;

            if (segmentlookup[segmentid] < 0)
            {
                segmentlookup[segmentid] = static_cast<int>(segmentcells.size());
                segmentids.push_back(segmentid);
                segmentcells.emplace_back();
            }

            segmentcells[segmentlookup[segmentid]].emplace_back(x, y);
        }
    }

    vector<int> order(segmentids.size(), 0);
    for (size_t index = 0; index < order.size(); index++)
        order[index] = static_cast<int>(index);

    for (int index = static_cast<int>(order.size()) - 1; index > 0; index--)
        swap(order[index], order[random(0, index)]);

    const int trenchcount = min(static_cast<int>(order.size()), clamp(static_cast<int>(order.size()) / 3 + 1, 2, 6));

    for (int trenchindex = 0; trenchindex < trenchcount; trenchindex++)
    {
        const vector<pair<int, int>>& cells = segmentcells[order[trenchindex]];
        if (cells.empty())
            continue;

        const int branchcount = clamp(static_cast<int>(cells.size()) / 48, 1, 3);

        for (int branch = 0; branch < branchcount; branch++)
        {
            const pair<int, int>& anchor = cells[random(0, static_cast<int>(cells.size()) - 1)];
            const float maxsegmentdistance = static_cast<float>(random(18, 38));
            const int radius = random(2, 4);

            for (const pair<int, int>& cell : cells)
            {
                const int dx = wrapdistx(world, cell.first, anchor.first);
                const int dy = abs(cell.second - anchor.second);
                const float segmentdistance = sqrtf(static_cast<float>(dx * dx + dy * dy));

                if (segmentdistance > maxsegmentdistance)
                    continue;

                const float branchstrength = 1.0f - clamp01(segmentdistance / maxsegmentdistance);

                for (int ox = -radius; ox <= radius; ox++)
                {
                    for (int oy = -radius; oy <= radius; oy++)
                    {
                        const int nx = wrapx(world, cell.first + ox);
                        const int ny = cell.second + oy;

                        if (ny < 0 || ny > height || world.sea(nx, ny) == 0)
                            continue;

                        const float kerneldistance = sqrtf(static_cast<float>(ox * ox + oy * oy));
                        if (kerneldistance > static_cast<float>(radius) + 0.5f)
                            continue;

                        const float kernelfalloff = 1.0f - clamp01(kerneldistance / (static_cast<float>(radius) + 0.5f));
                        const float proximity = 1.0f - clamp01(static_cast<float>(world.tectonicboundarydistance(nx, ny)) / 16.0f);
                        const float subsidence = clamp01(world.tectonicsubsidencetendency(nx, ny));
                        const float convergence = clamp01(static_cast<float>(world.tectonicconvergence(nx, ny)) / 255.0f);
                        const float randomness = 0.82f + static_cast<float>(random(0, 35)) / 100.0f;
                        const int trenchdepth = static_cast<int>(roundf((1200.0f + 2200.0f * branchstrength) * kernelfalloff * randomness * (0.45f + subsidence * 0.30f + convergence * 0.25f)));

                        if (trenchdepth <= 0)
                            continue;

                        const int target = max(1, min(world.nom(nx, ny), sealevel - 450 - trenchdepth));
                        if (target >= world.nom(nx, ny))
                            continue;

                        world.setnom(nx, ny, target);
                        world.setgeologicregime(nx, ny, GeologicRegime::trench_adjacent);
                        markfeaturecell(scratch, nx, ny);
                    }
                }
            }
        }
    }
}

void applyvolcanicupliftkernel(planet& world, GenerationScratch& scratch, int centrex, int centrey, int radius, int maximumuplift, GeologicRegime regime)
{
    for (int dx = -radius; dx <= radius; dx++)
    {
        for (int dy = -radius; dy <= radius; dy++)
        {
            const int nx = wrapx(world, centrex + dx);
            const int ny = centrey + dy;

            if (ny < 0 || ny > world.height())
                continue;

            const float distance = sqrtf(static_cast<float>(dx * dx + dy * dy));
            if (distance > static_cast<float>(radius) + 0.5f)
                continue;

            const float falloff = 1.0f - clamp01(distance / (static_cast<float>(radius) + 0.5f));
            const int uplift = static_cast<int>(roundf(static_cast<float>(maximumuplift) * falloff * falloff));
            if (uplift <= 0)
                continue;

            const int target = min(world.maxelevation() - 1, world.nom(nx, ny) + uplift);
            if (target <= world.nom(nx, ny))
                continue;

            world.setnom(nx, ny, target);
            if (regime != GeologicRegime::stable)
                world.setgeologicregime(nx, ny, regime);
            markfeaturecell(scratch, nx, ny);
        }
    }
}

void metadatavolcanoes(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    clearfeaturemask(scratch);
    rebuildshelves(world, scratch);
    reseedterrainpass(world, 0x4115);

    vector<vector<int>> volcanodirection = buildterrainfractal(world, 0x4116, world.maxelevation());
    vector<pair<int, int>> selectedseeds;

    auto farenough = [&](int x, int y, int minimumdistance)
    {
        for (const pair<int, int>& seed : selectedseeds)
        {
            const int dx = wrapdistx(world, x, seed.first);
            const int dy = abs(y - seed.second);
            if (dx * dx + dy * dy < minimumdistance * minimumdistance)
                return false;
        }

        return true;
    };

    struct VolcanoCandidate
    {
        int x = 0;
        int y = 0;
        float score = 0.0f;
    };

    vector<VolcanoCandidate> conecandidates;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.tectonicboundarytype(x, y) != BoundaryType::convergent)
                continue;

            const float proximity = 1.0f - clamp01(static_cast<float>(world.tectonicboundarydistance(x, y)) / 12.0f);
            const float uplift = clamp01(world.tectonicuplifttendency(x, y));
            const float convergence = clamp01(static_cast<float>(world.tectonicconvergence(x, y)) / 255.0f);
            const float history = clamp01(world.tectonicboundaryhistory(x, y));
            const float score = proximity * (0.35f + uplift * 0.35f + convergence * 0.20f + history * 0.10f);

            if (score < 0.34f)
                continue;

            conecandidates.push_back({ x, y, score });
        }
    }

    for (int index = static_cast<int>(conecandidates.size()) - 1; index > 0; index--)
        swap(conecandidates[index], conecandidates[random(0, index)]);

    const int area = (width + 1) * (height + 1);
    const int desiredcones = clamp(area / 200000, 2, 10);

    for (const VolcanoCandidate& candidate : conecandidates)
    {
        if (static_cast<int>(selectedseeds.size()) >= desiredcones)
            break;

        if (!farenough(candidate.x, candidate.y, 24))
            continue;

        const int peakheight = static_cast<int>(roundf(1800.0f + candidate.score * 3800.0f));
        applyvolcanicupliftkernel(world, scratch, candidate.x, candidate.y, world.sea(candidate.x, candidate.y) ? 2 : 1, static_cast<int>(roundf(240.0f + candidate.score * 520.0f)), GeologicRegime::convergent_arc);
        createisolatedvolcano(world, candidate.x, candidate.y, scratch.shelves, volcanodirection, peakheight, true);
        selectedseeds.emplace_back(candidate.x, candidate.y);
    }

    auto placenexthotspot = [&](bool oceanpreferred, int targetcount)
    {
        for (int placed = 0, attempts = 0; placed < targetcount && attempts < 512; attempts++)
        {
            const int x = random(0, width);
            const int y = random(0, height);

            if (!farenough(x, y, 34))
                continue;

            if (world.tectonicboundarydistance(x, y) < 18)
                continue;

            if (oceanpreferred != (world.sea(x, y) == 1))
                continue;

            if (world.sea(x, y) == 1 && scratch.shelves[x][y])
                continue;

            const int radius = oceanpreferred ? random(4, 6) : random(3, 5);
            const int uplift = oceanpreferred ? random(260, 700) : random(180, 460);
            const int peakheight = oceanpreferred ? random(900, 2400) : random(700, 1600);
            applyvolcanicupliftkernel(world, scratch, x, y, radius, uplift, GeologicRegime::stable);
            createisolatedvolcano(world, x, y, scratch.shelves, volcanodirection, peakheight, false);
            selectedseeds.emplace_back(x, y);
            placed++;
        }
    };

    placenexthotspot(true, clamp(area / 280000 + 1, 1, 6));
    placenexthotspot(false, clamp(area / 520000, 0, 3));

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.volcano(x, y) == 0)
                continue;

            markfeaturecell(scratch, x, y);
            if (world.strato(x, y))
                world.setgeologicregime(x, y, GeologicRegime::convergent_arc);
        }
    }
}

void detectinlandseas(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    scratch.inlandSeaComponentIds = makeintgrid(0);
    scratch.inlandSeaComponents.clear();

    vector<vector<bool>> visited = makeboolgrid(false);
    vector<vector<pair<int, int>>> components(1);
    int largestcomponentid = 0;
    int largestcomponentarea = 0;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (visited[x][y] || world.sea(x, y) == 0)
                continue;

            queue<pair<int, int>> frontier;
            vector<pair<int, int>> cells;
            frontier.emplace(x, y);
            visited[x][y] = true;

            while (frontier.empty() == false)
            {
                const pair<int, int> current = frontier.front();
                frontier.pop();
                cells.push_back(current);

                const int directions[8][2] =
                {
                    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
                    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
                };

                for (const auto& direction : directions)
                {
                    const int nx = wrapx(world, current.first + direction[0]);
                    const int ny = current.second + direction[1];

                    if (ny < 0 || ny > height || visited[nx][ny] || world.sea(nx, ny) == 0)
                        continue;

                    visited[nx][ny] = true;
                    frontier.emplace(nx, ny);
                }
            }

            components.push_back(std::move(cells));

            if (static_cast<int>(components.back().size()) > largestcomponentarea)
            {
                largestcomponentarea = static_cast<int>(components.back().size());
                largestcomponentid = static_cast<int>(components.size() - 1);
            }
        }
    }

    for (size_t componentindex = 1; componentindex < components.size(); componentindex++)
    {
        if (static_cast<int>(componentindex) == largestcomponentid)
            continue;

        GenerationWaterComponent component;
        component.id = static_cast<int>(scratch.inlandSeaComponents.size() + 1);
        component.seedX = components[componentindex].front().first;
        component.seedY = components[componentindex].front().second;
        component.cellCount = static_cast<int>(components[componentindex].size());
        scratch.inlandSeaComponents.push_back(component);

        for (const pair<int, int>& cell : components[componentindex])
            scratch.inlandSeaComponentIds[cell.first][cell.second] = component.id;
    }
}

vector<pair<int, int>> findcomponentdrainpath(planet& world, const vector<vector<int>>& componentids, int componentid)
{
    const int width = world.width();
    const int height = world.height();
    vector<vector<bool>> visited = makeboolgrid(false);
    vector<vector<int>> parentx = makeintgrid(-1);
    vector<vector<int>> parenty = makeintgrid(-1);
    queue<pair<int, int>> frontier;
    pair<int, int> endpoint = { -1, -1 };

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (componentids[x][y] == componentid)
            {
                frontier.emplace(x, y);
                visited[x][y] = true;
                parentx[x][y] = x;
                parenty[x][y] = y;
            }
        }
    }

    const int directions[8][2] =
    {
        { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
        { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
    };

    while (frontier.empty() == false && endpoint.first < 0)
    {
        const pair<int, int> current = frontier.front();
        frontier.pop();

        for (const auto& direction : directions)
        {
            const int nx = wrapx(world, current.first + direction[0]);
            const int ny = current.second + direction[1];

            if (ny < 0 || ny > height || visited[nx][ny])
                continue;

            visited[nx][ny] = true;
            parentx[nx][ny] = current.first;
            parenty[nx][ny] = current.second;

            if (world.sea(nx, ny) == 1 && componentids[nx][ny] == 0)
            {
                endpoint = { nx, ny };
                break;
            }

            frontier.emplace(nx, ny);
        }
    }

    vector<pair<int, int>> path;

    if (endpoint.first < 0)
        return path;

    int x = endpoint.first;
    int y = endpoint.second;

    while (parentx[x][y] != x || parenty[x][y] != y)
    {
        path.emplace_back(x, y);
        const int px = parentx[x][y];
        const int py = parenty[x][y];
        x = px;
        y = py;
    }

    path.emplace_back(x, y);
    reverse(path.begin(), path.end());
    return path;
}

void applyinlandseapolicies(planet& world, GenerationScratch& scratch)
{
    if (scratch.inlandSeaComponents.empty())
        detectinlandseas(world, scratch);

    for (const GenerationWaterComponent& component : scratch.inlandSeaComponents)
    {
        if (component.policy == GenerationComponentPolicy::keep)
            continue;

        if (component.policy == GenerationComponentPolicy::fill)
        {
            for (int x = 0; x <= world.width(); x++)
            {
                for (int y = 0; y <= world.height(); y++)
                {
                    if (scratch.inlandSeaComponentIds[x][y] == component.id)
                        world.setnom(x, y, max(world.sealevel() + 1, world.nom(x, y)));
                }
            }
        }
        else
        {
            const vector<pair<int, int>> path = findcomponentdrainpath(world, scratch.inlandSeaComponentIds, component.id);

            for (const pair<int, int>& cell : path)
                world.setnom(cell.first, cell.second, max(1, world.sealevel() - 10));
        }
    }

    rebuildshelves(world, scratch);
    detectinlandseas(world, scratch);
}

void detectbasincomponents(planet& world, GenerationScratch& scratch)
{
    const int width = world.width();
    const int height = world.height();
    scratch.basinComponentIds = makeintgrid(0);
    scratch.basinComponents.clear();
    vector<vector<bool>> visited = makeboolgrid(false);
    int nextid = 1;

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (visited[x][y] || world.basinclass(x, y) != BasinClass::endorheic)
                continue;

            queue<pair<int, int>> frontier;
            frontier.emplace(x, y);
            visited[x][y] = true;
            GenerationBasinComponent component;
            component.id = nextid;
            component.seedX = x;
            component.seedY = y;
            component.basinClass = BasinClass::endorheic;

            while (frontier.empty() == false)
            {
                const pair<int, int> current = frontier.front();
                frontier.pop();
                component.cellCount++;
                scratch.basinComponentIds[current.first][current.second] = nextid;

                const int directions[8][2] =
                {
                    { 0, -1 }, { 1, -1 }, { 1, 0 }, { 1, 1 },
                    { 0, 1 }, { -1, 1 }, { -1, 0 }, { -1, -1 }
                };

                for (const auto& direction : directions)
                {
                    const int nx = wrapx(world, current.first + direction[0]);
                    const int ny = current.second + direction[1];

                    if (ny < 0 || ny > height || visited[nx][ny] || world.basinclass(nx, ny) != BasinClass::endorheic)
                        continue;

                    visited[nx][ny] = true;
                    frontier.emplace(nx, ny);
                }
            }

            scratch.basinComponents.push_back(component);
            nextid++;
        }
    }
}

void applybasinpolicies(planet& world, GenerationScratch& scratch)
{
    if (scratch.basinComponents.empty())
        detectbasincomponents(world, scratch);

    for (const GenerationBasinComponent& component : scratch.basinComponents)
    {
        if (component.policy == GenerationComponentPolicy::keep)
            continue;

        if (component.policy == GenerationComponentPolicy::fill)
        {
            int targetelevation = 0;

            for (int x = 0; x <= world.width(); x++)
            {
                for (int y = 0; y <= world.height(); y++)
                {
                    if (scratch.basinComponentIds[x][y] == component.id)
                        targetelevation = max(targetelevation, world.nom(x, y));
                }
            }

            for (int x = 0; x <= world.width(); x++)
            {
                for (int y = 0; y <= world.height(); y++)
                {
                    if (scratch.basinComponentIds[x][y] == component.id)
                        world.setnom(x, y, targetelevation);
                }
            }
        }
        else
        {
            const vector<pair<int, int>> path = findcomponentdrainpath(world, scratch.basinComponentIds, component.id);

            for (const pair<int, int>& cell : path)
                world.setnom(cell.first, cell.second, max(1, world.sealevel() - 10));
        }
    }

    detectbasincomponents(world, scratch);
}

void applyterraintexturing(planet& world, GenerationScratch& scratch)
{
    clearfeaturemask(scratch);
    vector<vector<int>> fractal = buildterrainfractal(world, 0x4111, world.maxelevation());
    flip(fractal, world.width(), world.height(), 1, 1);
    shift(fractal, world.width(), world.height(), random(1, world.width()));

    const int sealevel = world.sealevel();
    const int midpoint = world.maxelevation() / 2;

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
        {
            const int oldvalue = world.nom(x, y);
            const int noise = fractal[x][y] - midpoint;
            const float strength = world.sea(x, y) ? 0.14f : 0.10f;
            int newvalue = oldvalue + static_cast<int>(roundf(static_cast<float>(noise) * strength));

            if (world.sea(x, y))
                newvalue = min(newvalue, sealevel - 10);
            else
                newvalue = max(newvalue, sealevel + 5);

            newvalue = clamp(newvalue, 1, world.maxelevation() - 1);

            if (newvalue == oldvalue)
                continue;

            world.setnom(x, y, newvalue);
            markfeaturecell(scratch, x, y);
        }
    }
}

bool applyfastlemmountains(planet& world)
{
    vector<vector<int>> fractal = buildterrainfractal(world, 0x4112, 12750);
    return generatefastlemmountains(world, fractal);
}

void rebuildterrainroughness(planet& world)
{
    vector<vector<int>> roughness = buildterrainfractal(world, 0x4113, world.maxelevation());

    for (int x = 0; x <= world.width(); x++)
    {
        for (int y = 0; y <= world.height(); y++)
            world.setroughness(x, y, static_cast<float>(roughness[x][y]));
    }
}

void applytemperaturestage(planet& world)
{
    const int maxelev = world.maxelevation();
    vector<vector<int>> fractal = makeintgrid();
    reseedterrainpass(world, 0x5002);
    createfractal(fractal, world.width(), world.height(), 8, 0.2f, static_cast<float>(random(1, 4)), 1, maxelev, 0, 0);
    warp(fractal, world.width(), world.height(), maxelev, random(40, 80), 0);
    createtemperaturemap(world, fractal);
    world.syncseasonalclimatefromlegacy();
}

void applyrainfallstage(planet& world, const GenerationExecutionContext& context)
{
    const int maxelev = world.maxelevation();
    vector<vector<int>> fractal = makeintgrid();
    reseedterrainpass(world, 0x5004);
    createfractal(fractal, world.width(), world.height(), 8, 0.2f, 0.6f, 1, maxelev, 0, 0);
    warp(fractal, world.width(), world.height(), maxelev, random(40, 80), 0);
    warp(fractal, world.width(), world.height(), maxelev, random(40, 80), 0);
    getlandandseatotals(world);
    createrainmap(world, fractal, world.landtotal(), world.seatotal(), context.smalllake, context.landshape);
}

void applyriverandbasinstage(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    if (context.dorivers == false)
        return;

    const int width = world.width();
    const int height = world.height();
    int saltlakesplaced = 0;
    scratch.saltLakeMap = makethreeintgrid();
    scratch.noLake = makeintgrid();
    scratch.basinSeeds = makeintgrid();

    createrivermap(world, scratch.mountaindrainage);

    if (context.dolakes || world.seatotal() == 0)
    {
        createsaltlakes(world, saltlakesplaced, scratch.saltLakeMap, scratch.noLake, scratch.basinSeeds, context.smalllake);
        addlandnoise(world);
        depressionfill(world);

        for (int x = 0; x <= width; x++)
        {
            for (int y = 0; y <= height; y++)
            {
                world.setriverdir(x, y, 0);
                world.setriverjan(x, y, 0);
                world.setriverjul(x, y, 0);
            }
        }

        createrivermap(world, scratch.mountaindrainage);
    }

    world.setmaxriverflow();
    detectbasincomponents(world, scratch);
}

void applylakestage(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    if (context.dolakes == false)
    {
        if (world.seatotal() == 0)
            convertsaltlakes(world, scratch.saltLakeMap);

        return;
    }

    convertsaltlakes(world, scratch.saltLakeMap);
    createlakemap(world, scratch.noLake, context.smalllake, context.largelake);
    createriftlakemap(world, scratch.noLake);
    world.setmaxriverflow();
}

void applyaridfeatures(planet& world, const GenerationExecutionContext& context)
{
    createergs(world, context.smalllake, context.largelake, context.landshape);
    createsaltpans(world, context.smalllake, context.largelake);
}

void applydeltawetlandstage(planet& world, const GenerationExecutionContext& context)
{
    if (context.dodeltas)
    {
        createriverdeltas(world);
        checkrivers(world);
    }

    createwetlands(world, context.smalllake);
    removeexcesswetlands(world);
    refineroughnessmap(world);
    removesealakes(world);
    connectlakes(world);
    checkpoleclimates(world);
    world.syncseasonalclimatefromlegacy();
    applyseasonaltemperaturelapse(world);
    createclimatemap(world);
    createbiomemap(world);
    exportclimatevalidationreport(world);
}

bool executegenerationstage(GenerationStageId stageId, planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui, const GenerationExecutionContext& context, string* errormessage)
{
    const GenerationStageDefinition& stagedefinition = getgenerationstage(stageId);

    if (beginworldgenstep(stagedefinition.label) == false)
        return true;

    clearfeaturemask(scratch);

    switch (stageId)
    {
    case GenerationStageId::plate_tectonics:
        world.clear();
        applyplatetectonicssimulation(world, scratch.shelves, maketectonicsimulationoptions(ui));
        rebuildshelves(world, scratch);
        ui.seaLevel = world.sealevel();
        rebuildheightdistribution(world, scratch);
        return true;

    case GenerationStageId::terraforming:
        applyterraforming(world, scratch, ui);
        return true;

    case GenerationStageId::inland_seas:
        applyinlandseapolicies(world, scratch);
        return true;

    case GenerationStageId::archipelagos:
    {
        vector<vector<bool>> removedland = makeboolgrid(false);
        makearchipelagos(world, removedland, context.landshape);
        return true;
    }

    case GenerationStageId::ocean_cleanup:
        widenchannels(world);
        loweroceans(world);
        rebuildshelves(world, scratch);
        return true;

    case GenerationStageId::coastline_refinement:
        removestraights(world);
        return true;

    case GenerationStageId::ocean_depth_refinement:
        metadataoceanrefinement(world, scratch);
        return true;

    case GenerationStageId::tectonic_trenches:
        metadataactivetrenches(world, scratch);
        return true;

    case GenerationStageId::tectonic_volcanoes:
        metadatavolcanoes(world, scratch);
        return true;

    case GenerationStageId::terrain_texturing:
        applyterraintexturing(world, scratch);
        return true;

    case GenerationStageId::fastlem_mountains:
        if (applyfastlemmountains(world) == false && errormessage != nullptr)
            *errormessage = "FastLEM mountains did not produce enough ridge candidates for this seed.";
        return true;

    case GenerationStageId::mountain_bases:
    {
        vector<vector<bool>> defaultmountains = makeboolgrid(false);
        vector<vector<bool>>& mountainsok = context.okmountains != nullptr ? *context.okmountains : defaultmountains;
        raisemountainbases(world, scratch.mountaindrainage, mountainsok);
        return true;
    }

    case GenerationStageId::terrain_smoothing:
        smoothland(world, 2);
        return true;

    case GenerationStageId::canyon_uplift:
        createextraelev(world);
        return true;

    case GenerationStageId::depression_fill:
        depressionfill(world);
        addlandnoise(world);
        depressionfill(world);
        return true;

    case GenerationStageId::coastline_adjustment:
        normalisecoasts(world, 13, 11, 4);
        normalisecoasts(world, 13, 11, 4);
        return true;

    case GenerationStageId::island_check:
        checkislands(world);
        extendnoshade(world);
        return true;

    case GenerationStageId::terrain_roughness:
        getlandandseatotals(world);
        rebuildterrainroughness(world);
        return true;

    case GenerationStageId::global_temperature:
        applytemperaturestage(world);
        return true;

    case GenerationStageId::ocean_currents:
        if (world.seatotal() > 0)
            createoceancurrentmap(world);
        return true;

    case GenerationStageId::sea_surface_temperatures:
        if (world.seatotal() > 0)
            createsurfacetemperaturemap(world);
        return true;

    case GenerationStageId::pressure:
        if (world.seatotal() > 0)
            createpressuremap(world);
        return true;

    case GenerationStageId::winds:
        if (world.seatotal() > 0)
        {
            createvectorwindmap(world);
            updatehorsebeltsfrompressure(world);
        }
        return true;

    case GenerationStageId::sea_ice_and_tides:
        if (world.seatotal() > 0)
        {
            vector<vector<int>> fractal = buildterrainfractal(world, 0x5003, world.maxelevation());
            createseaicemap(world, fractal);
            createtidalmap(world);
        }
        return true;

    case GenerationStageId::rainfall:
        applyrainfallstage(world, context);
        return true;

    case GenerationStageId::fjords:
        if (world.seatotal() > 0)
            addfjordmountains(world);
        return true;

    case GenerationStageId::rivers_and_basins:
        applyriverandbasinstage(world, scratch, context);
        return true;

    case GenerationStageId::basin_editor:
        applybasinpolicies(world, scratch);
        if (context.dorivers)
            createrivermap(world, scratch.mountaindrainage);
        return true;

    case GenerationStageId::lakes:
        applylakestage(world, scratch, context);
        return true;

    case GenerationStageId::post_river_fastlem:
        broadenfastlemterrainfromrivers(world);
        return true;

    case GenerationStageId::mountain_temperature_lapse:
        checkpoleclimates(world);
        world.syncseasonalclimatefromlegacy();
        applyseasonaltemperaturelapse(world);
        return true;

    case GenerationStageId::climates_and_biomes:
        createclimatemap(world);
        createbiomemap(world);
        return true;

    case GenerationStageId::arid_features:
        applyaridfeatures(world, context);
        return true;

    case GenerationStageId::deltas_wetlands_and_roughness:
        applydeltawetlandstage(world, context);
        return true;

    case GenerationStageId::finalize_layers:
        finalizegeneratedworld(world, scratch, context);
        return true;
    }

    return false;
}
}

namespace
{
void advanceworkbenchstage(GenerationSessionState& session)
{
    if (session.currentStageIndex < getgenerationstages().size())
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

void seedcommittedbasincomponents(GenerationSessionState& session)
{
    if (session.committedScratch.basinComponents.empty() && session.previewAvailable && session.previewScratch.basinComponents.empty() == false)
    {
        session.committedScratch.basinComponentIds = session.previewScratch.basinComponentIds;
        session.committedScratch.basinComponents = session.previewScratch.basinComponents;
    }
}

const char* policylabel(GenerationComponentPolicy policy)
{
    switch (policy)
    {
    case GenerationComponentPolicy::keep:
        return "Keep";
    case GenerationComponentPolicy::fill:
        return "Fill";
    case GenerationComponentPolicy::drain:
        return "Drain";
    }

    return "Keep";
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
    const long debugseed = stageid == GenerationStageId::plate_tectonics ? static_cast<long>(session.ui.tectonicSeed) : committedworld.seed();
    ScopedWorkbenchDebugRun debugrun(debugseed, session.ui.tectonicCycleCount, session.ui.tectonicPlateCount);
    session.previewWorld = committedworld;
    if (stageid == GenerationStageId::plate_tectonics)
        session.previewWorld.setseed(debugseed);
    session.previewScratch = session.committedScratch;
    session.previewAvailable = true;

    if (stageusesmanualpreview(stageid))
    {
        if (stageid == GenerationStageId::inland_seas && session.previewScratch.inlandSeaComponents.empty())
            detectinlandseas(session.previewWorld, session.previewScratch);
        else if (stageid == GenerationStageId::basin_editor && session.previewScratch.basinComponents.empty())
            detectbasincomponents(session.previewWorld, session.previewScratch);

        if (isterraingenerationstage(stageid) && session.ui.terrainRenderPreset != WorkbenchTerrainRenderPreset::custom)
            applyterrainrenderpreset(session.ui, session.previewWorld, session.ui.terrainRenderPreset);

        return true;
    }

    session.previewAvailable = executegenerationstage(stageid, session.previewWorld, session.previewScratch, session.ui, context, errormessage);
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
        if (executegenerationstage(getcurrentgenerationstage(session).id, session.previewWorld, session.previewScratch, session.ui, context, errormessage) == false)
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

bool runremaininggenerationstages(planet& world, GenerationScratch& scratch, GenerationWorkbenchUiState& ui, const GenerationExecutionContext& context, size_t startStageIndex, string* errormessage)
{
    ScopedWorkbenchDebugRun debugrun(world.seed(), ui.tectonicCycleCount, ui.tectonicPlateCount);
    const vector<GenerationStageDefinition>& stages = getgenerationstages();

    for (size_t index = startStageIndex; index < stages.size(); index++)
    {
        if (executegenerationstage(stages[index].id, world, scratch, ui, context, errormessage) == false)
            return false;
    }

    return true;
}

void finalizegeneratedworld(planet& world, GenerationScratch& scratch, const GenerationExecutionContext& context)
{
    generatephysicalworldlayers(world, scratch.shelves);
    generatesocialworld(world, context.socialoptions);

    if (context.appendclimateworkbook && appendclimatebenchmarkworkbook(world) == false)
        updatereport("Climate workbook benchmark update failed");
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

    if (stage.id == GenerationStageId::plate_tectonics)
    {
        result.discardPreviewRequested |= ImGui::InputInt("Seed", &session.ui.tectonicSeed);
        if (session.ui.tectonicSeed < 0)
            session.ui.tectonicSeed = 0;
        settooltipifhovered("Seed used for plate tectonics. The same seed yields the same result.");

        ImGui::SameLine();
        if (ImGui::Button("Recompute"))
            result.recomputeRequested = true;
        settooltipifhovered("Rerun the plate tectonics stage with the current seed and parameters.");

        result.discardPreviewRequested |= ImGui::InputInt("Plate cycles", &session.ui.tectonicCycleCount);
        if (session.ui.tectonicCycleCount < 1)
            session.ui.tectonicCycleCount = 1;
        settooltipifhovered("Number of plate-tectonics cycles to run. Upstream default is 2.");

        result.discardPreviewRequested |= ImGui::InputInt("Cycle step limit", &session.ui.tectonicCycleStepLimit);
        if (session.ui.tectonicCycleStepLimit < 0)
            session.ui.tectonicCycleStepLimit = 0;
        settooltipifhovered("Maximum simulation updates per cycle before restart. Set 0 to disable the limit. Upstream default is 600.");

        result.discardPreviewRequested |= ImGui::InputInt("Plate count", &session.ui.tectonicPlateCount);
        if (session.ui.tectonicPlateCount < 1)
            session.ui.tectonicPlateCount = 1;
        settooltipifhovered("Number of tectonic plates to initialize. Upstream default is 10.");

        result.discardPreviewRequested |= ImGui::Checkbox("Override sea level (m)", &session.ui.tectonicUseSeaLevelMeters);
        settooltipifhovered("Enable the --sea-level-m override for the metric coastline threshold.");

        ImGui::BeginDisabled(!session.ui.tectonicUseSeaLevelMeters);
        result.discardPreviewRequested |= ImGui::InputInt("Sea level (m)", &session.ui.tectonicSeaLevelMeters);
        session.ui.tectonicSeaLevelMeters = clamp(session.ui.tectonicSeaLevelMeters, 0, 65535);
        settooltipifhovered("Metric coastline threshold in the range [0, 65535].");
        ImGui::EndDisabled();

        int aggregationoverlapabs = session.ui.tectonicAggregationOverlapAbsolute >= 0
            ? session.ui.tectonicAggregationOverlapAbsolute
            : defaultplatetectonicsaggregationoverlapabs(displayedworld.width() + 1, displayedworld.height() + 1);
        if (ImGui::InputInt("Aggregation overlap abs", &aggregationoverlapabs))
        {
            session.ui.tectonicAggregationOverlapAbsolute = max(0, aggregationoverlapabs);
            result.discardPreviewRequested = true;
        }
        settooltipifhovered("Absolute overlap needed to aggregate continents. Default is max(64, width*height/1000).");

        result.discardPreviewRequested |= ImGui::SliderFloat("Aggregation overlap rel", &session.ui.tectonicAggregationOverlapRelative, 0.0f, 1.0f, "%.2f");
        settooltipifhovered("Relative overlap needed to aggregate continents in the range [0, 1].");
        result.discardPreviewRequested |= ImGui::SliderFloat("Folding ratio", &session.ui.tectonicFoldingRatio, 0.0f, 1.0f, "%.2f");
        settooltipifhovered("Fraction of overlapping continental crust converted into uplift in the range [0, 1].");
        result.discardPreviewRequested |= ImGui::InputInt("Erosion period", &session.ui.tectonicErosionPeriod);
        if (session.ui.tectonicErosionPeriod < 1)
            session.ui.tectonicErosionPeriod = 1;
        settooltipifhovered("Simulation updates between erosion passes. Upstream default is 60.");

        result.discardPreviewRequested |= ImGui::InputFloat("Erosion strength", &session.ui.tectonicErosionStrength, 0.1f, 1.0f, "%.2f");
        if (session.ui.tectonicErosionStrength < 0.0f)
            session.ui.tectonicErosionStrength = 0.0f;
        settooltipifhovered("Erosion strength multiplier. Set 0 to disable erosion.");

        result.discardPreviewRequested |= ImGui::InputFloat("Landmass rotation", &session.ui.tectonicLandmassRotation, 0.05f, 0.25f, "%.2f");
        if (session.ui.tectonicLandmassRotation < 0.0f)
            session.ui.tectonicLandmassRotation = 0.0f;
        settooltipifhovered("Visible crust rotation multiplier. Set 0 to disable landmass rotation.");

        result.discardPreviewRequested |= ImGui::InputFloat("Rotation strength", &session.ui.tectonicRotationStrength, 0.05f, 0.25f, "%.2f");
        if (session.ui.tectonicRotationStrength < 0.0f)
            session.ui.tectonicRotationStrength = 0.0f;
        settooltipifhovered("Angular plate motion multiplier.");

        result.discardPreviewRequested |= ImGui::SliderFloat("Subduction strength", &session.ui.tectonicSubductionStrength, 0.0f, 1.0f, "%.2f");
        settooltipifhovered("Oceanic crust removal strength during subduction in the range [0, 1].");
        result.discardPreviewRequested |= ImGui::InputFloat("Divergent carve", &session.ui.tectonicDivergentCarveStrength, 0.005f, 0.01f, "%.3f");
        if (session.ui.tectonicDivergentCarveStrength < 0.0f)
            session.ui.tectonicDivergentCarveStrength = 0.0f;
        settooltipifhovered("Downward carve strength for regenerated divergent crust. Upstream default is 0.015.");
        result.discardPreviewRequested |= ImGui::InputFloat("Delta time (Myr)", &session.ui.tectonicDeltaTimeMyr, 0.1f, 1.0f, "%.3f");
        if (session.ui.tectonicDeltaTimeMyr <= 0.0f)
            session.ui.tectonicDeltaTimeMyr = 0.001f;
        settooltipifhovered("Geological time advanced by each simulation update, in millions of years.");

        if (session.previewAvailable == false)
            ImGui::TextWrapped("Press Recompute to preview plate tectonics with the current seed and parameters.");
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
    else if (stage.id == GenerationStageId::basin_editor)
    {
        seedcommittedbasincomponents(session);
        vector<GenerationBasinComponent>& visiblecomponents = session.previewAvailable ? session.previewScratch.basinComponents : session.committedScratch.basinComponents;

        ImGui::TextUnformatted("Endorheic basin components");

        if (visiblecomponents.empty())
            ImGui::TextUnformatted("Preview this stage to detect components.");

        for (GenerationBasinComponent& component : visiblecomponents)
        {
            ImGui::PushID(component.id);
            ImGui::Text("Basin %d | cells %d", component.id, component.cellCount);
            ImGui::SameLine();
            if (ImGui::SmallButton(policylabel(component.policy)))
            {
                cyclepolicy(component.policy);

                for (GenerationBasinComponent& committedcomponent : session.committedScratch.basinComponents)
                {
                    if (committedcomponent.id == component.id)
                    {
                        committedcomponent.policy = component.policy;
                        break;
                    }
                }

                for (GenerationBasinComponent& previewcomponent : session.previewScratch.basinComponents)
                {
                    if (previewcomponent.id == component.id)
                    {
                        previewcomponent.policy = component.policy;
                        break;
                    }
                }

                result.displayChanged = true;
            }
            ImGui::PopID();
        }
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
