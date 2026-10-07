#include "pipeline/world_recipes.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <map>
#include <set>
#include <stdexcept>

#include "world/planet_fields.hpp"
#include "core/world/field_access.hpp"
#include "pipeline/recipe_parameters.hpp"

using worldmodel::Component;
using worldmodel::FieldId;
using worldmodel::TableId;

namespace
{
// World properties set before any recipe runs (seed, planet body, climate
// parameters, grid and datum).
std::vector<FieldId> worldpropertyinputs()
{
    return {
        FieldId::seed, FieldId::size_class, FieldId::radius, FieldId::width_index, FieldId::height_index, FieldId::prograde_rotation,
        FieldId::axial_tilt, FieldId::eccentricity, FieldId::perihelion, FieldId::gravity, FieldId::lunar_strength,
        FieldId::lapse_rate, FieldId::north_polar_adjust, FieldId::south_polar_adjust, FieldId::average_temperature,
        FieldId::north_polar_temperature, FieldId::south_polar_temperature, FieldId::equatorial_temperature,
        FieldId::water_pickup, FieldId::river_factor, FieldId::river_land_reduce, FieldId::estuary_limit,
        FieldId::glacial_temperature, FieldId::glacier_temperature, FieldId::mountain_reduce, FieldId::climate_type_count,
        FieldId::max_elevation, FieldId::sea_level
    };
}

// Imported terrain supplies the whole Geosphere state.
std::vector<FieldId> importedterraininputs()
{
    std::vector<FieldId> inputs = worldpropertyinputs();
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        if (worldmodel::fieldspec(id).component == Component::geosphere && std::find(inputs.begin(), inputs.end(), id) == inputs.end())
            inputs.push_back(id);
    }
    return inputs;
}

// The stages from first to last in registry order (inclusive).
std::vector<GenerationStageId> stagesbetween(GenerationStageId first, GenerationStageId last)
{
    std::vector<GenerationStageId> steps;
    const auto& stages = getgenerationstages();
    for (std::size_t index = getgenerationstageindex(first); index <= getgenerationstageindex(last); index++)
        steps.push_back(stages[index].id);
    return steps;
}

// Steps of the given stages, each labelled with a suffix (none: the process id) and the same overrides.
void appendsteps(std::vector<RecipeStep>& steps, const std::vector<GenerationStageId>& stages,
    std::string_view suffix = {}, const std::vector<RecipeOverride>& overrides = {})
{
    for (GenerationStageId stage : stages)
        steps.push_back(recipestep(stage, suffix, overrides));
}

std::vector<RecipeStep> labelledsteps(const std::vector<GenerationStageId>& stages)
{
    std::vector<RecipeStep> steps;
    appendsteps(steps, stages);
    return steps;
}

// Climate-context settings a step can override (cli/climate_options.cpp sets them for the run).
constexpr std::string_view kClimateResolution = "climate.resolution";
constexpr std::string_view kClimateCouplingIterations = "climate.coupling_iterations";

bool readinteger(const RecipeOverride& entry, int& value, std::string& error)
{
    std::size_t parsed = 0;
    try
    {
        value = std::stoi(entry.value, &parsed);
    }
    catch (const std::exception&)
    {
        parsed = 0;
    }
    if (parsed == 0 || parsed != entry.value.size())
    {
        error = entry.name + " takes a whole number, not " + entry.value;
        return false;
    }
    return true;
}

// Checks an override's name and value with the bounds of the matching command-line
// option (cli/climate_options.cpp) or the parameter's type; `number` receives a
// climate-context value. Changes nothing.
bool checkoverride(const RecipeOverride& entry, std::string& error, int* number = nullptr)
{
    int value = 0;
    if (entry.name == kClimateResolution)
    {
        if (!readinteger(entry, value, error))
            return false;
        if (value < 8 || value % 4 != 0 || value > ARRAYWIDTH / 2)
        {
            error = entry.name + " must be a multiple of 4 between 8 and " + std::to_string(ARRAYWIDTH / 2) + ", not " + entry.value;
            return false;
        }
    }
    else if (entry.name == kClimateCouplingIterations)
    {
        if (!readinteger(entry, value, error))
            return false;
        if (value < 1 || value > 100)
        {
            error = entry.name + " must be between 1 and 100, not " + entry.value;
            return false;
        }
    }
    else
        return checkrecipeparameter(entry.name, entry.value, &error);
    if (number != nullptr)
        *number = value;
    return true;
}

// Applies one step's overrides and puts back what they replaced.
class StepOverrides
{
public:
    explicit StepOverrides(planet& world) : world_(world) {}
    ~StepOverrides() { restore(); }
    StepOverrides(const StepOverrides&) = delete;
    StepOverrides& operator=(const StepOverrides&) = delete;

    bool apply(const std::vector<RecipeOverride>& overrides, std::vector<RecipeOverride>& applied, std::string& error)
    {
        for (const RecipeOverride& entry : overrides)
        {
            int number = 0;
            if (!checkoverride(entry, error, &number))
                return false;
            if (entry.name == kClimateResolution)
            {
                climatecontext::Context& climate = world_.climatesimulation();
                const int columns = number;
                if (!savedresolution_)
                    resolution_ = climate.resolution;
                savedresolution_ = true;
                // As --climate-resolution without --climate-hydrology-resolution.
                climate.resolution = { columns, columns, 2 * columns, columns / 2 };
                applied.push_back({ entry.name, std::to_string(columns) });
            }
            else if (entry.name == kClimateCouplingIterations)
            {
                climatecontext::Context& climate = world_.climatesimulation();
                const int iterations = number;
                if (!savedcoupling_)
                    coupling_ = climate.couplingIterationLimit;
                savedcoupling_ = true;
                climate.couplingIterationLimit = iterations;
                applied.push_back({ entry.name, std::to_string(iterations) });
            }
            else
            {
                std::string previous;
                getrecipeparameter(entry.name, previous);
                if (!setrecipeparameter(entry.name + "=" + entry.value, &error))
                    return false;
                parameters_.push_back({ entry.name, previous });
                std::string value;
                getrecipeparameter(entry.name, value);
                applied.push_back({ entry.name, value });
            }
        }
        if (!parameters_.empty() && !validaterecipeparameters(&error))
            return false;
        return true;
    }

    void restore()
    {
        // In reverse, so that a name overridden twice gets its first value back.
        for (auto entry = parameters_.rbegin(); entry != parameters_.rend(); ++entry)
            setrecipeparameter(entry->name + "=" + entry->value);
        parameters_.clear();
        if (savedresolution_)
            world_.climatesimulation().resolution = resolution_;
        if (savedcoupling_)
            world_.climatesimulation().couplingIterationLimit = coupling_;
        savedresolution_ = savedcoupling_ = false;
    }

private:
    planet& world_;
    std::vector<RecipeOverride> parameters_;      // values before the step
    climatecontext::Resolution resolution_;
    int coupling_ = 0;
    bool savedresolution_ = false;
    bool savedcoupling_ = false;
};

bool isfrozen(const WorldRecipe& recipe, Component component)
{
    return std::find(recipe.frozen.begin(), recipe.frozen.end(), component) != recipe.frozen.end();
}

struct WorldHashes
{
    std::vector<std::uint64_t> fields;
    std::vector<std::uint64_t> tables;
};

WorldHashes hashworld(const planet& world)
{
    WorldHashes hashes;
    hashes.fields.resize(worldmodel::kFieldCount);
    hashes.tables.resize(worldmodel::kTableCount);
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
        hashes.fields[index] = hashplanetfield(world, static_cast<FieldId>(index));
    for (std::size_t index = 0; index < worldmodel::kTableCount; index++)
        hashes.tables[index] = hashplanettable(world, static_cast<TableId>(index));
    return hashes;
}
}

RecipeStep recipestep(GenerationStageId stage, std::string_view suffix, std::vector<RecipeOverride> overrides)
{
    RecipeStep step;
    step.stage = stage;
    step.label = getworldprocess(stage).id;
    if (!suffix.empty())
        step.label += "@" + std::string(suffix);
    step.overrides = std::move(overrides);
    return step;
}

bool recipehasstage(const WorldRecipe& recipe, GenerationStageId stage)
{
    return std::any_of(recipe.steps.begin(), recipe.steps.end(), [&](const RecipeStep& step) { return step.stage == stage; });
}

WorldRecipe proceduralworldrecipe()
{
    WorldRecipe recipe;
    recipe.name = "procedural-world";
    recipe.description = "Full generation from a seed: tectonics and the terrain passes; a coarse climate (32 columns, 2 "
                         "coupling iterations) for the fluvial erosion, which routes itself; the flow "
                         "routing and lake levels on the eroded terrain; the full climate on it at the run's climate "
                         "settings, with those lakes and pans; the flow routing again, rivers from it (hydrology.rivers) "
                         "and lakes from a water balance in its depressions (hydrology.lake_levels); the remaining stages, "
                         "layers.finalize once at the end (society with --social)";
    // The full climate recomputes everything on the eroded terrain except one deliberate
    // warm start: climate.pressure and climate.winds read the sea-surface temperatures the
    // coarse climate left in the world (pressuresurfacetemperature, climate_coupling_internal.hpp)
    // until climate.ocean_currents rewrites them.
    using S = GenerationStageId;
    const std::vector<GenerationStageId> climate = stagesbetween(S::global_temperature, S::rainfall);
    appendsteps(recipe.steps, stagesbetween(S::geosphere_initialization, S::terrain_roughness));
    appendsteps(recipe.steps, climate, "coarse", { { "climate.resolution", "32" }, { "climate.coupling_iterations", "2" } });
    appendsteps(recipe.steps, { S::fluvial_erosion });
    // The lakes and salt pans of the eroded terrain before the full climate, so that it sees them
    // (erosion plan H3: an enclosed below-sea basin is a lake or a pan, not sea); both run again
    // after it. A repeated process needs a label on every step.
    appendsteps(recipe.steps, { S::flow_routing, S::lake_levels }, "before-climate");
    appendsteps(recipe.steps, climate, "full");
    appendsteps(recipe.steps, { S::flow_routing }, "after-climate");
    appendsteps(recipe.steps, { S::rivers });
    appendsteps(recipe.steps, { S::lake_levels }, "after-climate");
    appendsteps(recipe.steps, { S::mountain_temperature_lapse, S::climates_and_biomes, S::arid_features,
                                S::deltas_wetlands_and_roughness, S::finalize_layers });
    recipe.inputs = worldpropertyinputs();
    // hydrology.lake_levels reads the imported lake mask: all zero in a procedural world.
    recipe.inputs.push_back(FieldId::imported_lake);
    return recipe;
}

WorldRecipe importedworldrecipe()
{
    WorldRecipe recipe;
    recipe.name = "imported-world";
    recipe.description = "Imported terrain, then the remaining terrain passes, climate, the flow routing, "
                         "rivers from it (hydrology.rivers) and lakes from a water balance in its depressions "
                         "(hydrology.lake_levels), and layers; the imported terrain is not eroded";
    using S = GenerationStageId;
    // terrain.depression_fill is left out: the lake balance owns the depressions of the imported terrain.
    std::vector<GenerationStageId> terrainandclimate = stagesbetween(S::mountain_bases, S::rainfall);
    terrainandclimate.erase(std::remove(terrainandclimate.begin(), terrainandclimate.end(), S::depression_fill),
                            terrainandclimate.end());
    appendsteps(recipe.steps, terrainandclimate);
    // The routing on the final terrain, then the rivers and lakes that read it.
    appendsteps(recipe.steps, { S::flow_routing, S::rivers, S::lake_levels });
    appendsteps(recipe.steps, { S::mountain_temperature_lapse, S::climates_and_biomes, S::arid_features,
                                S::deltas_wetlands_and_roughness, S::finalize_layers });
    recipe.inputs = importedterraininputs();
    // hydrology.lake_levels reads the imported lake mask (--import-lakes), written by the import before the recipe.
    if (std::find(recipe.inputs.begin(), recipe.inputs.end(), FieldId::imported_lake) == recipe.inputs.end())
        recipe.inputs.push_back(FieldId::imported_lake);
    return recipe;
}

WorldRecipe earthclimaterecipe()
{
    WorldRecipe recipe;
    recipe.name = "earth-climate";
    recipe.description = "Climate processes only on imported terrain; Geosphere is frozen";
    for (const GenerationStageDefinition& stage : getgenerationstages())
        if (isclimatefeedbackstage(stage.id) && stage.id != GenerationStageId::terrain_roughness)
            recipe.steps.push_back(recipestep(stage.id));
    recipe.frozen = { Component::geosphere };
    recipe.inputs = importedterraininputs();
    return recipe;
}

WorldRecipe earthhydrologyrecipe()
{
    WorldRecipe recipe = earthclimaterecipe();
    recipe.name = "earth-hydrology";
    recipe.description = "earth-climate with the flow routing and the lake levels after the rainfall, on the imported "
                         "terrain and its imported lakes (--import-lakes); Geosphere is frozen";
    const auto rainfall = std::find_if(recipe.steps.begin(), recipe.steps.end(),
        [](const RecipeStep& step) { return step.stage == GenerationStageId::rainfall; });
    const auto routing = recipe.steps.insert(rainfall == recipe.steps.end() ? rainfall : rainfall + 1,
        recipestep(GenerationStageId::flow_routing));
    recipe.steps.insert(routing + 1, recipestep(GenerationStageId::lake_levels));
    // hydrology.lake_levels reads the imported lake mask, which the import writes before the recipe.
    recipe.inputs.push_back(FieldId::imported_lake);
    return recipe;
}

WorldRecipe tectonicsrecipe()
{
    WorldRecipe recipe;
    recipe.name = "tectonics";
    recipe.description = "The tectonics module alone: starting Geosphere, then plate tectonics; the engine runs on the world grid";
    recipe.steps = labelledsteps({ GenerationStageId::geosphere_initialization, GenerationStageId::plate_tectonics });
    recipe.inputs = worldpropertyinputs();
    return recipe;
}

bool recipeimportsterrain(const WorldRecipe& recipe)
{
    return recipe.name == "imported-world" || recipe.name == "earth-climate" || recipe.name == "earth-hydrology";
}

std::vector<std::string> worldrecipenames()
{
    return { "procedural-world", "imported-world", "earth-climate", "earth-hydrology", "tectonics" };
}

bool findworldrecipe(std::string_view name, WorldRecipe& recipe)
{
    if (name == "procedural-world")
        recipe = proceduralworldrecipe();
    else if (name == "imported-world")
        recipe = importedworldrecipe();
    else if (name == "earth-climate")
        recipe = earthclimaterecipe();
    else if (name == "earth-hydrology")
        recipe = earthhydrologyrecipe();
    else if (name == "tectonics")
        recipe = tectonicsrecipe();
    else
        return false;
    return true;
}

std::vector<RecipeIssue> validaterecipe(const WorldRecipe& recipe)
{
    std::vector<RecipeIssue> issues;
    std::set<FieldId> available(recipe.inputs.begin(), recipe.inputs.end());
    std::set<TableId> availabletables(recipe.inputtables.begin(), recipe.inputtables.end());

    // Labels are unique, and a label that is a process id belongs to that process's only step.
    std::map<GenerationStageId, int> occurrences;
    for (const RecipeStep& step : recipe.steps)
        occurrences[step.stage]++;
    std::set<std::string> labels;
    for (std::size_t index = 0; index < recipe.steps.size(); index++)
    {
        const RecipeStep& step = recipe.steps[index];
        const int position = static_cast<int>(index);
        if (step.label.empty() || !labels.insert(step.label).second)
            issues.push_back({ position, true, "label '" + step.label + "' is empty or not unique" });
        if (const WorldProcessDefinition* named = findworldprocess(step.label))
            if (named->stage != step.stage || occurrences[step.stage] > 1)
                issues.push_back({ position, true, "label " + step.label + " is a process id; a repeated process needs a suffix on every step" });
        // Override values are checked here, before any step runs.
        for (const RecipeOverride& entry : step.overrides)
        {
            std::string problem;
            if (!checkoverride(entry, problem))
                issues.push_back({ position, true, step.label + " override: " + problem });
        }
    }

    for (std::size_t index = 0; index < recipe.steps.size(); index++)
    {
        const int step = static_cast<int>(index);
        const WorldProcessDefinition* process = nullptr;
        try
        {
            process = &getworldprocess(recipe.steps[index].stage);
        }
        catch (const std::exception&)
        {
            issues.push_back({ step, true, "Stage has no process definition" });
            continue;
        }

        for (FieldId read : process->reads)
            if (available.count(read) == 0)
                issues.push_back({ step, false, std::string(process->id) + " reads " + worldmodel::qualifiedname(read) + " before any step writes it" });
        for (TableId read : process->tablereads)
            if (availabletables.count(read) == 0)
                issues.push_back({ step, false, std::string(process->id) + " reads table " + worldmodel::qualifiedname(read) + " before any step writes it" });

        for (FieldId write : process->writes)
        {
            if (isfrozen(recipe, worldmodel::fieldspec(write).component))
                issues.push_back({ step, true, std::string(process->id) + " writes " + worldmodel::qualifiedname(write) + " in a frozen component" });
            available.insert(write);
        }
        for (TableId write : process->tablewrites)
        {
            if (isfrozen(recipe, worldmodel::tablespec(write).component))
                issues.push_back({ step, true, std::string(process->id) + " writes table " + worldmodel::qualifiedname(write) + " in a frozen component" });
            availabletables.insert(write);
        }
    }
    return issues;
}

const WorldProcessDefinition* findtimedevolver(const WorldRecipe& recipe, std::string_view label, std::string* reason)
{
    const WorldProcessDefinition* process = findworldprocess(label);
    for (std::size_t step = 0; process == nullptr && step < recipe.steps.size(); step++)
        if (label == recipe.steps[step].label)
            process = &getworldprocess(recipe.steps[step].stage);
    if (process != nullptr && process->time.source != nullptr)
        return process;
    if (reason != nullptr)
        *reason = std::string(label) + " is neither an evolver of the time model (" + timedevolverids() +
            ") nor the label of one of its steps in recipe " + recipe.name;
    return nullptr;
}

bool findrecipestep(const WorldRecipe& recipe, std::string_view text, std::size_t& index, std::string* reason)
{
    const auto fail = [&](const std::string& message)
    {
        if (reason != nullptr)
            *reason = message;
        return false;
    };
    const std::string missing = "recipe " + recipe.name + " has no step " + std::string(text);
    if (!text.empty() && text.find_first_not_of("0123456789") == std::string_view::npos)
    {
        index = static_cast<std::size_t>(std::stoul(std::string(text)));
        return index < recipe.steps.size() || fail(missing);
    }
    for (std::size_t step = 0; step < recipe.steps.size(); step++)
        if (text == recipe.steps[step].label)
        {
            index = step;
            return true;
        }
    std::vector<std::size_t> matches;
    for (std::size_t step = 0; step < recipe.steps.size(); step++)
        if (text == getworldprocess(recipe.steps[step].stage).id)
            matches.push_back(step);
    if (matches.empty())
        return fail(missing);
    if (matches.size() > 1)
    {
        std::string labels;
        for (std::size_t match : matches)
            labels += (labels.empty() ? "" : ", ") + recipe.steps[match].label;
        return fail("process " + std::string(text) + " occurs " + std::to_string(matches.size()) + " times in recipe " +
            recipe.name + "; name one step by its label: " + labels);
    }
    index = matches.front();
    return true;
}

std::string recipesteppathname(std::string_view label)
{
    std::string name(label);
    for (char& character : name)
    {
        const unsigned char code = static_cast<unsigned char>(character);
        if (!std::isalnum(code) && character != '.' && character != '_' && character != '@' && character != '-')
            character = '_';
    }
    return name;
}

bool parserecipesteps(const WorldRecipe& recipe, std::string_view text, RecipeSteps& steps, std::string* reason)
{
    const auto fail = [&](const std::string& message)
    {
        if (reason != nullptr)
            *reason = message;
        return false;
    };
    if (recipe.steps.empty())
        return fail("recipe " + recipe.name + " has no steps");
    const std::size_t separator = text.find("..");
    const std::string_view first = text.substr(0, separator);
    const std::string_view last = separator == std::string_view::npos ? first : text.substr(separator + 2);
    steps.first = 0;
    steps.last = recipe.steps.size() - 1;
    if (!first.empty() && !findrecipestep(recipe, first, steps.first, reason))
        return false;
    if (!last.empty() && !findrecipestep(recipe, last, steps.last, reason))
        return false;
    if (steps.last < steps.first)
        return fail("step " + std::string(last) + " comes before " + std::string(first));
    return true;
}

std::vector<RecipeIssue> validaterecipesteps(const WorldRecipe& recipe, const RecipeRunOptions& options)
{
    std::vector<RecipeIssue> issues;
    if (recipe.steps.empty())
        return issues;
    const std::size_t first = options.selectSteps ? options.steps.first : 0;
    const std::size_t last = options.selectSteps ? options.steps.last : recipe.steps.size() - 1;

    // What the full recipe has before each step, and what this run has.
    std::set<FieldId> expected(recipe.inputs.begin(), recipe.inputs.end());
    std::set<TableId> expectedtables(recipe.inputtables.begin(), recipe.inputtables.end());
    std::set<FieldId> present = options.fromBundle ? std::set<FieldId>(options.loadedfields.begin(), options.loadedfields.end()) : expected;
    std::set<TableId> presenttables = options.fromBundle ? std::set<TableId>(options.loadedtables.begin(), options.loadedtables.end()) : expectedtables;
    if (options.fromBundle)
    {
        // A bundle is a whole world: whatever it supplied was written before.
        expected.insert(present.begin(), present.end());
        expectedtables.insert(presenttables.begin(), presenttables.end());
    }

    for (std::size_t index = 0; index <= last && index < recipe.steps.size(); index++)
    {
        const WorldProcessDefinition& process = getworldprocess(recipe.steps[index].stage);
        const int step = static_cast<int>(index);
        if (index >= first)
        {
            for (FieldId read : process.reads)
                if (expected.count(read) != 0 && present.count(read) == 0)
                    issues.push_back({ step, true, std::string(process.id) + " reads " + worldmodel::qualifiedname(read) +
                        ", which an earlier step outside this run provides and a bundle cannot load" });
            for (TableId read : process.tablereads)
                if (expectedtables.count(read) != 0 && presenttables.count(read) == 0)
                    issues.push_back({ step, true, std::string(process.id) + " reads table " + worldmodel::qualifiedname(read) +
                        ", which an earlier step outside this run provides and a bundle cannot load" });
        }
        expected.insert(process.writes.begin(), process.writes.end());
        expectedtables.insert(process.tablewrites.begin(), process.tablewrites.end());
        if (index >= first)
        {
            present.insert(process.writes.begin(), process.writes.end());
            presenttables.insert(process.tablewrites.begin(), process.tablewrites.end());
        }
    }
    return issues;
}

RecipeRunResult runrecipe(const WorldRecipe& recipe, planet& world, GenerationScratch& scratch,
    WorldStageParameters& parameters, const GenerationExecutionContext& context, const RecipeRunOptions& options)
{
    RecipeRunResult result;
    std::vector<RecipeIssue> issues = validaterecipe(recipe);
    if (options.selectSteps || options.fromBundle)
    {
        const std::vector<RecipeIssue> stepissues = validaterecipesteps(recipe, options);
        issues.insert(issues.end(), stepissues.begin(), stepissues.end());
    }
    for (const RecipeIssue& issue : issues)
    {
        if (issue.error)
        {
            result.error = "Recipe " + recipe.name + " step " + std::to_string(issue.step) + ": " + issue.message;
            return result;
        }
    }
    result.provenance.fields = options.sourceprovenance.fields;
    result.provenance.tables = options.sourceprovenance.tables;

    WorldHashes before;
    if (options.auditAccess)
        before = hashworld(world);

    const std::size_t firststep = options.selectSteps ? options.steps.first : 0;
    const std::size_t laststep = options.selectSteps ? options.steps.last + 1 : recipe.steps.size();
    for (std::size_t index = firststep; index < laststep && index < recipe.steps.size(); index++)
    {
        const RecipeStep& step = recipe.steps[index];
        const GenerationStageId stage = step.stage;
        const WorldProcessDefinition& process = getworldprocess(stage);
        RecipeStepRecord record;
        record.index = static_cast<int>(index);
        record.stage = stage;
        record.process = process.id;
        record.label = step.label;

        // The step's own settings, set before the read tracking starts and restored before anything else runs.
        StepOverrides overrides(world);
        std::string overrideerror;
        if (!overrides.apply(step.overrides, record.overrides, overrideerror))
        {
            result.error = record.label + " failed: override: " + overrideerror;
            result.steps.push_back(record);
            return result;
        }

        worldmodel::access::resetreads();
        std::string message;
        const auto start = std::chrono::steady_clock::now();
        bool executed = false;
        try
        {
            executed = executegenerationstage(stage, world, scratch, parameters, context, &message);
        }
        catch (const std::exception& error)
        {
            // A process rejects its parameters or a checkpoint fails: the run stops here.
            message = error.what();
        }
        record.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();

        if (worldmodel::access::tracking)
        {
            record.reads = worldmodel::access::fieldsread();
            record.tablereads = worldmodel::access::tablesread();
            record.drewrandom = worldmodel::access::drewrandom();
        }
        overrides.restore();

        if (!executed)
        {
            result.error = record.label + " failed" + (message.empty() ? std::string() : ": " + message);
            result.steps.push_back(record);
            return result;
        }
        if (!message.empty())
            record.violations.push_back("note: " + message);

        if (options.auditAccess)
        {
            const WorldHashes after = hashworld(world);
            record.measured = true;
            for (std::size_t field = 0; field < worldmodel::kFieldCount; field++)
                if (after.fields[field] != before.fields[field])
                {
                    record.writes.push_back(static_cast<FieldId>(field));
                    record.writehashes.push_back(after.fields[field]);
                }
            for (std::size_t table = 0; table < worldmodel::kTableCount; table++)
                if (after.tables[table] != before.tables[table])
                    record.tablewrites.push_back(static_cast<TableId>(table));
            before = after;

            for (FieldId write : record.writes)
            {
                const worldmodel::FieldSpec& spec = worldmodel::fieldspec(write);
                if (spec.component != Component::presentation && !declareswrite(process, write))
                    record.violations.push_back("undeclared write " + worldmodel::qualifiedname(write));
                if (isfrozen(recipe, spec.component))
                    record.violations.push_back("write to frozen " + worldmodel::qualifiedname(write));
            }
            for (TableId write : record.tablewrites)
            {
                if (!declareswrite(process, write))
                    record.violations.push_back("undeclared table write " + worldmodel::qualifiedname(write));
                if (isfrozen(recipe, worldmodel::tablespec(write).component))
                    record.violations.push_back("write to frozen table " + worldmodel::qualifiedname(write));
            }
        }
        else
        {
            record.writes = process.writes;
            record.tablewrites = process.tablewrites;
        }

        for (FieldId write : record.writes)
            result.provenance.fields[write] = record.label;
        for (TableId write : record.tablewrites)
            result.provenance.tables[write] = record.label;
        for (const std::string& violation : record.violations)
            if (violation.rfind("note: ", 0) != 0)
                result.violations++;
        result.steps.push_back(std::move(record));
        if (options.afterstep && !options.afterstep(result, result.steps.back(), result.error))
            return result;
    }

    result.completed = true;
    return result;
}

std::vector<runbundle::StepRecord> bundlesteps(const RecipeRunResult& result)
{
    std::vector<runbundle::StepRecord> steps;
    for (const RecipeStepRecord& record : result.steps)
    {
        runbundle::StepRecord step;
        step.index = record.index;
        step.process = record.process;
        step.label = record.label;
        for (const RecipeOverride& entry : record.overrides)
            step.overrides.emplace_back(entry.name, entry.value);
        step.seconds = record.seconds;
        for (FieldId write : record.writes)
            step.writes.push_back(worldmodel::qualifiedname(write));
        for (TableId write : record.tablewrites)
            step.writes.push_back(worldmodel::qualifiedname(write));
        steps.push_back(std::move(step));
    }
    return steps;
}
