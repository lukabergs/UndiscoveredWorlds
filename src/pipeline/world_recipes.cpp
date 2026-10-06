#include "pipeline/world_recipes.hpp"

#include <algorithm>
#include <chrono>
#include <set>
#include <stdexcept>

#include "world/planet_fields.hpp"
#include "core/world/field_access.hpp"

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

std::vector<GenerationStageId> stagesfrom(GenerationStageId first)
{
    std::vector<GenerationStageId> steps;
    const auto& stages = getgenerationstages();
    for (std::size_t index = getgenerationstageindex(first); index < stages.size(); index++)
        steps.push_back(stages[index].id);
    return steps;
}

// The erosion stages run only in eroded-world; the other recipes leave them out.
std::vector<GenerationStageId> withouterosion(std::vector<GenerationStageId> steps)
{
    steps.erase(std::remove_if(steps.begin(), steps.end(), [](GenerationStageId stage)
    {
        return stage == GenerationStageId::flow_routing || stage == GenerationStageId::fluvial_erosion;
    }), steps.end());
    return steps;
}

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

WorldRecipe proceduralworldrecipe()
{
    WorldRecipe recipe;
    recipe.name = "procedural-world";
    recipe.description = "Full generation from a seed: tectonics, terrain passes, climate, hydrology, layers, optional society";
    recipe.steps = withouterosion(stagesfrom(GenerationStageId::geosphere_initialization));
    recipe.inputs = worldpropertyinputs();
    return recipe;
}

WorldRecipe importedworldrecipe()
{
    WorldRecipe recipe;
    recipe.name = "imported-world";
    recipe.description = "Imported terrain, then the remaining terrain passes, climate, hydrology and layers";
    recipe.steps = withouterosion(stagesfrom(GenerationStageId::mountain_bases));
    recipe.inputs = importedterraininputs();
    return recipe;
}

WorldRecipe erodedworldrecipe()
{
    WorldRecipe recipe;
    recipe.name = "eroded-world";
    recipe.description = "procedural-world with flow routing and fluvial erosion after the climate coupling and before the fjords; "
                         "the legacy hydrology stages still run after them (Phase A)";
    recipe.steps = stagesfrom(GenerationStageId::geosphere_initialization);
    recipe.inputs = worldpropertyinputs();
    return recipe;
}

WorldRecipe earthclimaterecipe()
{
    WorldRecipe recipe;
    recipe.name = "earth-climate";
    recipe.description = "Climate processes only on imported terrain; Geosphere is frozen";
    for (const GenerationStageDefinition& stage : getgenerationstages())
        if (isclimatefeedbackstage(stage.id) && stage.id != GenerationStageId::terrain_roughness)
            recipe.steps.push_back(stage.id);
    recipe.frozen = { Component::geosphere };
    recipe.inputs = importedterraininputs();
    return recipe;
}

WorldRecipe tectonicsrecipe()
{
    WorldRecipe recipe;
    recipe.name = "tectonics";
    recipe.description = "The tectonics module alone: starting Geosphere, then plate tectonics; the engine runs on the world grid";
    recipe.steps = { GenerationStageId::geosphere_initialization, GenerationStageId::plate_tectonics };
    recipe.inputs = worldpropertyinputs();
    return recipe;
}

bool recipeimportsterrain(const WorldRecipe& recipe)
{
    return recipe.name == "imported-world" || recipe.name == "earth-climate";
}

std::vector<std::string> worldrecipenames()
{
    return { "procedural-world", "imported-world", "eroded-world", "earth-climate", "tectonics" };
}

bool findworldrecipe(std::string_view name, WorldRecipe& recipe)
{
    if (name == "procedural-world")
        recipe = proceduralworldrecipe();
    else if (name == "imported-world")
        recipe = importedworldrecipe();
    else if (name == "eroded-world")
        recipe = erodedworldrecipe();
    else if (name == "earth-climate")
        recipe = earthclimaterecipe();
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

    for (std::size_t index = 0; index < recipe.steps.size(); index++)
    {
        const int step = static_cast<int>(index);
        const WorldProcessDefinition* process = nullptr;
        try
        {
            process = &getworldprocess(recipe.steps[index]);
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

bool findrecipestep(const WorldRecipe& recipe, std::string_view text, std::size_t& index)
{
    if (!text.empty() && text.find_first_not_of("0123456789") == std::string_view::npos)
    {
        index = static_cast<std::size_t>(std::stoul(std::string(text)));
        return index < recipe.steps.size();
    }
    for (std::size_t step = 0; step < recipe.steps.size(); step++)
        if (text == getworldprocess(recipe.steps[step]).id)
        {
            index = step;
            return true;
        }
    return false;
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
    if (!first.empty() && !findrecipestep(recipe, first, steps.first))
        return fail("recipe " + recipe.name + " has no step " + std::string(first));
    if (!last.empty() && !findrecipestep(recipe, last, steps.last))
        return fail("recipe " + recipe.name + " has no step " + std::string(last));
    if (steps.last < steps.first)
        return fail("step " + std::string(last) + " comes before " + std::string(first));
    return true;
}

namespace
{
// State a stage leaves for a later one in GenerationScratch, outside the world
// fields. Until it becomes world state (world-model migration step 5), a run
// cannot start between the two.
struct InRunState
{
    GenerationStageId from;
    GenerationStageId to;
    const char* what;
};

constexpr InRunState kInRunState[] =
{
    { GenerationStageId::mountain_bases, GenerationStageId::rivers_and_basins, "mountain drainage" },
    { GenerationStageId::mountain_bases, GenerationStageId::basin_editor, "mountain drainage" },
    { GenerationStageId::rivers_and_basins, GenerationStageId::lakes, "salt-lake and no-lake masks" },
};
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
        const WorldProcessDefinition& process = getworldprocess(recipe.steps[index]);
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
            for (const InRunState& state : kInRunState)
            {
                if (state.to != recipe.steps[index])
                    continue;
                for (std::size_t earlier = 0; earlier < first; earlier++)
                    if (recipe.steps[earlier] == state.from)
                        issues.push_back({ step, true, std::string(process.id) + " needs the " + state.what + " of " +
                            getworldprocess(state.from).id + " from the same run" });
            }
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
        const GenerationStageId stage = recipe.steps[index];
        const WorldProcessDefinition& process = getworldprocess(stage);
        RecipeStepRecord record;
        record.index = static_cast<int>(index);
        record.stage = stage;
        record.process = process.id;
        record.label = getgenerationstage(stage).label;

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
        if (!executed)
        {
            result.error = record.label + " failed" + (message.empty() ? std::string() : ": " + message);
            result.steps.push_back(record);
            return result;
        }
        if (!message.empty())
            record.violations.push_back("note: " + message);

        if (worldmodel::access::tracking)
        {
            record.reads = worldmodel::access::fieldsread();
            record.tablereads = worldmodel::access::tablesread();
            record.drewrandom = worldmodel::access::drewrandom();
        }

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
        step.seconds = record.seconds;
        for (FieldId write : record.writes)
            step.writes.push_back(worldmodel::qualifiedname(write));
        for (TableId write : record.tablewrites)
            step.writes.push_back(worldmodel::qualifiedname(write));
        steps.push_back(std::move(step));
    }
    return steps;
}
