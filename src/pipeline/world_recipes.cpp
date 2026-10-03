#include "world_recipes.hpp"

#include <algorithm>
#include <chrono>
#include <set>
#include <stdexcept>

#include "planet_fields.hpp"
#include "world/field_access.hpp"

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
        FieldId::seed, FieldId::size_class, FieldId::width_index, FieldId::height_index, FieldId::prograde_rotation,
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
    recipe.steps = stagesfrom(GenerationStageId::geosphere_initialization);
    recipe.inputs = worldpropertyinputs();
    return recipe;
}

WorldRecipe importedworldrecipe()
{
    WorldRecipe recipe;
    recipe.name = "imported-world";
    recipe.description = "Imported terrain, then the remaining terrain passes, climate, hydrology and layers";
    recipe.steps = stagesfrom(GenerationStageId::mountain_bases);
    recipe.inputs = importedterraininputs();
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

std::vector<std::string> worldrecipenames()
{
    return { "procedural-world", "imported-world", "earth-climate" };
}

bool findworldrecipe(std::string_view name, WorldRecipe& recipe)
{
    if (name == "procedural-world")
        recipe = proceduralworldrecipe();
    else if (name == "imported-world")
        recipe = importedworldrecipe();
    else if (name == "earth-climate")
        recipe = earthclimaterecipe();
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

RecipeRunResult runrecipe(const WorldRecipe& recipe, planet& world, GenerationScratch& scratch,
    WorldStageParameters& parameters, const GenerationExecutionContext& context, const RecipeRunOptions& options)
{
    RecipeRunResult result;
    for (const RecipeIssue& issue : validaterecipe(recipe))
    {
        if (issue.error)
        {
            result.error = "Recipe " + recipe.name + " step " + std::to_string(issue.step) + ": " + issue.message;
            return result;
        }
    }

    WorldHashes before;
    if (options.auditAccess)
        before = hashworld(world);

    for (std::size_t index = 0; index < recipe.steps.size(); index++)
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
        const bool executed = executegenerationstage(stage, world, scratch, parameters, context, &message);
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
