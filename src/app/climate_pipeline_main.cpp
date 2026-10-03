// Console entry point for climate on imported terrain. Runs the earth-climate
// recipe on a climateworld with no planet: terrain is imported straight into the
// climate inputs, and the executable links only climate, I/O and validation.

#include <windows.h>
#include <shellapi.h>

#include <chrono>
#include <exception>
#include <filesystem>
#include <memory>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "app_environment.hpp"
#include "climate_artifacts.hpp"
#include "climate_benchmark_outputs.hpp"
#include "climate_options.hpp"
#include "climate_stages.hpp"
#include "climate_world.hpp"
#include "climate_world_fields.hpp"
#include "generation_progress.hpp"
#include "generation_stage_registry.hpp"
#include "generation_tuning.hpp"
#include "raster_imports.hpp"
#include "run_bundle.hpp"
#include "world_generation_debug.hpp"
#include "world_processes.hpp"

using namespace std;

void exportclimatevalidationreport(climateworld& world);
bool recordclimatebenchmarkrun(climateworld& world, const string& information, const climatebenchmarkmapselection& maps, int* runid);

namespace
{
// Same default grid as world_pipeline climate-benchmark.
constexpr int kDefaultWidth = 512;
// Imported terrain uses this elevation range, as in runimportedworldrecipe.
constexpr int kImportedMaxElevation = 200000;

struct ClimatePipelineOptions
{
    bool benchmark = false;
    long seed = 1;
    int width = kDefaultWidth;
    string importLandPath;
    string importSeaPath;
    string bundlePath;
    string benchmarkInformation;
    bool mapsSpecified = false;
    climatebenchmarkmapselection maps = climatebenchmarkmapselection::defaultselection();
    ClimateRunOptions climate;
    vector<string> arguments;
};

void printusage()
{
    cout << "climate_pipeline <command> [options]\n";
    cout << "Commands:\n";
    cout << "  run         Climate on imported terrain (default: the Earth benchmark terrain); unregistered.\n";
    cout << "  benchmark   As run, then record a numbered Earth climate benchmark run.\n";
    cout << "Options:\n";
    cout << "  --seed <number>         Seed of the climate noise streams (default 1).\n";
    cout << "  --resolution <width>    Even world width; height is width/2 (default 512).\n";
    cout << "  --import-land <tiff>    Land height map; default from configs/app.env for this width.\n";
    cout << "  --import-sea <tiff>     Sea depth map; default from configs/app.env for this width.\n";
    cout << "  --bundle <directory>    Write the climate fields and inputs as a run bundle.\n";
    cout << "  --benchmark-info <text> Change summary stored with a benchmark run.\n";
    cout << "  --map <id[,id...]>      Benchmark maps (as world_pipeline); --map none writes none.\n";
    printclimateoptionsusage(cout);
}

bool parse(const vector<string>& arguments, ClimatePipelineOptions& options)
{
    options.arguments = arguments;
    for (size_t index = 0; index < arguments.size(); index++)
    {
        const string& argument = arguments[index];
        if (index + 1 >= arguments.size())
        {
            cerr << "Missing value for " << argument << '\n';
            return false;
        }
        const string value = arguments[++index];
        string reason;
        const ClimateOptionResult result = parseclimateoption(argument, value, options.climate, &reason);
        if (result == ClimateOptionResult::invalid)
        {
            cerr << "Invalid value for " << argument << ": " << reason << '\n';
            return false;
        }
        if (result == ClimateOptionResult::parsed)
            continue;
        try
        {
            if (argument == "--seed")
                options.seed = stol(value);
            else if (argument == "--resolution")
            {
                options.width = stoi(value);
                if (options.width < 2 || options.width % 2 != 0 || options.width > ARRAYWIDTH)
                    throw invalid_argument("resolution");
            }
            else if (argument == "--import-land")
                options.importLandPath = value;
            else if (argument == "--import-sea")
                options.importSeaPath = value;
            else if (argument == "--bundle")
                options.bundlePath = value;
            else if (argument == "--benchmark-info")
                options.benchmarkInformation = value;
            else if (argument == "--map" || argument == "--maps")
            {
                if (!options.mapsSpecified)
                {
                    options.maps.clear();
                    options.mapsSpecified = true;
                }
                string failure;
                if (!addclimatebenchmarkmapargument(options.maps, value, &failure))
                {
                    cerr << failure << '\n';
                    return false;
                }
            }
            else
            {
                cerr << "Unknown argument: " << argument << '\n';
                return false;
            }
        }
        catch (const exception&)
        {
            cerr << "Invalid value for " << argument << ": " << value << '\n';
            return false;
        }
    }
    return true;
}

// Terrain owned by this run; climateinputs points into it.
struct ImportedTerrain
{
    ImportedTerrain(int columns, int rows) : nom(columns, rows), zeroshort(columns, rows), zeroint(columns, rows) {}
    worldgrid<int> nom;
    worldgrid<short> zeroshort;
    worldgrid<int> zeroint;
};

climateinputs earthinputs(const ClimatePipelineOptions& options, const ImportedTerrain& terrain)
{
    using namespace tuning::worlddefaults;
    climateinputs inputs;
    inputs.seed = options.seed;
    inputs.rotation = rotation;
    inputs.tilt = tilt;
    inputs.eccentricity = eccentricity;
    inputs.perihelion = perihelion;
    inputs.gravity = gravity;
    inputs.lunar = lunar;
    inputs.tempdecrease = tempDecrease;
    inputs.northpolaradjust = northPolarAdjust;
    inputs.southpolaradjust = southPolarAdjust;
    inputs.averagetemp = averageTemp;
    inputs.waterpickup = waterPickup;
    inputs.maxelevation = kImportedMaxElevation;
    inputs.sealevel = seaLevel;
    inputs.rows = terrain.nom.rows();
    inputs.nom = &terrain.nom[0][0];
    inputs.mountainheights = &terrain.zeroshort[0][0];
    inputs.extraelev = &terrain.zeroshort[0][0];
    inputs.craterrims = &terrain.zeroshort[0][0];
    inputs.cratercentres = &terrain.zeroshort[0][0];
    inputs.volcano = &terrain.zeroshort[0][0];
    inputs.oceanridgeheight = &terrain.zeroshort[0][0];
    inputs.lakesurface = &terrain.zeroint[0][0];
    inputs.riftlakesurface = &terrain.zeroint[0][0];
    inputs.special = &terrain.zeroshort[0][0];
    return inputs;
}

void writebundle(const climateworld& world, const ClimatePipelineOptions& options, const vector<runbundle::StepRecord>& steps)
{
    // Declared provenance: the last step whose declared writes include the field.
    map<string, string> writtenby;
    for (const runbundle::StepRecord& step : steps)
        for (const string& field : step.writes)
            writtenby[field] = step.label;

    runbundle::Writer writer(options.bundlePath);
    for (worldmodel::FieldId id : climateworldfields(true))
    {
        worldmodel::FieldData data;
        if (!readclimateworldfield(world, id, data, true))
            continue;
        const worldmodel::FieldSpec& spec = worldmodel::fieldspec(id);
        runbundle::FieldEntry entry;
        entry.id = worldmodel::qualifiedname(id);
        entry.component = worldmodel::componentname(spec.component);
        entry.units = spec.units;
        entry.description = spec.description;
        entry.layout = worldmodel::layoutname(spec.layout);
        entry.dtype = worldmodel::elementname(spec.element);
        entry.shape = data.shape;
        entry.data = data.bytes.data();
        entry.bytes = data.bytes.size();
        if (const auto found = writtenby.find(entry.id); found != writtenby.end())
            entry.writtenby = found->second;
        writer.addfield(entry);
    }
    for (worldmodel::TableId id : { worldmodel::TableId::weather_state, worldmodel::TableId::atmosphere_levels })
    {
        runbundle::TableEntry table;
        table.id = worldmodel::qualifiedname(id);
        table.component = worldmodel::componentname(worldmodel::tablespec(id).component);
        table.description = worldmodel::tablespec(id).description;
        readclimateworldtable(world, id, table.text, table.extension);
        writer.addtable(table);
    }

    runbundle::Metadata metadata;
    metadata.seed = options.seed;
    metadata.columns = world.width() + 1;
    metadata.rows = world.height() + 1;
    metadata.recipe = "earth-climate";
    metadata.command = options.arguments;
    metadata.steps = steps;
    metadata.properties.push_back({ "executable", "climate_pipeline" });
    metadata.properties.push_back({ "executable_built", string(__DATE__) + " " + __TIME__ });
    metadata.properties.push_back({ "access", "declared" });
    writer.finish(metadata);
}
}

int main()
{
    int argc = 0;
    LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
    vector<string> arguments;
    for (int index = 1; argv != nullptr && index < argc; index++)
        arguments.push_back(filesystem::path(argv[index]).string());
    if (argv != nullptr)
        LocalFree(argv);

    if (arguments.empty() || arguments[0] == "--help" || arguments[0] == "-h" || arguments[0] == "help")
    {
        printusage();
        return arguments.empty() ? 1 : 0;
    }

    ClimatePipelineOptions options;
    if (arguments[0] == "benchmark")
        options.benchmark = true;
    else if (arguments[0] != "run")
    {
        cerr << "Unknown command: " << arguments[0] << '\n';
        printusage();
        return 1;
    }

    ensureworkingdirectory();
    reloadappenvironment();
    if (!parse(vector<string>(arguments.begin() + 1, arguments.end()), options))
    {
        printusage();
        return 1;
    }

    const int columns = options.width;
    const int rows = options.width / 2;
    const AppEnvironmentConfig& environment = getappenvironment();
    const filesystem::path landpath = options.importLandPath.empty()
        ? resolutionvariant(environment.earthBenchmarkLandPath, columns, rows) : filesystem::path(options.importLandPath);
    const filesystem::path seapath = options.importSeaPath.empty()
        ? resolutionvariant(environment.earthBenchmarkSeaPath, columns, rows) : filesystem::path(options.importSeaPath);

    climateworld world;
    world.resize(columns - 1, rows - 1);
    auto terrain = std::make_unique<ImportedTerrain>(columns, rows);
    climateinputs inputs = earthinputs(options, *terrain);
    const auto setnom = [&terrain](int x, int y, int value) { terrain->nom[x][y] = value; };
    if (!importlandheightmapcells(landpath.string(), columns - 1, rows - 1, inputs.sealevel, nullptr, setnom))
    {
        cerr << "Failed to import land map: " << landpath.string() << '\n';
        return 1;
    }
    if (!importseadepthmapcells(seapath.string(), columns - 1, rows - 1, inputs.sealevel, nullptr, setnom))
    {
        cerr << "Failed to import sea map: " << seapath.string() << '\n';
        return 1;
    }
    world.bindinputs(inputs);
    countclimatelandandsea(world);
    applyclimateoptions(options.climate, world.climatesimulation());

    WorldGenerationDebugOptions debugoptions;
    debugoptions.logToProfilingWorkbook = false;
    updatereport("Climate on imported terrain, seed " + to_string(options.seed) + ":");
    updatereport("");
    beginworldgendebugrun(options.seed, &debugoptions);
    begintimedreporting();

    vector<runbundle::StepRecord> steps;
    ClimateStageContext context;
    for (const GenerationStageDefinition& stage : getgenerationstages())
    {
        if (!isclimateonlystage(stage.id))
            continue;
        const auto start = chrono::steady_clock::now();
        beginworldgenstep(stage.label);
        executeclimatestage(stage.id, world, context);
        runbundle::StepRecord record;
        record.index = static_cast<int>(steps.size());
        record.process = getworldprocess(stage.id).id;
        record.label = stage.label;
        record.seconds = chrono::duration<double>(chrono::steady_clock::now() - start).count();
        for (worldmodel::FieldId write : getworldprocess(stage.id).writes)
            record.writes.push_back(worldmodel::qualifiedname(write));
        steps.push_back(record);
    }

    exportclimatevalidationreport(world);
    endtimedreporting();
    endworldgendebugrun();

    if (options.benchmark)
    {
        int runid = 0;
        if (!recordclimatebenchmarkrun(world, options.benchmarkInformation, options.maps, &runid))
        {
            cerr << "Failed to record climate benchmark artifacts.\n";
            return 1;
        }
        updatereport("Climate benchmark run ID: " + to_string(runid));
        updatereport("Climate benchmark log: " + environment.climateBenchmarkRunLogPath.string());
    }
    else if (options.mapsSpecified && !exportproceduralclimateartifacts(world, options.maps))
    {
        cerr << "Failed to write climate maps.\n";
        return 1;
    }

    if (!options.bundlePath.empty())
    {
        try
        {
            writebundle(world, options, steps);
        }
        catch (const exception& error)
        {
            cerr << "Failed to write run bundle: " << error.what() << '\n';
            return 1;
        }
        updatereport("Run bundle: " + options.bundlePath);
    }
    updatereport("Climate completed.");
    return 0;
}
