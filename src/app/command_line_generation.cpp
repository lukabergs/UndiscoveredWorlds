// Windows headers first: later project headers bring std names into scope.
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include "command_line_generation.hpp"
#include "climate_options.hpp"

#include <algorithm>
#include <functional>
#include <cmath>
#include <iostream>
#include <memory>

#include "app_environment.hpp"
#include "classes.hpp"
#include "climate_artifacts.hpp"
#include "functions.hpp"
#include "world_stages.hpp"
#include "map_imports.hpp"
#include "physical_layers.hpp"
#include "planet_bundle.hpp"
#include "region.hpp"
#include "world_recipes.hpp"
#include "world_generation_debug.hpp"

using namespace std;

void printcommandlineusage()
{
    cout << "world_pipeline <command> [options]\n";
    cout << "Commands:\n";
    cout << "  run                 Run a recipe (--recipe) or part of it; nothing is registered.\n";
    cout << "  generate            Generate a procedural world, or climate on imported terrain with --climate-only.\n";
    cout << "  climate-benchmark   Run the earth-climate recipe on the Earth terrain and record a numbered run.\n";
    cout << "Options of run:\n";
    cout << "  --recipe <name>         procedural-world, imported-world, earth-climate or tectonics. Imported recipes use the\n";
    cout << "                           Earth terrain of configs/app.env at --resolution (default 512) unless\n";
    cout << "                           --import-land and --import-sea or --from-bundle are given.\n";
    cout << "  --from-bundle <dir>     Start from the world of a run bundle or checkpoint. Its seed and grid apply;\n";
    cout << "                           other options are not inherited.\n";
    cout << "  --steps <a..b>          Run steps a to b of the recipe (indices or process ids; a.. and ..b allowed).\n";
    cout << "  --only <process>        Run one step.\n";
    cout << "  --checkpoint-at <step|Myr[,...]>  Write a bundle after these steps (process ids) and at these\n";
    cout << "                           times of the tectonics step (numbers, in Myr).\n";
    cout << "  --checkpoint-every <Myr>  Write a bundle at this interval of the tectonics step.\n";
    cout << "  --checkpoint-dir <dir>  Checkpoint directory (default <bundle>-checkpoints); one bundle per checkpoint.\n";
    cout << "  --duration <Myr>        Tectonic duration (as --tectonic-duration-myr).\n";
    cout << "  --mineral-history <on|off>  Record mineral deposits during tectonics (default on).\n";
    cout << "                           The tectonics recipe runs the engine on the world grid (--resolution, default\n";
    cout << "                           from configs/tectonics-baseline.json) with the default planet body.\n";
    cout << "Options:\n";
    cout << "  --climate-only          With generate: recompute climate on imported terrain without modifying it.\n";
    cout << "  --seed <number>         Use a fixed world seed (default 1).\n";
    cout << "  --save <path>           Save the generated world as a .uww file.\n";
    cout << "  --bundle <directory>    Write every world field and table as a run bundle (new or empty directory).\n";
    cout << "  --recipe <name>         procedural-world, imported-world or earth-climate (default from the command).\n";
    cout << "  --audit-access          Hash every field around each step and report undeclared writes.\n";
    cout << "  --access-report <file>  Write measured per-step access as JSON (reads need the tracking build).\n";
    cout << "  --tectonic-duration-myr <Myr>  Set continuous tectonic duration in Myr (0 uses recipe).\n";
    cout << "  --plates <count>        Set tectonic plate count.\n";
    cout << "  --world-width <pixels>  Set generated world width in pixels.\n";
    cout << "  --world-height <pixels> Set generated world height in pixels.\n";
    printclimateoptionsusage(cout);
    printtectonicoptionsusage(cout);
    cout << "  --resolution <width>    Set an even global width and derive height as width/2.\n";
    cout << "  --reference-precip <csv> Compare against a precipitation reference grid.\n";
    cout << "  --import-land <tiff>    Import a land height map instead of generating terrain.\n";
    cout << "  --import-sea <tiff>     Import a sea depth map instead of generating terrain.\n";
    cout << "  --benchmark-info <text> Store a change summary in the benchmark run log.\n";
    cout << "  --map <id[,id...]>      Replace default benchmark maps; repeat to add maps.\n";
    cout << "                           Default: koppen,temp,precip,precip-tif; Jan/Jul surface/upper speed, LIC, particles; SST, ocean currents.\n";
    cout << "                           Seasons: jan, apr, jul, oct. Use --map all or --map none.\n";
    cout << "  --print-climate-relative-error Print per-climate relative error against the Earth benchmark counts.\n";
    cout << "  --no-rivers             Disable river generation.\n";
    cout << "  --no-lakes              Disable lake generation.\n";
    cout << "  --no-deltas             Disable delta generation.\n";
    cout << "  --social                Enable social generation.\n";
    cout << "  --social-mode <mode>    static or historical.\n";
    cout << "  --no-prehistory         Disable prehistory pass in historical mode.\n";
    cout << "  --history-years <years> Historical simulation span.\n";
    cout << "  --no-profile-log        Skip profiling workbook logging.\n";
    cout << "  --help                  Show this message.\n";
}

bool parsecommandlineoptions(const vector<string>& arguments, CommandLineGenerationOptions& options)
{
    const int argc = static_cast<int>(arguments.size());
    const auto cleanup = []() {};
    options.commandArguments = arguments;

    for (int index = 0; index < argc; index++)
    {
        const string& argument = arguments[index];

        if (argument == "--help" || argument == "-h" || argument == "/?")
        {
            options.showHelp = true;
            continue;
        }

        if (argument == "--print-climate-relative-error")
        {
            options.printClimateRelativeError = true;
            continue;
        }

        if (argument == "--climate-only")
        {
            options.climateOnly = true;
            continue;
        }

        if (argument == "--audit-access")
        {
            options.auditAccess = true;
            continue;
        }
        if (argument == "--no-rivers")
        {
            options.rivers = false;
            continue;
        }

        if (argument == "--no-lakes")
        {
            options.lakes = false;
            continue;
        }

        if (argument == "--no-deltas")
        {
            options.deltas = false;
            continue;
        }

        if (argument == "--social")
        {
            options.socialEnabled = true;
            continue;
        }

        if (argument == "--no-prehistory")
        {
            options.usePrehistory = false;
            continue;
        }

        if (argument == "--no-profile-log")
        {
            options.logToProfilingWorkbook = false;
            continue;
        }

        if (isclimateoption(argument) ||
            argument == "--map-resolution" || argument == "--wind-map-resolution" || argument == "--seed" || argument == "--save" || argument == "--bundle" || argument == "--recipe" || argument == "--from-bundle" || argument == "--steps" || argument == "--only" || argument == "--checkpoint-at" || argument == "--checkpoint-every" || argument == "--checkpoint-dir" || argument == "--duration" || argument == "--mineral-history" || argument == "--access-report" || argument == "--tectonic-duration-myr" || argument == "--plates" || argument == "--world-width" || argument == "--world-height" || argument == "--resolution" || argument == "--map" || argument == "--maps" || argument == "--reference-precip" || argument == "--import-land" || argument == "--import-sea" || argument == "--benchmark-info" || argument == "--social-mode" || argument == "--history-years")
        {
            if (index + 1 >= argc)
            {
                cerr << "Missing value for " << argument << '\n';
                cleanup();
                return false;
            }

            const string value = arguments[++index];

            try
            {
                string climatereason;
                const ClimateOptionResult climateresult = parseclimateoption(argument, value, options.climate, &climatereason);
                if (climateresult == ClimateOptionResult::invalid)
                    throw invalid_argument(climatereason);
                if (climateresult == ClimateOptionResult::parsed)
                {
                }
                else if (argument == "--seed")
                {
                    options.seed = stol(value);
                    options.hasSeed = true;
                }
                else if (argument == "--bundle")
                {
                    options.bundlePath = value;
                }
                else if (argument == "--recipe")
                {
                    WorldRecipe recipe;
                    if (!findworldrecipe(value, recipe))
                        throw invalid_argument("recipe");
                    options.recipeName = value;
                }
                else if (argument == "--from-bundle")
                {
                    options.fromBundlePath = value;
                }
                else if (argument == "--steps" || argument == "--only")
                {
                    if (!options.stepsText.empty() || value.empty())
                        throw invalid_argument("steps");
                    options.stepsText = argument == "--only" ? value + ".." + value : value;
                }
                else if (argument == "--checkpoint-at")
                {
                    size_t start = 0;
                    while (start <= value.size())
                    {
                        const size_t comma = min(value.find(',', start), value.size());
                        if (comma > start)
                        {
                            // A number is a time of the tectonics step, anything else a step.
                            const string item = value.substr(start, comma - start);
                            size_t parsed = 0;
                            double timemyr = 0.0;
                            try { timemyr = stod(item, &parsed); } catch (const exception&) { parsed = 0; }
                            if (parsed == item.size())
                                options.checkpointTimesMyr.push_back(timemyr);
                            else
                                options.checkpointAt.push_back(item);
                        }
                        start = comma + 1;
                    }
                }
                else if (argument == "--checkpoint-every")
                {
                    options.checkpointEveryMyr = stod(value);
                    if (!(options.checkpointEveryMyr > 0.0))
                        throw invalid_argument("checkpoint-every");
                }
                else if (argument == "--duration")
                {
                    options.plateTectonicsDurationMyr = stod(value);
                    if (!(options.plateTectonicsDurationMyr > 0.0))
                        throw invalid_argument("duration");
                }
                else if (argument == "--mineral-history")
                {
                    if (value != "on" && value != "off")
                        throw invalid_argument("mineral-history");
                    options.mineralHistory = value == "on";
                }
                else if (argument == "--checkpoint-dir")
                {
                    options.checkpointDirectory = value;
                }
                else if (argument == "--access-report")
                {
                    options.accessReportPath = value;
                }
                else if (argument == "--tectonic-duration-myr")
                {
                    options.plateTectonicsDurationMyr = max(0.0, stod(value));
                }
                else if (argument == "--plates")
                {
                    options.plateTectonicsPlateCount = max(1, stoi(value));
                }
                else if (argument == "--world-width")
                {
                    options.worldWidth = stoi(value);
                    options.hasWorldWidth = true;
                }
                else if (argument == "--world-height")
                {
                    options.worldHeight = stoi(value);
                    options.hasWorldHeight = true;
                }
                else if (argument == "--resolution")
                {
                    options.resolution = stoi(value);
                    if (options.resolution < 2 || options.resolution % 2 != 0)
                        throw invalid_argument("resolution");
                    options.hasResolution = true;
                }
                else if (argument == "--map-resolution" || argument == "--wind-map-resolution")
                {
                    // Retired: every map now has the world grid's size. Accepted
                    // so archived experiment commands still replay.
                    cerr << "Ignoring " << argument << ": exported maps use the world grid (--resolution).\n";
                }
                else if (argument == "--map" || argument == "--maps")
                {
                    if (!options.benchmarkMapsSpecified)
                    {
                        options.benchmarkMaps.clear();
                        options.benchmarkMapsSpecified = true;
                    }
                    string failure;
                    if (!addclimatebenchmarkmapargument(
                            options.benchmarkMaps, value, &failure))
                    {
                        cerr << failure << '\n';
                        cleanup();
                        return false;
                    }
                }
                else if (argument == "--reference-precip")
                {
                    options.referencePath = value;
                    options.hasReferencePath = true;
                }
                else if (argument == "--import-land")
                {
                    options.importLandPath = value;
                    options.hasImportLandPath = true;
                }
                else if (argument == "--import-sea")
                {
                    options.importSeaPath = value;
                    options.hasImportSeaPath = true;
                }
                else if (argument == "--benchmark-info")
                {
                    options.benchmarkInformation = value;
                }
                else if (argument == "--social-mode")
                {
                    if (value == "static")
                        options.socialMode = SocialGenerationOptions::Mode::static_ex_nihilo;
                    else if (value == "historical")
                        options.socialMode = SocialGenerationOptions::Mode::historical;
                    else
                    {
                        cerr << "Invalid value for --social-mode: " << value << '\n';
                        cleanup();
                        return false;
                    }
                }
                else if (argument == "--history-years")
                {
                    options.historyYears = std::max(10, stoi(value));
                }
                else
                {
                    options.savePath = value;
                    options.hasSavePath = true;
                }
            }
            catch (const exception&)
            {
                cerr << "Invalid value for " << argument << ": " << value << '\n';
                cleanup();
                return false;
            }

            continue;
        }

        // Every other --key value pair is a tectonic scenario option, or unknown.
        if (index + 1 < argc)
        {
            try
            {
                if (parsetectonicoption(options.tectonicScenario, argument, arguments[index + 1]))
                {
                    options.tectonicScenarioGiven = true;
                    index++;
                    continue;
                }
            }
            catch (const exception& error)
            {
                cerr << error.what() << '\n';
                cleanup();
                return false;
            }
        }

        cerr << "Unknown argument: " << argument << '\n';
        cleanup();
        return false;
    }

    cleanup();
    return true;
}

filesystem::path commandlinevalidationdirectory(long seed)
{
    filesystem::path outputroot = getappenvironment().climateValidationDirectory;

    if (outputroot.empty())
        outputroot = filesystem::current_path();

    return outputroot / ("seed_" + to_string(seed));
}

SocialGenerationOptions makesocialoptions(bool enabled, SocialGenerationOptions::Mode mode, bool useprehistory, int historyyears)
{
    SocialGenerationOptions options;
    options.enabled = enabled;
    options.mode = mode;
    options.usePrehistory = mode == SocialGenerationOptions::Mode::historical ? useprehistory : false;
    options.historyYears = std::max(10, historyyears);
    return options;
}

namespace
{
// Step notes, access violations and, on request, the per-step access report.
bool reportrecipe(const WorldRecipe& recipe, const RecipeRunResult& result, const string& accessreportpath)
{
    for (const RecipeStepRecord& step : result.steps)
        for (const string& violation : step.violations)
            updatereport(step.label + ": " + violation);
    if (result.violations > 0)
        updatereport("Recipe " + recipe.name + ": " + to_string(result.violations) + " access violation(s)");

    if (accessreportpath.empty())
        return true;
    ofstream out(accessreportpath, ios::binary);
    out << "{\n  \"schema\": \"uw-access-report/v1\",\n  \"recipe\": " << runbundle::jsonstring(recipe.name)
        << ",\n  \"tracking_reads\": " << (worldmodel::access::tracking ? "true" : "false") << ",\n  \"steps\": [";
    auto list = [&](const auto& ids)
    {
        string text = "[";
        for (size_t index = 0; index < ids.size(); index++)
            text += (index ? ", " : "") + runbundle::jsonstring(worldmodel::qualifiedname(ids[index]));
        return text + "]";
    };
    for (size_t index = 0; index < result.steps.size(); index++)
    {
        const RecipeStepRecord& step = result.steps[index];
        out << (index ? ",\n" : "\n") << "    {\"index\": " << step.index << ", \"process\": " << runbundle::jsonstring(step.process)
            << ", \"label\": " << runbundle::jsonstring(step.label) << ", \"seconds\": " << step.seconds
            << ", \"measured\": " << (step.measured ? "true" : "false") << ", \"drew_random\": " << (step.drewrandom ? "true" : "false")
            << ", \"writes\": " << list(step.writes) << ", \"table_writes\": " << list(step.tablewrites)
            << ", \"reads\": " << list(step.reads) << ", \"table_reads\": " << list(step.tablereads);
        // Hashes after the step: two runs diverge at the first step whose hashes differ.
        out << ", \"write_hashes\": {";
        for (size_t write = 0; write < step.writehashes.size(); write++)
        {
            char hex[17];
            snprintf(hex, sizeof(hex), "%016llx", static_cast<unsigned long long>(step.writehashes[write]));
            out << (write ? ", " : "") << runbundle::jsonstring(worldmodel::qualifiedname(step.writes[write])) << ": \"" << hex << '"';
        }
        out << "}, \"violations\": [";
        for (size_t violation = 0; violation < step.violations.size(); violation++)
            out << (violation ? ", " : "") << runbundle::jsonstring(step.violations[violation]);
        out << "]}";
    }
    out << "\n  ]\n}\n";
    return static_cast<bool>(out);
}
}

namespace
{
// Runs a recipe on terrain already in the world (generated, imported or loaded
// from a bundle). A run through the last step of a climate recipe writes the
// climate validation report.
RecipeRunResult runrecipeonterrain(planet& world, const WorldRecipe& recipe, GenerationScratch& scratch,
    WorldStageParameters& parameters, const GenerationExecutionContext& context, const RecipeRunOptions& runoptions)
{
    // Geosphere is frozen in climate recipes; keep the explicit elevation check.
    const bool frozenterrain = find(recipe.frozen.begin(), recipe.frozen.end(), worldmodel::Component::geosphere) != recipe.frozen.end();
    vector<int> inputterrain;
    if (frozenterrain) for (int y=0;y<=world.height();++y) for (int x=0;x<=world.width();++x)
        inputterrain.push_back(world.nom(x,y));

    RecipeRunResult result = runrecipe(recipe, world, scratch, parameters, context, runoptions);
    if (!result.completed)
        return result;

    if (frozenterrain) {
        if (!runoptions.selectSteps || runoptions.steps.last + 1 == recipe.steps.size())
            exportclimatevalidationreport(world.climate());
        size_t cell=0;
        for (int y=0;y<=world.height();++y) for (int x=0;x<=world.width();++x)
            if (world.nom(x,y)!=inputterrain[cell++]) throw runtime_error("Climate-only stage modified terrain");
        updatereport("Climate-only input terrain preserved exactly after import.");
    }
    return result;
}
}

RecipeRunResult runimportedworldrecipe(planet& world, const WorldRecipe& recipe, const GenerationExecutionContext& context,
    const RecipeRunOptions& runoptions)
{
    world.setmaxelevation(200000);
    // Land and sea counts describe the imported terrain; climate stages read them.
    getlandandseatotals(world);

    GenerationScratch scratch;
    initializegenerationscratch(scratch);
    makecontinentalshelves(world, world.continentalshelves(), 4);

    WorldStageParameters parameters;
    parameters.tectonicDurationMyr = tectonicdurationmyr();
    parameters.seaLevel = world.sealevel();

    return runrecipeonterrain(world, recipe, scratch, parameters, context, runoptions);
}

void completeimportedworldgeneration(planet& world, bool dorivers, bool dolakes, bool dodeltas, bool appendclimateworkbook, const SocialGenerationOptions& socialoptions, boolshapetemplate smalllake[], boolshapetemplate largelake[], boolshapetemplate landshape[], vector<vector<bool>>& okmountains, const ImportedClimateMaps* importedclimate, bool climateOnly)
{
    GenerationExecutionContext context;
    context.dorivers = dorivers;
    context.climateOnly = climateOnly;
    context.dolakes = dolakes;
    context.dodeltas = dodeltas;
    context.appendclimateworkbook = appendclimateworkbook;
    context.socialoptions = socialoptions;
    context.landshape = landshape;
    context.smalllake = smalllake;
    context.largelake = largelake;
    context.okmountains = &okmountains;
    context.importedClimate = importedclimate;

    const WorldRecipe recipe = climateOnly ? earthclimaterecipe() : importedworldrecipe();
    const RecipeRunResult result = runimportedworldrecipe(world, recipe, context, {});
    if (!result.completed)
    {
        if (climateOnly)
            throw runtime_error("Climate-only generation failed: " + result.error);
        updatereport(result.error);
    }
    reportrecipe(recipe, result, "");
}

int runcommandlineworldgeneration(const CommandLineGenerationOptions& options)
{
    auto world = std::make_unique<planet>();
    initialiseworld(*world);
    initialisemapcolours(*world);
    world->clear();

    const bool frombundle = !options.fromBundlePath.empty();
    const bool hasimport = options.hasImportLandPath || options.hasImportSeaPath;
    WorldRecipe recipe;
    if (!options.recipeName.empty())
        findworldrecipe(options.recipeName, recipe);
    if (options.runCommand)
    {
        if (options.recipeName.empty() || options.climateOnly)
        {
            cerr << "run requires --recipe (earth-climate replaces --climate-only).\n";
            return 1;
        }
        if (hasimport && (frombundle || !recipeimportsterrain(recipe)))
        {
            if (frombundle)
                cerr << "--from-bundle replaces imported terrain.\n";
            else
                cerr << "Recipe " << recipe.name << " cannot start from imported terrain.\n";
            return 1;
        }
    }
    else if (frombundle || !options.stepsText.empty() || !options.checkpointAt.empty() ||
        !options.checkpointTimesMyr.empty() || options.checkpointEveryMyr > 0.0)
    {
        cerr << "--from-bundle, --steps, --only and checkpoints belong to the run command.\n";
        return 1;
    }

    // run: a recipe on imported terrain imports the Earth terrain unless told otherwise.
    const bool importedworld = options.runCommand ? recipeimportsterrain(recipe) && !frombundle
        : options.earthClimateBenchmark || hasimport;
    // The tectonics module runs alone: the engine's grid is the world grid.
    const bool tectonicsmodule = recipe.name == "tectonics";
    if (options.climateOnly && !importedworld) {
        cerr << "--climate-only requires imported terrain.\n";
        return 1;
    }
    const long seed = options.hasSeed ? options.seed : 1;
    world->setseed(seed);

    fast_srand(world->seed());

    if (options.hasReferencePath)
        setreferenceprecipitationgridpath(options.referencePath);

    if (options.earthClimateBenchmark == false && (options.hasImportLandPath != options.hasImportSeaPath))
    {
        cerr << "Imported runs require both --import-land and --import-sea for generate.\n";
        return 1;
    }

    if (importedworld && (options.hasWorldWidth != options.hasWorldHeight))
    {
        cerr << "Imported runs require both --world-width and --world-height when either is specified.\n";
        return 1;
    }

    if (options.hasResolution && (options.hasWorldWidth || options.hasWorldHeight))
    {
        cerr << "Use either --resolution or --world-width/--world-height, not both.\n";
        return 1;
    }

    if (options.hasWorldWidth &&
        (options.worldWidth < 2 || options.worldWidth % 2 != 0 ||
         options.worldHeight != options.worldWidth / 2))
    {
        cerr << "World dimensions must be an even W by W/2.\n";
        return 1;
    }

    if (options.benchmarkMapsSpecified && !options.earthClimateBenchmark && !options.climateOnly && !options.runCommand)
    {
        cerr << "--map requires run, climate-benchmark or generate --climate-only.\n";
        return 1;
    }

    boolshapetemplate landshape[12];
    createlandshapetemplates(landshape);
    boolshapetemplate chainland[2];
    createchainlandtemplates(chainland);
    boolshapetemplate smalllake[12];
    createsmalllaketemplates(smalllake);
    boolshapetemplate largelake[10];
    createlargelaketemplates(largelake);

    vector<int> squareroot((MAXCRATERRADIUS * MAXCRATERRADIUS + MAXCRATERRADIUS + 1) * 24);

    for (int n = 1; n < squareroot.size(); n++)
        squareroot[n] = (int)sqrt(n);

    WorldGenerationDebugOptions debugoptions;
    debugoptions.visualizePlateTectonicsRealtime = false;
    debugoptions.logToProfilingWorkbook = options.logToProfilingWorkbook;
    debugoptions.plateTectonicsDurationMyr = options.plateTectonicsDurationMyr;
    debugoptions.plateTectonicsPlateCount = options.plateTectonicsPlateCount;

    if (options.earthClimateBenchmark || (options.runCommand && importedworld))
    {
        world->setwidth(std::clamp(
            options.hasResolution ? options.resolution :
                (options.hasWorldWidth ? options.worldWidth : kEarthBenchmarkWidth),
            2,
            ARRAYWIDTH) - 1);
    }
    else if (importedworld)
    {
        if (options.hasResolution)
        {
            world->setwidth(std::clamp(options.resolution, 2, ARRAYWIDTH) - 1);
        }
        else if (options.hasWorldWidth && options.hasWorldHeight)
        {
            world->setwidth(std::clamp(options.worldWidth, 2, ARRAYWIDTH) - 1);
        }
    }
    else if (!frombundle)
    {
        // A world draws its planet body from the seed; the tectonics module keeps the default body.
        if (!tectonicsmodule)
            changeworldproperties(*world);
        const int defaultwidth = tectonicsmodule ? static_cast<int>(options.tectonicScenario.width) : world->width() + 1;
        world->setwidth(std::clamp(options.hasResolution ? options.resolution :
            (options.hasWorldWidth ? options.worldWidth : defaultwidth), 2, ARRAYWIDTH) - 1);
    }

    PlanetBundleSource source;
    if (frombundle)
    {
        try
        {
            source = loadplanetbundle(*world, options.fromBundlePath);
        }
        catch (const exception& error)
        {
            cerr << "Failed to load run bundle: " << error.what() << '\n';
            return 1;
        }
        if (options.hasWorldWidth || (options.hasResolution && options.resolution != world->width() + 1))
        {
            cerr << "--from-bundle fixes the grid at " << world->width() + 1 << " columns.\n";
            return 1;
        }
        if (options.hasSeed)
            world->setseed(options.seed);
        fast_srand(world->seed());
    }

    // The fractal noise of the terrain passes does not terminate on smaller grids.
    if (options.runCommand && world->width() + 1 < 32)
    {
        cerr << "The world grid needs at least 32 columns.\n";
        return 1;
    }

    applyclimateoptions(options.climate, world->climatesimulation());
    const SocialGenerationOptions socialoptions = makesocialoptions(options.socialEnabled, options.socialMode, options.usePrehistory, options.historyYears);

    updatereport((importedworld ? "Generating imported world from seed: " : "Generating world from seed: ") + to_string(world->seed()) + ":");
    updatereport("");

    if (options.recipeName.empty())
        recipe = importedworld ? (options.climateOnly || options.earthClimateBenchmark ? earthclimaterecipe() : importedworldrecipe()) : proceduralworldrecipe();
    if (!frombundle && importedworld != recipeimportsterrain(recipe))
    {
        cerr << "Recipe " << recipe.name << (importedworld ? " cannot start from imported terrain.\n" : " requires imported terrain.\n");
        return 1;
    }
    RecipeRunOptions runoptions;
    runoptions.auditAccess = options.auditAccess;
    if (!options.stepsText.empty())
    {
        string reason;
        if (!parserecipesteps(recipe, options.stepsText, runoptions.steps, &reason))
        {
            cerr << "Invalid step selection: " << reason << '\n';
            return 1;
        }
        runoptions.selectSteps = true;
        if (runoptions.steps.first > 0 && !frombundle)
        {
            cerr << "A run that starts after the first step needs --from-bundle.\n";
            return 1;
        }
    }
    if (frombundle)
    {
        runoptions.fromBundle = true;
        runoptions.loadedfields = source.loadedfields;
        runoptions.loadedtables = source.loadedtables;
        runoptions.sourceprovenance = source.provenance;
    }

    // Every bundle of the run, final or checkpoint, is written the same way.
    const auto writebundle = [&](const filesystem::path& directory, const RecipeRunResult& result)
    {
        runbundle::Metadata metadata;
        metadata.seed = world->seed();
        metadata.columns = world->width() + 1;
        metadata.rows = world->height() + 1;
        metadata.recipe = recipe.name;
        metadata.command = options.commandArguments;
        metadata.steps = bundlesteps(result);
        metadata.properties.push_back({ "executable_built", string(__DATE__) + " " + __TIME__ });
        metadata.properties.push_back({ "command", options.runCommand ? "run" : options.earthClimateBenchmark ? "climate-benchmark" : "generate" });
        metadata.properties.push_back({ "access", options.auditAccess ? "measured" : "declared" });
        if (frombundle)
            metadata.properties.push_back({ "from_bundle", options.fromBundlePath });

        // What this run wrote replaces what the source bundle carried.
        vector<worldmodel::FieldId> writtenfields;
        vector<worldmodel::TableId> writtentables;
        for (const RecipeStepRecord& step : result.steps)
        {
            const WorldProcessDefinition& process = getworldprocess(step.stage);
            writtenfields.insert(writtenfields.end(), process.writes.begin(), process.writes.end());
            writtenfields.insert(writtenfields.end(), step.writes.begin(), step.writes.end());
            writtentables.insert(writtentables.end(), process.tablewrites.begin(), process.tablewrites.end());
            writtentables.insert(writtentables.end(), step.tablewrites.begin(), step.tablewrites.end());
        }
        writeplanetbundle(*world, directory, metadata, result.provenance, frombundle ? &source : nullptr, writtenfields, writtentables);
    };

    vector<size_t> checkpointsteps;
    for (const string& step : options.checkpointAt)
    {
        size_t index = 0;
        if (!findrecipestep(recipe, step, index))
        {
            cerr << "Recipe " << recipe.name << " has no step " << step << " to checkpoint.\n";
            return 1;
        }
        checkpointsteps.push_back(index);
    }
    // Checkpoints inside an evolver see the run as of the last finished step.
    RecipeRunResult finishedsteps;
    std::function<void(GenerationStageId, double)> evolvercheckpoint;
    if (!checkpointsteps.empty() || !options.checkpointTimesMyr.empty() || options.checkpointEveryMyr > 0.0)
    {
        if (options.checkpointDirectory.empty() && options.bundlePath.empty())
        {
            cerr << "Checkpoints need --checkpoint-dir or --bundle.\n";
            return 1;
        }
        const filesystem::path checkpointdirectory = options.checkpointDirectory.empty()
            ? filesystem::path(options.bundlePath + "-checkpoints") : filesystem::path(options.checkpointDirectory);
        evolvercheckpoint = [&, checkpointdirectory](GenerationStageId stage, double timemyr)
        {
            // The running step is listed with its declared writes: they are in the world as of this time.
            const WorldProcessDefinition& process = getworldprocess(stage);
            RecipeRunResult partial = finishedsteps;
            RecipeStepRecord record;
            record.index = static_cast<int>(find(recipe.steps.begin(), recipe.steps.end(), stage) - recipe.steps.begin());
            record.stage = stage;
            record.process = process.id;
            record.label = getgenerationstage(stage).label;
            record.writes = process.writes;
            record.tablewrites = process.tablewrites;
            for (worldmodel::FieldId write : record.writes)
                partial.provenance.fields[write] = record.label;
            for (worldmodel::TableId write : record.tablewrites)
                partial.provenance.tables[write] = record.label;
            partial.steps.push_back(record);

            // One bundle per time: <recipe index>-<process id>-t<Myr>.
            char time[32];
            snprintf(time, sizeof(time), "-t%08.3f", timemyr);
            const filesystem::path directory = checkpointdirectory /
                ((record.index < 10 ? "0" : "") + to_string(record.index) + "-" + record.process + time);
            writebundle(directory, partial);
            updatereport("Checkpoint: " + directory.string());
        };
        runoptions.afterstep = [&, checkpointdirectory](const RecipeRunResult& sofar, const RecipeStepRecord& step, string& error)
        {
            finishedsteps = sofar;
            if (find(checkpointsteps.begin(), checkpointsteps.end(), static_cast<size_t>(step.index)) == checkpointsteps.end())
                return true;
            // One bundle per step: <recipe index>-<process id>.
            const filesystem::path directory = checkpointdirectory / ((step.index < 10 ? "0" : "") + to_string(step.index) + "-" + step.process);
            try
            {
                writebundle(directory, sofar);
            }
            catch (const exception& failure)
            {
                error = string("checkpoint failed: ") + failure.what();
                return false;
            }
            updatereport("Checkpoint: " + directory.string());
            return true;
        };
    }
    RecipeRunResult runresult;

    beginworldgendebugrun(world->seed(), &debugoptions);
    begintimedreporting();

    if (importedworld)
    {
        const AppEnvironmentConfig& appenv = getappenvironment();
        const filesystem::path landpath = options.hasImportLandPath
            ? filesystem::path(options.importLandPath)
            : resolutionvariant(
                appenv.earthBenchmarkLandPath, world->width() + 1, world->height() + 1);
        const filesystem::path seapath = options.hasImportSeaPath
            ? filesystem::path(options.importSeaPath)
            : resolutionvariant(
                appenv.earthBenchmarkSeaPath, world->width() + 1, world->height() + 1);
        vector<vector<bool>> okmountains(ARRAYWIDTH, vector<bool>(ARRAYHEIGHT, false));

        if (importlandheightmap(*world, landpath.string()) == false)
        {
            endtimedreporting();
            endworldgendebugrun();
            cerr << "Failed to import land map: " << landpath.string() << '\n';
            return 1;
        }

        if (importseadepthmap(*world, seapath.string()) == false)
        {
            endtimedreporting();
            endworldgendebugrun();
            cerr << "Failed to import sea map: " << seapath.string() << '\n';
            return 1;
        }

        GenerationExecutionContext context;
        context.dorivers = options.rivers;
        context.climateOnly = recipe.name == "earth-climate";
        context.dolakes = options.lakes;
        context.dodeltas = options.deltas;
        context.socialoptions = socialoptions;
        context.landshape = landshape;
        context.smalllake = smalllake;
        context.largelake = largelake;
        context.okmountains = &okmountains;

        runresult = runimportedworldrecipe(*world, recipe, context, runoptions);
    }
    else
    {
        GenerationScratch scratch;
        initializegenerationscratch(scratch);

        WorldStageParameters parameters;
        parameters.tectonicDurationMyr = options.plateTectonicsDurationMyr;
        parameters.tectonicPlateCount = options.plateTectonicsPlateCount;
        parameters.seaLevel = world->sealevel();
        if (options.tectonicScenarioGiven)
            parameters.tectonicScenario = std::make_shared<const tectonics::scenario::Scenario>(options.tectonicScenario);
        parameters.tectonicWorldGrid = tectonicsmodule;
        parameters.tectonicMineralHistory = options.mineralHistory;
        parameters.tectonicCheckpointEveryMyr = options.checkpointEveryMyr;
        parameters.tectonicCheckpointsMyr = options.checkpointTimesMyr;

        GenerationExecutionContext context;
        context.evolverCheckpoint = evolvercheckpoint;
        context.dorivers = options.rivers;
        context.dolakes = options.lakes;
        context.dodeltas = options.deltas;
        context.appendclimateworkbook = false;
        context.socialoptions = socialoptions;
        context.landshape = landshape;
        context.chainland = chainland;
        context.smalllake = smalllake;
        context.largelake = largelake;
        context.climateOnly = recipe.name == "earth-climate";

        runresult = runrecipeonterrain(*world, recipe, scratch, parameters, context, runoptions);
    }

    endtimedreporting();
    endworldgendebugrun();

    if (!runresult.completed)
    {
        cerr << "Recipe " << recipe.name << " failed: " << runresult.error << '\n';
        if (recipe.name == "earth-climate" || options.runCommand)
            return 1;
    }
    if (!reportrecipe(recipe, runresult, options.accessReportPath))
    {
        cerr << "Failed to write access report: " << options.accessReportPath << '\n';
        return 1;
    }
    // An audited run with a violation still writes its outputs and then fails.
    const bool auditfailed = options.auditAccess && runresult.violations > 0;
    if (auditfailed)
        cerr <<"Access audit found " << runresult.violations << " undeclared or frozen write(s).\n";

    if (options.earthClimateBenchmark || options.climateOnly)
    {
        int benchmarkrunid = options.climateOnly ? 1 : 0;

        if ((options.climateOnly ? exportproceduralclimateartifacts(world->climate(), options.benchmarkMaps) : recordclimatebenchmarkrun(
                world->climate(),
                options.benchmarkInformation,
                options.benchmarkMaps,
                &benchmarkrunid)) == false)
        {
            cerr << "Failed to record climate benchmark artifacts.\n";
            return 1;
        }

        const AppEnvironmentConfig& appenv = getappenvironment();
        updatereport("Climate benchmark run ID: " + to_string(benchmarkrunid));
        if (!options.climateOnly) updatereport("Climate benchmark log: " + appenv.climateBenchmarkRunLogPath.string());
        for (const climatebenchmarkmaprequest& request : options.benchmarkMaps.requests())
        {
            if (!climatebenchmarkmapisreference(request.kind))
            {
                updatereport("Climate benchmark map: " +
                    ((request.kind == climatebenchmarkmapkind::precipitationtiff
                        ? appenv.climateBenchmarkImageDirectory.parent_path() / "fields"
                        : appenv.climateBenchmarkImageDirectory) /
                        climatebenchmarkmaprelativepath(request, benchmarkrunid)).string());
            }
        }
        updatereport("Climate benchmark map guide: " +
            (appenv.climateBenchmarkImageDirectory / "benchmark_map_guide.txt").string());
    }
    else if (options.runCommand && options.benchmarkMapsSpecified &&
        !exportproceduralclimateartifacts(world->climate(), options.benchmarkMaps))
    {
        cerr << "Failed to write climate maps.\n";
        return 1;
    }

    if (options.printClimateRelativeError)
        printclimaterelativeerrorreport(world->climate());

    updatereport("");
    updatereport("World generation completed.");
    updatereport("Validation output: " + commandlinevalidationdirectory(world->seed()).string());

    if (options.hasSavePath)
    {
        world->saveworld(options.savePath);
        updatereport("Saved world: " + options.savePath);
    }

    if (!options.bundlePath.empty())
    {
        try
        {
            writebundle(options.bundlePath, runresult);
        }
        catch (const exception& error)
        {
            cerr << "Failed to write run bundle: " << error.what() << '\n';
            return 1;
        }
        updatereport("Run bundle: " + options.bundlePath);
    }

    return auditfailed ? 1 : 0;
}

