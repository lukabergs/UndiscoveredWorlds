// Windows headers first: later project headers bring std names into scope.
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include "command_line_generation.hpp"
#include "climate_options.hpp"

#include <algorithm>
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
    cout << "  generate            Generate a procedural world, or climate on imported terrain with --climate-only.\n";
    cout << "  climate-benchmark   Run climate on the imported Earth benchmark terrain and record a numbered run.\n";
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
            argument == "--map-resolution" || argument == "--wind-map-resolution" || argument == "--seed" || argument == "--save" || argument == "--bundle" || argument == "--recipe" || argument == "--access-report" || argument == "--tectonic-duration-myr" || argument == "--plates" || argument == "--world-width" || argument == "--world-height" || argument == "--resolution" || argument == "--map" || argument == "--maps" || argument == "--reference-precip" || argument == "--import-land" || argument == "--import-sea" || argument == "--benchmark-info" || argument == "--social-mode" || argument == "--history-years")
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

    // Geosphere is frozen in climate recipes; keep the explicit elevation check.
    const bool frozenterrain = find(recipe.frozen.begin(), recipe.frozen.end(), worldmodel::Component::geosphere) != recipe.frozen.end();
    vector<int> inputterrain;
    if (frozenterrain) for (int y=0;y<=world.height();++y) for (int x=0;x<=world.width();++x)
        inputterrain.push_back(world.nom(x,y));

    RecipeRunResult result = runrecipe(recipe, world, scratch, parameters, context, runoptions);
    if (!result.completed)
        return result;

    if (frozenterrain) {
        exportclimatevalidationreport(world.climate());
        size_t cell=0;
        for (int y=0;y<=world.height();++y) for (int x=0;x<=world.width();++x)
            if (world.nom(x,y)!=inputterrain[cell++]) throw runtime_error("Climate-only stage modified terrain");
        updatereport("Climate-only input terrain preserved exactly after import.");
    }
    return result;
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

    const bool importedworld = options.earthClimateBenchmark || options.hasImportLandPath || options.hasImportSeaPath;
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

    if (options.benchmarkMapsSpecified && !options.earthClimateBenchmark && !options.climateOnly)
    {
        cerr << "--map requires climate-benchmark or generate --climate-only.\n";
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

    if (options.earthClimateBenchmark)
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
    else if (importedworld == false)
    {
        changeworldproperties(*world);
        world->setwidth(std::clamp(options.hasResolution ? options.resolution :
            (options.hasWorldWidth ? options.worldWidth : world->width() + 1), 2, ARRAYWIDTH) - 1);
    }

    applyclimateoptions(options.climate, world->climatesimulation());
    const SocialGenerationOptions socialoptions = makesocialoptions(options.socialEnabled, options.socialMode, options.usePrehistory, options.historyYears);

    updatereport((importedworld ? "Generating imported world from seed: " : "Generating world from seed: ") + to_string(world->seed()) + ":");
    updatereport("");

    WorldRecipe recipe;
    if (!options.recipeName.empty())
        findworldrecipe(options.recipeName, recipe);
    else
        recipe = importedworld ? (options.climateOnly ? earthclimaterecipe() : importedworldrecipe()) : proceduralworldrecipe();
    if (importedworld == (recipe.name == "procedural-world"))
    {
        cerr << "Recipe " << recipe.name << (importedworld ? " cannot start from imported terrain.\n" : " requires imported terrain.\n");
        return 1;
    }
    RecipeRunOptions runoptions;
    runoptions.auditAccess = options.auditAccess;
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

        GenerationExecutionContext context;
        context.dorivers = options.rivers;
        context.dolakes = options.lakes;
        context.dodeltas = options.deltas;
        context.appendclimateworkbook = false;
        context.socialoptions = socialoptions;
        context.landshape = landshape;
        context.chainland = chainland;
        context.smalllake = smalllake;
        context.largelake = largelake;

        runresult = runrecipe(recipe, *world, scratch, parameters, context, runoptions);
    }

    endtimedreporting();
    endworldgendebugrun();

    if (!runresult.completed)
    {
        cerr << "Recipe " << recipe.name << " failed: " << runresult.error << '\n';
        if (recipe.name == "earth-climate")
            return 1;
    }
    if (!reportrecipe(recipe, runresult, options.accessReportPath))
    {
        cerr << "Failed to write access report: " << options.accessReportPath << '\n';
        return 1;
    }
    if (options.auditAccess && runresult.violations > 0)
        cerr << "Access audit found " << runresult.violations << " undeclared or frozen write(s).\n";

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
        runbundle::Metadata metadata;
        metadata.seed = world->seed();
        metadata.columns = world->width() + 1;
        metadata.rows = world->height() + 1;
        metadata.recipe = recipe.name;
        metadata.command = options.commandArguments;
        metadata.steps = bundlesteps(runresult);
        metadata.properties.push_back({ "executable_built", string(__DATE__) + " " + __TIME__ });
        metadata.properties.push_back({ "command", options.earthClimateBenchmark ? "climate-benchmark" : "generate" });
        metadata.properties.push_back({ "access", options.auditAccess ? "measured" : "declared" });
        try
        {
            writeplanetbundle(*world, options.bundlePath, metadata, runresult.provenance);
        }
        catch (const exception& error)
        {
            cerr << "Failed to write run bundle: " << error.what() << '\n';
            return 1;
        }
        updatereport("Run bundle: " + options.bundlePath);
    }

    return 0;
}

