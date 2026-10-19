#include "cli/regional_command.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "io/regional_bundle.hpp"
#include "pipeline/recipe_parameters.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/landscape/refinement_grid.hpp"
#include "simulations/landscape/refinement_levels.hpp"
#include "simulations/landscape/regional_cube.hpp"
#include "simulations/landscape/regional_pass.hpp"
#include "simulations/regional/regional_inputs.hpp"
#include "world/planet.hpp"
#include "world/planet_bundle.hpp"
#include "world/world_setup.hpp"

using namespace std;

void printregionusage(ostream& out)
{
    out << "Options of region (a window of a world bundle at a scale, written as a regional bundle):\n";
    out << "  --from-bundle <dir>     The world bundle; its seed and grid apply (required).\n";
    out << "  --window <w>,<s>,<e>,<n>  The window in degrees, snapped outwards to whole global cells;\n";
    out << "                           west > east crosses the date line.\n";
    out << "  --pixels <x0>,<y0>,<x1>,<y1>  The window in global cells (pixels of the global map): x0,y0 the\n";
    out << "                           top-left cell (column from the west edge at -180 degrees, row from the north),\n";
    out << "                           x1,y1 the bottom-right cell, both inclusive; x0 > x1 crosses the date line.\n";
    out << "                           Exactly one of --window and --pixels.\n";
    out << "  --scale <k>             Fine cells per global cell edge, a whole number >= 1. A power of\n";
    out << "                           landscape.refinement.factor (4, 16, 64, ... at factor 4; at most\n";
    out << "                           refinement.stored_scale_max) runs the refinement hierarchy: log_factor k\n";
    out << "                           levels, each projected exactly onto its parent; any other k one pass.\n";
    out << "  --metres-per-pixel <m>  Instead of --scale (exactly one of the two): the pixel's north-south size on\n";
    out << "                           the ground. The hierarchy runs to the coarsest level whose cells are at most\n";
    out << "                           m and is box-averaged to m. Aligned (no resample) when m is a whole number n\n";
    out << "                           of that level's cells and n divides its k; otherwise resampled, partial\n";
    out << "                           pixels at the east and south edges dropped. The window block records level,\n";
    out << "                           factor, metres_per_pixel, resampled and the pixel grid. Refused at or above\n";
    out << "                           the global cell size and below the cells of refinement.stored_scale_max.\n";
    out << "  --bundle <dir>          The regional bundle (new or empty directory; required).\n";
    out << "  --set <name=value>, --parameters <file>  As for run (landscape.regional.* sets the pass,\n";
    out << "                           halo_cells the halo, for the hierarchy at least 2 x mean_passes + 7 +\n";
    out << "                           2 x landscape.refinement.smooth_passes); the manifest records every\n";
    out << "                           value. The stream power\n";
    out << "                           law (landscape.stream_power.erodibility, area_exponent) and the glacial\n";
    out << "                           cooling (landscape.glacial.cooling_c) are first set to the world bundle's\n";
    out << "                           recorded values; --set and --parameters override them.\n";
    out << "  --seed <number>         Must equal the world bundle's seed when given.\n";
}

namespace
{
struct RegionOptions
{
    bool showHelp = false;
    string fromBundlePath;
    string bundlePath;
    bool hasWindow = false;
    bool hasPixels = false;
    double window[4] = {};  // west, south, east, north
    int pixels[4] = {};     // x0, y0, x1, y1
    int scale = 0;
    double metresPerPixel = 0.0;  // 0: not given
    bool hasSeed = false;
    long seed = 0;
    vector<string> commandArguments;
};

bool parseinteger(const string& text, long long& value)
{
    if (text.empty())
        return false;
    size_t parsed = 0;
    try { value = stoll(text, &parsed); } catch (const exception&) { return false; }
    return parsed == text.size();
}

bool parsereal(const string& text, double& value)
{
    if (text.empty())
        return false;
    size_t parsed = 0;
    try { value = stod(text, &parsed); } catch (const exception&) { return false; }
    return parsed == text.size() && isfinite(value);
}

// Four comma-separated items.
bool splitfour(const string& text, string items[4])
{
    size_t start = 0;
    for (int index = 0; index < 4; index++)
    {
        const size_t comma = text.find(',', start);
        if ((comma == string::npos) != (index == 3))
            return false;
        items[index] = text.substr(start, comma == string::npos ? string::npos : comma - start);
        start = comma + 1;
    }
    return true;
}

// The stream power law of the solve is the world's: the erodibility and the area exponent
// the world bundle's manifest records in its "parameters" (a flat object of scalar values
// keyed by full name, recipeparametersjson) are set before --set and --parameters, which
// may still override them. Returns false on a value the parameter cannot take; a bundle
// without a readable manifest or without the values sets nothing (its load fails later
// or the defaults stay). `applied` names the values taken from the world.
bool applyworldstreampower(const vector<string>& arguments, vector<string>& applied)
{
    string bundle;
    for (size_t index = 0; index + 1 < arguments.size(); index++)
        if (arguments[index] == "--from-bundle")
            bundle = arguments[index + 1];
    if (bundle.empty())
        return true;
    ifstream file(filesystem::path(bundle) / "manifest.json", ios::binary);
    if (!file)
        return true;
    ostringstream content;
    content << file.rdbuf();
    const string text = content.str();
    const size_t key = text.find("\"parameters\"");
    const size_t open = key == string::npos ? string::npos : text.find('{', key);
    const size_t close = open == string::npos ? string::npos : text.find('}', open);
    if (close == string::npos)
        return true;
    // Phase E5 (E5.0): the glacial cooling too, so the present climate's ELA offset of the
    // regional pass (regional_inputs.cpp) is the run's.
    for (const char* name : {"landscape.stream_power.erodibility", "landscape.stream_power.area_exponent",
                             "landscape.glacial.cooling_c"})
    {
        const size_t found = text.find("\"" + string(name) + "\"", open);
        if (found == string::npos || found > close)
            continue;
        size_t start = text.find(':', found);
        if (start == string::npos || start > close)
            continue;
        start = text.find_first_not_of(" \t\r\n", start + 1);
        if (start == string::npos || start >= close)
            continue;
        const size_t end = text.find_first_of(",}\r\n", start);
        string value = text.substr(start, end - start);
        value.erase(value.find_last_not_of(" \t") + 1);
        string reason;
        if (!setrecipeparameter(string(name) + "=" + value, &reason))
        {
            cerr << "region: the world bundle's " << name << " = " << value << ": " << reason << '\n';
            return false;
        }
        applied.push_back(name);
    }
    return true;
}

bool parseregionoptions(const vector<string>& arguments, RegionOptions& options)
{
    options.commandArguments = arguments;
    for (size_t index = 0; index < arguments.size(); index++)
    {
        const string& argument = arguments[index];
        if (argument == "--help" || argument == "-h" || argument == "/?")
        {
            options.showHelp = true;
            continue;
        }
        if (argument == "--resolution" || argument == "--world-width" || argument == "--world-height")
        {
            cerr << "region: " << argument << " is not accepted; the grid is the world bundle's.\n";
            return false;
        }
        if (argument != "--from-bundle" && argument != "--bundle" && argument != "--window" && argument != "--pixels" &&
            argument != "--scale" && argument != "--metres-per-pixel" && argument != "--seed" && argument != "--set" &&
            argument != "--parameters")
        {
            cerr << "region: unknown option " << argument << '\n';
            return false;
        }
        if (index + 1 >= arguments.size())
        {
            cerr << "region: missing value for " << argument << '\n';
            return false;
        }
        const string& value = arguments[++index];
        if (argument == "--from-bundle")
            options.fromBundlePath = value;
        else if (argument == "--bundle")
            options.bundlePath = value;
        else if (argument == "--window" || argument == "--pixels")
        {
            if (options.hasWindow || options.hasPixels)
            {
                cerr << "region: give exactly one of --window and --pixels, once.\n";
                return false;
            }
            string items[4];
            bool valid = splitfour(value, items);
            for (int item = 0; valid && item < 4; item++)
            {
                if (argument == "--window")
                    valid = parsereal(items[item], options.window[item]);
                else
                {
                    long long cell = 0;
                    valid = parseinteger(items[item], cell) && cell >= 0 && cell <= 1000000;
                    options.pixels[item] = static_cast<int>(cell);
                }
            }
            if (!valid)
            {
                cerr << "region: " << argument << " takes four comma-separated "
                     << (argument == "--window" ? "numbers: west,south,east,north (degrees)"
                                                : "whole numbers >= 0: x0,y0,x1,y1 (global cells)")
                     << ", not " << value << '\n';
                return false;
            }
            (argument == "--window" ? options.hasWindow : options.hasPixels) = true;
        }
        else if (argument == "--scale")
        {
            long long scale = 0;
            if (!parseinteger(value, scale) || scale < 1 || scale > 4096)
            {
                cerr << "region: --scale takes a whole number from 1 to 4096, not " << value << '\n';
                return false;
            }
            options.scale = static_cast<int>(scale);
        }
        else if (argument == "--metres-per-pixel")
        {
            double metres = 0.0;
            if (!parsereal(value, metres) || !(metres > 0.0))
            {
                cerr << "region: --metres-per-pixel takes a positive number of metres, not " << value << '\n';
                return false;
            }
            options.metresPerPixel = metres;
        }
        else if (argument == "--seed")
        {
            long long seed = 0;
            if (!parseinteger(value, seed) || seed < numeric_limits<long>::min() || seed > numeric_limits<long>::max())
            {
                cerr << "region: --seed takes a whole number from " << numeric_limits<long>::min() << " to "
                     << numeric_limits<long>::max() << ", not " << value << '\n';
                return false;
            }
            options.seed = static_cast<long>(seed);
            options.hasSeed = true;
        }
        else
        {
            string reason;
            const bool accepted = argument == "--set" ? setrecipeparameter(value, &reason) : loadrecipeparameters(value, &reason);
            if (!accepted)
            {
                cerr << argument << ": " << reason << '\n';
                return false;
            }
        }
    }
    if (options.showHelp)
        return true;
    if (options.scale != 0 && options.metresPerPixel > 0.0)
    {
        cerr << "region: give one of --scale and --metres-per-pixel, not both.\n";
        return false;
    }
    const char* missing = options.fromBundlePath.empty() ? "--from-bundle"
                        : options.bundlePath.empty()     ? "--bundle"
                        : options.scale == 0 && !(options.metresPerPixel > 0.0) ? "--scale or --metres-per-pixel"
                        : !options.hasWindow && !options.hasPixels ? "--window or --pixels"
                                                         : nullptr;
    if (missing != nullptr)
    {
        cerr << "region requires " << missing << ".\n";
        return false;
    }
    return true;
}
}

landscape::RegionalWindow regionalwindowfromcells(int worldwidth, int scale, int halo, int x0, int y0, int x1, int y1)
{
    if (worldwidth < 4 || worldwidth % 2 != 0)
        throw invalid_argument("the world width must be even and at least 4");
    const int worldheight = worldwidth / 2;
    if (x0 < 0 || x1 < 0 || x0 >= worldwidth || x1 >= worldwidth)
        throw invalid_argument("columns must lie within 0.." + to_string(worldwidth - 1));
    if (y0 < 0 || y1 >= worldheight || y0 > y1)
        throw invalid_argument("rows must satisfy 0 <= y0 <= y1 <= " + to_string(worldheight - 1));
    const int columns = x0 <= x1 ? x1 - x0 + 1 : worldwidth - x0 + x1 + 1;
    if (columns + 2 * halo > worldwidth)
        throw invalid_argument("the window (" + to_string(columns) + " columns) and its halo of " + to_string(halo) +
                               " cells on each side are wider than the planet (" + to_string(worldwidth) + " columns)");
    // Every column from x0 > 0 round to x0 - 1: its west and east edges are one meridian,
    // which regionalwindow reads as a box of width 0; the whole circle starts at column 0.
    if (columns == worldwidth && x0 != 0)
        throw invalid_argument("a window of all " + to_string(worldwidth) + " columns is given as x0 = 0, x1 = " +
                               to_string(worldwidth - 1));
    // The cells' edges: column i spans -180 + 360 i / W .. -180 + 360 (i + 1) / W, row j
    // spans 90 - 180 j / H .. 90 - 180 (j + 1) / H; x0 > x1 gives west > east.
    const double west = -180.0 + 360.0 * x0 / worldwidth;
    const double east = -180.0 + 360.0 * (x1 + 1) / worldwidth;
    const double north = 90.0 - 180.0 * y0 / worldheight;
    const double south = 90.0 - 180.0 * (y1 + 1) / worldheight;
    const landscape::RegionalWindow window = landscape::regionalwindow(worldwidth, scale, halo, west, south, east, north);
    if (window.column0 != x0 || window.columns != columns || window.row0 != y0 || window.rows != y1 - y0 + 1)
        throw logic_error("regionalwindowfromcells: the edges snapped to other cells");
    return window;
}

namespace
{
// The manifest's run record and source block of a regional bundle.
regionalbundle::Metadata regionmetadata(const RegionOptions& options, const planet& world, const vector<string>& fromworld,
                                        const string& manifesthash, double seconds)
{
    regionalbundle::Metadata metadata;
    metadata.run.seed = world.seed();
    metadata.run.recipe = "region";
    metadata.run.command = options.commandArguments;
    runbundle::StepRecord step;
    step.index = 0;
    step.process = "landscape.regional_pass";
    step.label = "region";
    step.seconds = seconds;
    metadata.run.steps.push_back(step);
    metadata.run.parameters = recipeparametersjson();
    metadata.run.properties.push_back({ "executable_built", string(__DATE__) + " " + __TIME__ });
    metadata.run.properties.push_back({ "command", "region" });
    metadata.run.properties.push_back({ "from_bundle", options.fromBundlePath });
    // The values set from the world bundle before --set and --parameters ("parameters" holds
    // the values used).
    string fromworldnames;
    for (const string& name : fromworld)
        fromworldnames += (fromworldnames.empty() ? "" : ",") + name;
    metadata.run.properties.push_back({ "parameters_from_world", fromworldnames });
    metadata.sourcebundle = options.fromBundlePath;
    metadata.sourcemanifestsha256 = manifesthash;
    metadata.sourceseed = world.seed();
    metadata.sealevel = world.sealevel();
    return metadata;
}

// A cube world's window (T8 W2-W2f; sphere-mesh.md section 8): the box in the lat-lon export's pixels
// (4 N columns) selects the level-0 cells whose centres lie in it, the halo counts rings of level-0
// cells, the fields are the kept level cells, written as a cube-window bundle with the export of the box.
int runcuberegion(const RegionOptions& options, const planet& world, const landscape::RegionalGlobalInputs& inputs,
                  const landscape::RefinementParameters& refinement, const vector<string>& fromworld,
                  const string& manifesthash)
{
    const mesh::CellGrid& base = *inputs.grid;
    const int exportcolumns = 4 * base.facecells();
    if (options.metresPerPixel > 0.0)
    {
        cerr << "region: --metres-per-pixel is not ported to cube worlds; give --scale.\n";
        return 1;
    }
    if (tuning::landscape::regional::tributaryEnabled != 0)
    {
        cerr << "region: the tributaries (landscape.regional.tributary_enabled 1) are not ported to cube worlds.\n";
        return 1;
    }
    const int levels = landscape::levelofscale(options.scale, refinement.factor);
    const bool hierarchy = levels >= 1;
    const int halo = hierarchy ? std::max(tuning::landscape::regional::haloCells,
                                          landscape::refinementhalo(tuning::landscape::regional::meanPasses,
                                                                    refinement.smoothpasses, false))
                               : tuning::landscape::regional::haloCells;
    if (hierarchy && options.scale > refinement.storedscalemax)
    {
        cerr << "region: --scale " << options.scale << " is above landscape.refinement.stored_scale_max ("
             << refinement.storedscalemax << "), the finest scale written as a bundle.\n";
        return 1;
    }
    landscape::RegionalWindow window;
    vector<mesh::LevelKey> kept0;
    try
    {
        // The export box (no halo: on the cube the halo is rings of level-0 cells, not columns).
        window = options.hasPixels
            ? regionalwindowfromcells(exportcolumns, options.scale, 0, options.pixels[0], options.pixels[1], options.pixels[2],
                                      options.pixels[3])
            : landscape::regionalwindow(exportcolumns, options.scale, 0, options.window[0], options.window[1],
                                        options.window[2], options.window[3]);
        kept0 = landscape::cellsinpixelbox(base, exportcolumns, window.column0,
                                           window.row0, (window.column0 + window.columns - 1) % exportcolumns,
                                           window.row0 + window.rows - 1);
    }
    catch (const exception& error)
    {
        cerr << "region: " << (options.hasPixels ? "--pixels" : "--window") << ": " << error.what() << '\n';
        return 1;
    }
    if (kept0.empty())
    {
        cerr << "region: the box holds no cube cell centre; widen it.\n";
        return 1;
    }
    window.halo = halo;

    const landscape::RegionalParameters parameters = regionalparameters(inputs.seed);
    landscape::CubeRegionResult result;
    const auto started = chrono::steady_clock::now();
    try
    {
        result = hierarchy ? landscape::runcuberefinement(inputs, kept0, halo, levels, parameters, refinement)
                           : landscape::runcuberegionalpass(inputs, kept0, options.scale, halo, parameters);
    }
    catch (const exception& error)
    {
        cerr << "Regional pass failed: " << error.what() << '\n';
        return 1;
    }
    const double seconds = chrono::duration<double>(chrono::steady_clock::now() - started).count();
    result.fields.window = window;
    regionalbundle::Metadata metadata = regionmetadata(options, world, fromworld, manifesthash, seconds);
    if (hierarchy)
    {
        metadata.level = levels;
        metadata.factor = refinement.factor;
        metadata.levels = result.levels;
    }
    try
    {
        regionalbundle::writecubewindowbundle(options.bundlePath, result.fields, *result.kept, base.facecells(),
                                              hierarchy ? levels : 0, metadata);
    }
    catch (const exception& error)
    {
        cerr << "Failed to write regional bundle: " << error.what() << '\n';
        return 1;
    }
    const landscape::RegionalFields& fields = result.fields;
    cout << "Regional window (cube N " << base.facecells() << "): export columns " << window.column0 << " + " << window.columns
         << ", rows " << window.row0 << " + " << window.rows << " of " << exportcolumns << " x " << exportcolumns / 2
         << "; scale " << options.scale << ", halo " << halo << " rings; " << kept0.size() << " level-0 cells, "
         << result.kept->cellcount() << " kept cells, domain " << result.domaincells.back() << " cells\n";
    cout << "Regional pass: " << seconds << " s; block-mean error " << fields.blockmeanerror_m << " m, relief "
         << (fields.reliefconverged ? "converged" : "not converged") << ", " << fields.changedreceivers
         << " receivers changed, lowering mean " << fields.meanlowering_m << " m, largest " << fields.maxlowering_m << " m";
    if (fields.sharepasses >= 0)
        cout << "; land share met after " << fields.sharepasses << " passes";
    else if (fields.sharemisses > 0)
        cout << "; land share not met: " << fields.sharemisses << " blocks missing after every pass";
    cout << '\n';
    for (std::size_t index = 0; index < result.levels.size(); index++)
    {
        const landscape::LevelDiagnostics& level = result.levels[index];
        const landscape::ProjectionReport& projection = level.projection;
        cout << "Level " << level.level << " (" << level.cellsize_km << " km cells, kept " << result.keptcells[index]
             << ", domain " << result.domaincells[index] << " cells): " << level.seconds << " s; block-mean error "
             << level.blockmeanerror_m << " m after the presmoothing; parent-mean error " << projection.maxparentmeanerror_m
             << " m (bound " << level.meanbound_m << " m), global-mean error " << level.globalmeanerror_m
             << " m after the projection; lambda max " << level.smoothing.maxlambdabefore_m << " m before the "
             << level.smoothing.passes << " smoothing passes, " << projection.maxlambda_m << " m after; border step max "
             << projection.maxborderstep_m << " m, p99 " << projection.borderstep99_m << " m; forced " << projection.forcedparents
             << "; dried cells " << level.driedcells << ", dry land below the sea level " << level.landbelowsea
             << ", fold cells " << level.foldcells << ", kept cells draining to the domain ring " << level.ringbasincells << '\n';
    }
    cout << "Regional bundle: " << options.bundlePath << '\n';
    return 0;
}
}

int runregioncommand(const vector<string>& arguments)
{
    RegionOptions options;
    vector<string> fromworld;
    if (!applyworldstreampower(arguments, fromworld))
        return 1;
    if (!parseregionoptions(arguments, options))
    {
        printregionusage(cerr);
        return 1;
    }
    if (options.showHelp)
    {
        printregionusage(cout);
        return 0;
    }
    string parameterproblem;
    if (!validaterecipeparameters(&parameterproblem))
    {
        cerr << "Parameters: " << parameterproblem << '\n';
        return 1;
    }

    const filesystem::path source = options.fromBundlePath;
    auto world = make_unique<planet>();
    initialiseworld(*world);
    initialisemapcolours(*world);
    world->clear();
    string manifesthash;
    try
    {
        manifesthash = regionalbundle::filesha256(source / "manifest.json");
        loadplanetbundle(*world, source);
    }
    catch (const exception& error)
    {
        cerr << "Failed to load world bundle: " << error.what() << '\n';
        return 1;
    }
    if (options.hasSeed && options.seed != world->seed())
    {
        cerr << "region: --seed " << options.seed << " is not the world bundle's seed " << world->seed() << ".\n";
        return 1;
    }

    const int worldwidth = world->width() + 1;
    const landscape::RefinementParameters refinement = refinementparameters();
    const landscape::RegionalGlobalInputs inputs = gatherregionalinputs(*world);
    if (inputs.grid != nullptr && inputs.grid->kind() == mesh::GridKind::cube)
        return runcuberegion(options, *world, inputs, refinement, fromworld, manifesthash);
    const bool resample = options.metresPerPixel > 0.0;
    if (resample)
    {
        // The source level of the resample (multiscale plan S2.4): the coarsest whose cells are
        // at most the pixel, computed and no finer; at most the finest stored level.
        int finest = 0;
        while (landscape::levelscale(refinement.factor, finest + 1) <= refinement.storedscalemax)
            finest++;
        int source = -1;
        try
        {
            source = landscape::resamplesourcelevel(options.metresPerPixel, inputs.radius_m, worldwidth,
                                                    refinement.factor, finest);
        }
        catch (const exception& error)
        {
            cerr << "region: --metres-per-pixel: " << error.what() << '\n';
            return 1;
        }
        const double global_m = 1000.0 * landscape::levelcellsize_km(inputs.radius_m, worldwidth, refinement.factor, 0);
        if (source == 0)
        {
            cerr << "region: --metres-per-pixel " << options.metresPerPixel << " is not finer than the global cells ("
                 << global_m << " m): the world bundle's own grid serves it.\n";
            return 1;
        }
        if (source < 0)
        {
            cerr << "region: --metres-per-pixel " << options.metresPerPixel << " is finer than the finest stored level ("
                 << finest << ", k = " << landscape::levelscale(refinement.factor, finest) << ", cells of "
                 << global_m / static_cast<double>(landscape::levelscale(refinement.factor, finest))
                 << " m; landscape.refinement.stored_scale_max " << refinement.storedscalemax << ").\n";
            return 1;
        }
        options.scale = static_cast<int>(landscape::levelscale(refinement.factor, source));
    }
    // A power of the refinement factor runs the hierarchy (multiscale plan S2.3); any other
    // scale the single pass. The hierarchy's halo also covers the smoothing before each
    // projection (S2.7), the lake flood (S2.8), the forms and the channel-cell mean (S3.8):
    // refinementhalo, 2 x mean_passes + 7 + 2 x smooth_passes global cells, 31 at the defaults;
    // the single pass keeps regional.halo_cells.
    const int levels = landscape::levelofscale(options.scale, refinement.factor);
    const bool hierarchy = levels >= 1;
    // With the tributaries on, both rules widen by landscape::tributaryreach (multiscale plan S3.3).
    const bool tributaries = tuning::landscape::regional::tributaryEnabled != 0;
    const int halo = hierarchy ? std::max(tuning::landscape::regional::haloCells,
                                          landscape::refinementhalo(tuning::landscape::regional::meanPasses,
                                                                    refinement.smoothpasses, tributaries))
                               : (tributaries ? std::max(tuning::landscape::regional::haloCells,
                                                         landscape::regionalhalo(tuning::landscape::regional::meanPasses, true))
                                              : tuning::landscape::regional::haloCells);
    landscape::RegionalWindow window;
    try
    {
        window = options.hasPixels
            ? regionalwindowfromcells(worldwidth, options.scale, halo, options.pixels[0], options.pixels[1],
                                      options.pixels[2], options.pixels[3])
            : landscape::regionalwindow(worldwidth, options.scale, halo, options.window[0], options.window[1],
                                        options.window[2], options.window[3]);
    }
    catch (const exception& error)
    {
        cerr << "region: " << (options.hasPixels ? "--pixels" : "--window") << ": " << error.what() << '\n';
        return 1;
    }

    if (hierarchy && options.scale > refinement.storedscalemax)
    {
        cerr << "region: --scale " << options.scale << " is above landscape.refinement.stored_scale_max ("
             << refinement.storedscalemax << "), the finest scale written as a bundle.\n";
        return 1;
    }

    const landscape::RegionalParameters parameters = regionalparameters(inputs.seed);
    landscape::RegionalFields fields;
    vector<landscape::LevelDiagnostics> leveldiagnostics;
    landscape::ResampledRegionalFields resampled;
    const auto started = chrono::steady_clock::now();
    try
    {
        if (hierarchy)
        {
            landscape::LevelTarget target;
            target.window = window;
            target.level = levels;
            landscape::RefinementResult refined = landscape::runrefinement(inputs, target, parameters, refinement);
            if (resample)
            {
                resampled = landscape::resampleregionalfields(refined, inputs.radius_m, options.metresPerPixel);
                fields = std::move(resampled.fields);
            }
            else
                fields = std::move(refined.fields);
            leveldiagnostics = std::move(refined.levels);
        }
        else
            fields = landscape::runregionalpass(inputs, window, parameters);
    }
    catch (const exception& error)
    {
        cerr << "Regional pass failed: " << error.what() << '\n';
        return 1;
    }
    const double seconds = chrono::duration<double>(chrono::steady_clock::now() - started).count();

    regionalbundle::Metadata metadata = regionmetadata(options, *world, fromworld, manifesthash, seconds);
    if (hierarchy)
    {
        metadata.level = levels;
        metadata.factor = refinement.factor;
        metadata.levels = leveldiagnostics;
    }
    if (resample)
    {
        metadata.metresperpixel = options.metresPerPixel;
        metadata.resampled = resampled.report.resampled();
        metadata.pixelcolumns = resampled.width;
        metadata.pixelrows = resampled.height;
        metadata.pixelcells = resampled.report.pixelcells;
    }
    try
    {
        regionalbundle::writeregionalbundle(options.bundlePath, fields, metadata);
    }
    catch (const exception& error)
    {
        cerr << "Failed to write regional bundle: " << error.what() << '\n';
        return 1;
    }

    cout << "Regional window: columns " << window.column0 << " + " << window.columns << ", rows " << window.row0 << " + "
         << window.rows << " of " << worldwidth << " x " << window.worldheight << "; scale " << window.scale << ", halo "
         << window.halo << "; fine grid " << window.finewidth() << " x " << window.fineheight() << '\n';
    if (resample)
    {
        const landscape::ResampleReport& report = resampled.report;
        cout << "Resample: " << options.metresPerPixel << " m per pixel from level " << report.sourcelevel << " ("
             << report.sourcecell_m << " m cells, " << report.pixelcells << " cells per pixel edge); "
             << resampled.width << " x " << resampled.height << " pixels, "
             << (report.aligned ? "aligned (box means of whole cells, exact)" : "resampled (overlap-weighted, not exact)")
             << '\n';
    }
    cout << "Regional pass: " << seconds << " s; block-mean error " << fields.blockmeanerror_m << " m, relief "
         << (fields.reliefconverged ? "converged" : "not converged") << ", " << fields.changedreceivers
         << " receivers changed, lowering mean " << fields.meanlowering_m << " m, largest " << fields.maxlowering_m << " m";
    if (fields.sharepasses >= 0)
        cout << "; land share met after " << fields.sharepasses << " passes";
    else if (fields.sharemisses > 0)
        cout << "; land share not met: " << fields.sharemisses << " blocks missing after every pass";
    cout << '\n';
    for (const landscape::LevelDiagnostics& level : leveldiagnostics)
    {
        const landscape::ProjectionReport& projection = level.projection;
        cout << "Level " << level.level << " (" << level.cellsize_km << " km cells, domain " << level.window.domaincolumns
             << " x " << level.window.domainrows << "): " << level.seconds << " s; block-mean error "
             << level.blockmeanerror_m << " m after the presmoothing; parent-mean error " << projection.maxparentmeanerror_m
             << " m (bound " << level.meanbound_m << " m), global-mean error " << level.globalmeanerror_m
             << " m after the projection; lambda max " << level.smoothing.maxlambdabefore_m << " m before the "
             << level.smoothing.passes << " smoothing passes, " << projection.maxlambda_m << " m after (step "
             << level.smoothing.maxlambdastepbefore_m << " -> " << projection.maxlambdastep_m << " m), border step max "
             << projection.maxborderstep_m << " m, p99 " << projection.borderstep99_m << " m; beds moved in "
             << projection.bedsmovedparents << " parents, forced " << projection.forcedparents << " (largest "
             << projection.maxforcedmove_m << " m); dried cells " << level.driedcells << ", shores below lakes "
             << level.shoresbelowlake << ", dry land below the sea level " << level.landbelowsea
             << ", kept cells draining to the domain ring " << level.ringbasincells << ", unreached by a catchment "
             << level.unreachedcells << ", sink junctions pinned " << level.sinkjunctions << ", diagonal links "
             << level.diagonallinks << " (through the corner " << level.cornerlinks << "), pieces merged "
             << level.mergedpieces << ", dry tails " << level.drytails << ", dry channel cells in a lake "
             << level.drychannelinlake << ", channel cells between lakes " << level.channelbetweenlakes
             << ", wet channel cells dried above their water " << level.driedchannels
             << ", wet channel cells above their surface " << level.wetchannelsabove << " (largest "
             << level.wetchannelexcess_m << " m)";
        if (level.tributarycells > 0)  // multiscale plan S3.3 (tributaries on)
            cout << "; tributaries: " << level.tributaries << " grown, " << level.tributarycells << " cells, "
                 << level.tributarybedsabove << " beds above the start surface";
        cout << '\n';
    }
    cout << "Regional bundle: " << options.bundlePath << '\n';
    return 0;
}
