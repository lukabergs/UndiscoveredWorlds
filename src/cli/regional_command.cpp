#include "cli/regional_command.hpp"

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

#include "io/regional_bundle.hpp"
#include "pipeline/recipe_parameters.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
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
    out << "  --scale <k>             Fine cells per global cell edge, a whole number >= 1 (required).\n";
    out << "  --bundle <dir>          The regional bundle (new or empty directory; required).\n";
    out << "  --set <name=value>, --parameters <file>  As for run (landscape.regional.* sets the pass,\n";
    out << "                           halo_cells the halo); the manifest records every value. The stream power\n";
    out << "                           law (landscape.stream_power.erodibility, area_exponent) is first set to\n";
    out << "                           the world bundle's recorded values; --set and --parameters override it.\n";
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
    for (const char* name : {"landscape.stream_power.erodibility", "landscape.stream_power.area_exponent"})
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
            argument != "--scale" && argument != "--seed" && argument != "--set" && argument != "--parameters")
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
    const char* missing = options.fromBundlePath.empty() ? "--from-bundle"
                        : options.bundlePath.empty()     ? "--bundle"
                        : options.scale == 0             ? "--scale"
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
    const int halo = tuning::landscape::regional::haloCells;
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

    const landscape::RegionalGlobalInputs inputs = gatherregionalinputs(*world);
    const landscape::RegionalParameters parameters = regionalparameters(inputs.seed);
    landscape::RegionalFields fields;
    const auto started = chrono::steady_clock::now();
    try
    {
        fields = landscape::runregionalpass(inputs, window, parameters);
    }
    catch (const exception& error)
    {
        cerr << "Regional pass failed: " << error.what() << '\n';
        return 1;
    }
    const double seconds = chrono::duration<double>(chrono::steady_clock::now() - started).count();

    regionalbundle::Metadata metadata;
    metadata.run.seed = world->seed();
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
    metadata.sourceseed = world->seed();
    metadata.sealevel = world->sealevel();
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
    cout << "Regional pass: " << seconds << " s; block-mean error " << fields.blockmeanerror_m << " m, relief "
         << (fields.reliefconverged ? "converged" : "not converged") << ", " << fields.changedreceivers
         << " receivers changed, lowering mean " << fields.meanlowering_m << " m, largest " << fields.maxlowering_m << " m\n";
    cout << "Regional bundle: " << options.bundlePath << '\n';
    return 0;
}
