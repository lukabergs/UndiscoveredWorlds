// Windows headers first: later headers bring std names into scope.
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include "core/app_environment.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <unordered_map>

using namespace std;

namespace
{
AppEnvironmentConfig environmentconfig;
bool environmentloaded = false;

string trimcopy(const string& value)
{
    size_t start = 0;
    size_t end = value.size();

    while (start < end && isspace(static_cast<unsigned char>(value[start])))
        start++;
    while (end > start && isspace(static_cast<unsigned char>(value[end - 1])))
        end--;

    return value.substr(start, end - start);
}

string getenvoverride(const char* name, const unordered_map<string, string>& filevalues)
{
    const char* envvalue = getenv(name);

    if (envvalue != nullptr && envvalue[0] != '\0')
        return envvalue;

    const auto found = filevalues.find(name);

    if (found != filevalues.end())
        return found->second;

    return "";
}

void applypathoverride(filesystem::path& target, const string& value)
{
    if (!value.empty())
        target = value;
}

void loadappenvironmentinternal()
{
    environmentconfig = AppEnvironmentConfig{};

    unordered_map<string, string> filevalues;
    ifstream infile("configs/app.env");
    string line;

    while (getline(infile, line))
    {
        line = trimcopy(line);

        if (line.empty() || line[0] == '#')
            continue;

        const size_t equals = line.find('=');

        if (equals == string::npos)
            continue;

        string key = trimcopy(line.substr(0, equals));
        string value = trimcopy(line.substr(equals + 1));

        if (!value.empty() && value.front() == '"' && value.back() == '"' && value.size() >= 2)
            value = value.substr(1, value.size() - 2);

        filevalues[key] = value;
    }

    applypathoverride(environmentconfig.defaultWorldDirectory, getenvoverride("UW_DEFAULT_WORLD_DIR", filevalues));
    applypathoverride(environmentconfig.defaultAppearanceDirectory, getenvoverride("UW_DEFAULT_UWS_DIR", filevalues));
    applypathoverride(environmentconfig.defaultImageDirectory, getenvoverride("UW_DEFAULT_IMAGE_DIR", filevalues));
    applypathoverride(environmentconfig.terrainTemplatePath, getenvoverride("UW_TERRAIN_TEMPLATES", filevalues));
    applypathoverride(environmentconfig.climateDiagnosticsDirectory, getenvoverride("UW_CLIMATE_DIAGNOSTICS_DIR", filevalues));
    applypathoverride(environmentconfig.earthBenchmarkLandPath, getenvoverride("UW_EARTH_BENCHMARK_LAND", filevalues));
    applypathoverride(environmentconfig.earthBenchmarkSeaPath, getenvoverride("UW_EARTH_BENCHMARK_SEA", filevalues));

    if (!environmentconfig.defaultAppearanceDirectory.empty() && environmentconfig.defaultAppearanceDirectory != ".")
        filesystem::create_directories(environmentconfig.defaultAppearanceDirectory);
    if (!environmentconfig.defaultImageDirectory.empty() && environmentconfig.defaultImageDirectory != ".")
        filesystem::create_directories(environmentconfig.defaultImageDirectory);
    if (!environmentconfig.climateDiagnosticsDirectory.empty())
        filesystem::create_directories(environmentconfig.climateDiagnosticsDirectory);

    environmentloaded = true;
}
}

const AppEnvironmentConfig& getappenvironment()
{
    if (environmentloaded == false)
        loadappenvironmentinternal();

    return environmentconfig;
}

void reloadappenvironment()
{
    loadappenvironmentinternal();
}

filesystem::path resolutionvariant(
    const filesystem::path& source,
    int columns,
    int rows)
{
    const filesystem::path candidate = source.parent_path() /
        (source.stem().string() + "_" + to_string(columns) + "x" +
            to_string(rows) + source.extension().string());
    return filesystem::exists(candidate) ? candidate : source;
}

void ensureworkingdirectory()
{
    if (GetFileAttributesA("version.txt") != INVALID_FILE_ATTRIBUTES)
        return;

    char exepath[MAX_PATH] = { 0 };
    GetModuleFileNameA(NULL, exepath, MAX_PATH);

    filesystem::path directory = filesystem::path(exepath).parent_path();
    while (!directory.empty())
    {
        if (filesystem::is_regular_file(directory / "version.txt"))
        {
            filesystem::current_path(directory);
            return;
        }
        const filesystem::path parent = directory.parent_path();
        if (parent == directory)
            break;
        directory = parent;
    }
}
