#include "app_environment.hpp"

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
    ifstream infile("app.env");
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
    applypathoverride(environmentconfig.profilingWorkbookPath, getenvoverride("UW_PROFILING_WORKBOOK", filevalues));

    if (!environmentconfig.defaultAppearanceDirectory.empty() && environmentconfig.defaultAppearanceDirectory != ".")
        filesystem::create_directories(environmentconfig.defaultAppearanceDirectory);
    if (!environmentconfig.defaultImageDirectory.empty() && environmentconfig.defaultImageDirectory != ".")
        filesystem::create_directories(environmentconfig.defaultImageDirectory);
    if (!environmentconfig.profilingWorkbookPath.empty() && environmentconfig.profilingWorkbookPath.has_parent_path())
        filesystem::create_directories(environmentconfig.profilingWorkbookPath.parent_path());

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
