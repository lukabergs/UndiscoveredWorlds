// The headless entry point: recipe runs, world generation and the Earth climate
// benchmark. It never opens a window; the desktop application is separate.

#include <windows.h>
#include <shellapi.h>

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include "app_environment.hpp"
#include "command_line_generation.hpp"

using namespace std;

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
        printcommandlineusage();
        return arguments.empty() ? 1 : 0;
    }

    CommandLineGenerationOptions options;
    const string command = arguments[0];
    if (command == "climate-benchmark")
    {
        options.earthClimateBenchmark = true;
        options.logToProfilingWorkbook = false;
    }
    else if (command == "run")
    {
        options.runCommand = true;
        options.logToProfilingWorkbook = false;
    }
    else if (command != "generate")
    {
        cerr << "Unknown command: " << command << '\n';
        printcommandlineusage();
        return 1;
    }
    options.run = true;

    ensureworkingdirectory();
    reloadappenvironment();

    if (!parsecommandlineoptions(vector<string>(arguments.begin() + 1, arguments.end()), options))
    {
        printcommandlineusage();
        return 1;
    }
    if (options.showHelp)
    {
        printcommandlineusage();
        return 0;
    }
    return runcommandlineworldgeneration(options);
}
