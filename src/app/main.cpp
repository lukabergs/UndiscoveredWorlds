
#include <iostream>

#include <windows.h>
#include <shellapi.h>

#include "app/application.hpp"

int main()
{
    ensureworkingdirectory();
    reloadappenvironment();

    int argc = 0;
    if (LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc))
        LocalFree(argv);
    if (argc > 1)
    {
        cerr << "UndiscoveredWorlds is the desktop application and takes no arguments.\n"
                "Use world_pipeline.exe for headless generation and climate benchmarks.\n";
        return 1;
    }

    Application application;
    application.run();

    return 0;
}
