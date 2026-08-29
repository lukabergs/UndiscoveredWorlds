#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <windows.h>

#include "app_environment.hpp"
#include "functions.hpp"
#include "generation_tuning.hpp"
#include "world_generation_debug.hpp"

using namespace std;

namespace
{
const WorldGenerationDebugOptions* activeworldgenerationoptions = nullptr;
long activeworldgenerationseed = 0;
vector<pair<string, double>> activeworldgenerationtimings;
function<void()> worldgenerationvisualizationcallback;

bool runprocessandwait(const wstring& commandline)
{
    STARTUPINFOW startupinfo{};
    PROCESS_INFORMATION processinfo{};
    startupinfo.cb = sizeof(startupinfo);
    startupinfo.dwFlags = STARTF_USESHOWWINDOW;
    startupinfo.wShowWindow = SW_HIDE;

    vector<wchar_t> mutablecommand(commandline.begin(), commandline.end());
    mutablecommand.push_back(L'\0');

    if (CreateProcessW(nullptr, mutablecommand.data(), nullptr, nullptr, FALSE, CREATE_NO_WINDOW, nullptr, nullptr, &startupinfo, &processinfo) == FALSE)
        return false;

    WaitForSingleObject(processinfo.hProcess, INFINITE);

    DWORD exitcode = 1;
    GetExitCodeProcess(processinfo.hProcess, &exitcode);

    CloseHandle(processinfo.hThread);
    CloseHandle(processinfo.hProcess);

    return exitcode == 0;
}

void appendprofilingrow(long seed, const vector<pair<string, double>>& timings)
{
    const filesystem::path workbookpath = getappenvironment().profilingWorkbookPath;

    if (filesystem::exists(workbookpath) == false)
        return;

    const filesystem::path temproot = filesystem::temp_directory_path();
    const filesystem::path datapath = temproot / "uw_profiling_log_input.txt";
    const filesystem::path scriptpath = temproot / "uw_profiling_log_excel.ps1";

    {
        ofstream datafile(datapath);

        if (!datafile.is_open())
            return;

        datafile << seed << '\n';

        for (const auto& timing : timings)
        {
            datafile << timing.first << '\t';
            datafile << fixed << setprecision(2) << timing.second << '\n';
        }
    }

    {
        ofstream scriptfile(scriptpath);

        if (!scriptfile.is_open())
            return;

        scriptfile << "param([string]$WorkbookPath, [string]$DataPath)\n";
        scriptfile << "$ErrorActionPreference = 'Stop'\n";
        scriptfile << "$lines = Get-Content -Path $DataPath\n";
        scriptfile << "if ($lines.Count -lt 1) { exit 0 }\n";
        scriptfile << "function Normalize-StepLabel([string]$label) {\n";
        scriptfile << "    if ($null -eq $label) { return '' }\n";
        scriptfile << "    $normalized = $label.ToLowerInvariant().Trim()\n";
        scriptfile << "    $normalized = [System.Text.RegularExpressions.Regex]::Replace($normalized, '\\s+', ' ')\n";
        scriptfile << "    if ($normalized -eq 'calculating ocean rainfal') { return 'calculating ocean rainfall' }\n";
        scriptfile << "    return $normalized\n";
        scriptfile << "}\n";
        scriptfile << "$excel = New-Object -ComObject Excel.Application\n";
        scriptfile << "$excel.Visible = $false\n";
        scriptfile << "$excel.DisplayAlerts = $false\n";
        scriptfile << "$workbook = $excel.Workbooks.Open($WorkbookPath)\n";
        scriptfile << "$sheet = $workbook.Worksheets.Item('LOG')\n";
        scriptfile << "$row = 2\n";
        scriptfile << "while ($true) {\n";
        scriptfile << "    $value = $sheet.Cells.Item($row, 1).Value2\n";
        scriptfile << "    if ($null -eq $value -or [string]::IsNullOrWhiteSpace([string]$value)) { break }\n";
        scriptfile << "    $row++\n";
        scriptfile << "}\n";
        scriptfile << "$culture = [System.Globalization.CultureInfo]::InvariantCulture\n";
        scriptfile << "$sheet.Cells.Item($row, 1).Value2 = [double]::Parse($lines[0], $culture)\n";
        scriptfile << "$headers = @{}\n";
        scriptfile << "$column = 2\n";
        scriptfile << "while ($true) {\n";
        scriptfile << "    $header = $sheet.Cells.Item(1, $column).Value2\n";
        scriptfile << "    if ($null -eq $header -or [string]::IsNullOrWhiteSpace([string]$header)) { break }\n";
        scriptfile << "    $normalizedHeader = Normalize-StepLabel([string]$header)\n";
        scriptfile << "    if (-not $headers.ContainsKey($normalizedHeader)) { $headers[$normalizedHeader] = New-Object System.Collections.ArrayList }\n";
        scriptfile << "    [void]$headers[$normalizedHeader].Add($column)\n";
        scriptfile << "    $column++\n";
        scriptfile << "}\n";
        scriptfile << "$usedColumns = New-Object System.Collections.Generic.HashSet[int]\n";
        scriptfile << "for ($index = 1; $index -lt $lines.Count; $index++) {\n";
        scriptfile << "    if ([string]::IsNullOrWhiteSpace($lines[$index])) { continue }\n";
        scriptfile << "    $parts = $lines[$index].Split(\"`t\", 2)\n";
        scriptfile << "    if ($parts.Count -ne 2) { continue }\n";
        scriptfile << "    $normalizedLabel = Normalize-StepLabel($parts[0])\n";
        scriptfile << "    if (-not $headers.ContainsKey($normalizedLabel)) { continue }\n";
        scriptfile << "    foreach ($targetColumn in $headers[$normalizedLabel]) {\n";
        scriptfile << "        if ($usedColumns.Add([int]$targetColumn)) {\n";
        scriptfile << "            $sheet.Cells.Item($row, [int]$targetColumn).Value2 = [double]::Parse($parts[1], $culture)\n";
        scriptfile << "            break\n";
        scriptfile << "        }\n";
        scriptfile << "    }\n";
        scriptfile << "}\n";
        scriptfile << "$workbook.Save()\n";
        scriptfile << "$workbook.Close($true)\n";
        scriptfile << "$excel.Quit()\n";
        scriptfile << "[void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($sheet)\n";
        scriptfile << "[void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($workbook)\n";
        scriptfile << "[void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($excel)\n";
    }

    wstring commandline = L"powershell.exe -NoProfile -ExecutionPolicy Bypass -File \"";
    commandline += scriptpath.wstring();
    commandline += L"\" -WorkbookPath \"";
    commandline += workbookpath.wstring();
    commandline += L"\" -DataPath \"";
    commandline += datapath.wstring();
    commandline += L"\"";

    runprocessandwait(commandline);
}
}

bool beginworldgenstep(const char* label)
{
    updatereport(label);
    return true;
}

void beginworldgendebugrun(long seed, const WorldGenerationDebugOptions* options)
{
    activeworldgenerationseed = seed;
    activeworldgenerationoptions = options;
    activeworldgenerationtimings.clear();
}

void endworldgendebugrun()
{
    if (activeworldgenerationoptions != nullptr && activeworldgenerationoptions->logToProfilingWorkbook)
        appendprofilingrow(activeworldgenerationseed, activeworldgenerationtimings);

    activeworldgenerationtimings.clear();
    activeworldgenerationoptions = nullptr;
    activeworldgenerationseed = 0;
    worldgenerationvisualizationcallback = nullptr;
}

void onworldgenstepcompleted(const string& label, double elapsedms)
{
    if (activeworldgenerationoptions == nullptr)
        return;

    activeworldgenerationtimings.emplace_back(label, elapsedms);
}

void setworldgenvisualizationcallback(function<void()> callback)
{
    worldgenerationvisualizationcallback = std::move(callback);
}

void clearworldgenvisualizationcallback()
{
    worldgenerationvisualizationcallback = nullptr;
}

bool hasworldgenvisualizationcallback()
{
    return static_cast<bool>(worldgenerationvisualizationcallback);
}

void requestworldgenvisualization()
{
    if (worldgenerationvisualizationcallback)
        worldgenerationvisualizationcallback();
}

bool isworldgendebugrunactive()
{
    return activeworldgenerationoptions != nullptr;
}

long worldgenerationdebugseed()
{
    return activeworldgenerationseed;
}

int platetectonicscyclecount()
{
    if (activeworldgenerationoptions == nullptr)
        return static_cast<int>(tuning::terrain::platetectonics::cycleCount);

    if (activeworldgenerationoptions->plateTectonicsCycleCount < 1)
        return 1;

    return activeworldgenerationoptions->plateTectonicsCycleCount;
}

int platetectonicsplatecount()
{
    if (activeworldgenerationoptions == nullptr)
        return static_cast<int>(tuning::terrain::platetectonics::plateCount);

    if (activeworldgenerationoptions->plateTectonicsPlateCount < 1)
        return 1;

    return activeworldgenerationoptions->plateTectonicsPlateCount;
}

#include "climate_validation.cpp"
