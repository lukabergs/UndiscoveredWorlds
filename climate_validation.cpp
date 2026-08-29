#include "app_environment.hpp"
#include "functions.hpp"

#include <windows.h>

#include <array>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace
{
struct zonalstats
{
    int cells = 0;
    int landcells = 0;
    int oceancells = 0;
    double annualrain = 0.0;
    double januaryrain = 0.0;
    double julyrain = 0.0;
    double landannualrain = 0.0;
    double oceanannualrain = 0.0;
    double januarypressure = 0.0;
    double julypressure = 0.0;
    double januaryuwind = 0.0;
    double januaryvwind = 0.0;
    double julyuwind = 0.0;
    double julyvwind = 0.0;
    double januarysst = 0.0;
    double julysst = 0.0;
    double januarycurrentu = 0.0;
    double januarycurrentv = 0.0;
    double julycurrentu = 0.0;
    double julycurrentv = 0.0;
    double januaryevaporation = 0.0;
    double julyevaporation = 0.0;
    double januarymoisture = 0.0;
    double julymoisture = 0.0;
};

struct comparisonmetrics
{
    bool referencefound = false;
    bool dimensionsmatch = false;
    int comparedcells = 0;
    double simulatedmean = 0.0;
    double referencemean = 0.0;
    double meanbias = 0.0;
    double meanabsoluteerror = 0.0;
    double rmse = 0.0;
    double correlation = 0.0;
    double tropicalmeanbias = 0.0;
};

struct climatespatialmetrics
{
    bool referencefound = false;
    bool dimensionsmatch = false;
    long long comparedcells = 0;
    long long exactmatches = 0;
    long long groupmatches = 0;
    double exactaccuracy = 0.0;
    double groupaccuracy = 0.0;
    double kappa = 0.0;
};

constexpr int benchmarkcolourdistancelimit = 65;
constexpr int benchmarkcolourdistancelimitsquared = benchmarkcolourdistancelimit * benchmarkcolourdistancelimit;

constexpr array<array<int, 3>, 31> benchmarkclimatecolours =
{
    array<int, 3>{ 0, 0, 254 }, { 1, 119, 255 }, { 70, 169, 250 }, { 70, 169, 250 },
    { 249, 15, 0 }, { 251, 150, 149 }, { 245, 163, 1 }, { 254, 219, 99 },
    { 255, 255, 0 }, { 198, 199, 1 }, { 184, 184, 114 }, { 138, 255, 162 },
    { 86, 199, 112 }, { 30, 150, 66 }, { 192, 254, 109 }, { 76, 255, 93 },
    { 19, 203, 74 }, { 255, 8, 245 }, { 204, 3, 192 }, { 154, 51, 144 },
    { 153, 100, 146 }, { 172, 178, 249 }, { 91, 121, 213 }, { 78, 83, 175 },
    { 54, 3, 130 }, { 0, 255, 245 }, { 32, 200, 250 }, { 0, 126, 125 },
    { 0, 69, 92 }, { 178, 178, 178 }, { 104, 104, 104 }
};

sf::Color benchmarkclimatecolour(short climate)
{
    if (climate < 1 || climate > static_cast<short>(benchmarkclimatecolours.size()))
        return sf::Color::Black;

    const auto& colour = benchmarkclimatecolours[climate - 1];
    return sf::Color(colour[0], colour[1], colour[2]);
}

short nearestbenchmarkclimate(const sf::Color& pixel, int& distancesquared)
{
    short nearest = 0;
    distancesquared = 3 * 255 * 255;

    for (short climate = 1; climate <= static_cast<short>(benchmarkclimatecolours.size()); climate++)
    {
        const auto& colour = benchmarkclimatecolours[climate - 1];
        const int red = static_cast<int>(pixel.r) - colour[0];
        const int green = static_cast<int>(pixel.g) - colour[1];
        const int blue = static_cast<int>(pixel.b) - colour[2];
        const int candidate = red * red + green * green + blue * blue;

        if (candidate < distancesquared)
        {
            nearest = climate;
            distancesquared = candidate;
        }
    }

    return nearest;
}

short comparableclimate(short climate)
{
    // Aw and As share a reference-map colour, so the image cannot distinguish them.
    return climate == 4 ? 3 : climate;
}

int climatemajorgroup(short climate)
{
    if (climate >= 1 && climate <= 4) return 0;
    if (climate >= 5 && climate <= 8) return 1;
    if (climate >= 9 && climate <= 17) return 2;
    if (climate >= 18 && climate <= 29) return 3;
    if (climate >= 30 && climate <= 31) return 4;
    return -1;
}

climatespatialmetrics compareclimatespatially(planet& world)
{
    climatespatialmetrics metrics;
    sf::Image reference;
    const filesystem::path referencepath = getappenvironment().earthKoppenImagePath;

    if (reference.loadFromFile(referencepath.string()) == false)
        return metrics;

    metrics.referencefound = true;
    const int width = world.width();
    const int height = world.height();
    const sf::Vector2u referencesize = reference.getSize();

    if (referencesize.x != static_cast<unsigned int>(width + 1) || referencesize.y != static_cast<unsigned int>(height + 1))
        return metrics;

    metrics.dimensionsmatch = true;
    array<array<long long, 31>, 31> confusion{};

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 1)
                continue;

            short simulated = comparableclimate(static_cast<short>(world.climate(x, y)));

            if (simulated < 1 || simulated > 31)
                continue;

            int distancesquared = 0;
            short expected = nearestbenchmarkclimate(reference.getPixel(x, y), distancesquared);

            if (distancesquared > benchmarkcolourdistancelimitsquared)
                continue;

            expected = comparableclimate(expected);
            metrics.comparedcells++;
            confusion[simulated - 1][expected - 1]++;

            if (simulated == expected)
                metrics.exactmatches++;

            if (climatemajorgroup(simulated) == climatemajorgroup(expected))
                metrics.groupmatches++;
        }
    }

    if (metrics.comparedcells <= 0)
        return metrics;

    metrics.exactaccuracy = static_cast<double>(metrics.exactmatches) / static_cast<double>(metrics.comparedcells);
    metrics.groupaccuracy = static_cast<double>(metrics.groupmatches) / static_cast<double>(metrics.comparedcells);

    array<long long, 31> simulatedtotals{};
    array<long long, 31> referencetotals{};

    for (int simulated = 0; simulated < 31; simulated++)
    {
        for (int expected = 0; expected < 31; expected++)
        {
            simulatedtotals[simulated] += confusion[simulated][expected];
            referencetotals[expected] += confusion[simulated][expected];
        }
    }

    double chanceagreement = 0.0;
    const double comparedsquared = static_cast<double>(metrics.comparedcells) * static_cast<double>(metrics.comparedcells);

    for (int climate = 0; climate < 31; climate++)
        chanceagreement += static_cast<double>(simulatedtotals[climate]) * static_cast<double>(referencetotals[climate]) / comparedsquared;

    if (chanceagreement < 1.0)
        metrics.kappa = (metrics.exactaccuracy - chanceagreement) / (1.0 - chanceagreement);

    return metrics;
}

double safeaverage(double total, int count)
{
    if (count <= 0)
        return 0.0;

    return total / static_cast<double>(count);
}

short climatefromcode(const string& code)
{
    for (short candidate = 1; candidate <= 31; candidate++)
    {
        if (getclimatecode(candidate) == code)
            return candidate;
    }

    return 0;
}

string csvescape(const string& value)
{
    if (value.find(',') == string::npos && value.find('"') == string::npos)
        return value;

    string escaped = "\"";

    for (char ch : value)
    {
        if (ch == '"')
            escaped += "\"\"";
        else
            escaped += ch;
    }

    escaped += '"';
    return escaped;
}

float latitudeforrow(int row, int height)
{
    if (height <= 0)
        return 0.0f;

    return 90.0f - (180.0f * static_cast<float>(row) / static_cast<float>(height));
}

filesystem::path climatevalidationoutputdirectory()
{
    filesystem::path outputroot = getappenvironment().profilingWorkbookPath.parent_path();

    if (outputroot.empty())
        outputroot = filesystem::current_path();

    outputroot /= "validation";
    outputroot /= "seed_" + to_string(worldgenerationdebugseed());
    filesystem::create_directories(outputroot);

    return outputroot;
}

vector<string> orderedclimatecodes()
{
    vector<string> codes;
    codes.reserve(31);

    for (short climate = 1; climate <= 31; climate++)
        codes.push_back(getclimatecode(climate));

    return codes;
}

vector<long long> collectsimulatedclimatecounts(planet& world)
{
    vector<long long> counts(31, 0);
    const int width = world.width();
    const int height = world.height();

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            if (world.sea(x, y) == 1)
                continue;

            const short climate = static_cast<short>(world.climate(x, y));

            if (climate >= 1 && climate <= 31)
                counts[climate - 1]++;
        }
    }

    return counts;
}

vector<long long> referenceclimatecounts()
{
    return
    {
        20086, 13831, 24368, 24368, 64970, 24332, 24110, 29277,
        5124, 3226, 6, 11798, 4443, 6, 19304, 12052, 17,
        758, 2188, 9114, 65, 3831, 6466, 12105, 1494,
        6413, 40215, 86356, 757, 48284, 215627
    };
}

string jsonescape(const string& value)
{
    ostringstream escaped;
    escaped << hex << setfill('0');

    for (const unsigned char ch : value)
    {
        switch (ch)
        {
        case '"': escaped << "\\\""; break;
        case '\\': escaped << "\\\\"; break;
        case '\b': escaped << "\\b"; break;
        case '\f': escaped << "\\f"; break;
        case '\n': escaped << "\\n"; break;
        case '\r': escaped << "\\r"; break;
        case '\t': escaped << "\\t"; break;
        default:
            if (ch < 0x20)
                escaped << "\\u" << setw(4) << static_cast<int>(ch);
            else
                escaped << static_cast<char>(ch);
            break;
        }
    }

    return escaped.str();
}

string climatebenchmarktimestamp()
{
    const time_t now = chrono::system_clock::to_time_t(chrono::system_clock::now());
    tm localtime{};
    localtime_s(&localtime, &now);

    ostringstream timestamp;
    timestamp << put_time(&localtime, "%Y%m%d%H%M%S");
    return timestamp.str();
}

int nextclimatebenchmarkrunid()
{
    const filesystem::path logpath = getappenvironment().climateBenchmarkRunLogPath;
    ifstream logfile(logpath);
    int maximumid = 1;

    if (!logfile.is_open())
        return maximumid + 1;

    const regex idpattern(R"(^\s*"id"\s*:\s*([0-9]+)\s*,?\s*$)");
    string line;

    while (getline(logfile, line))
    {
        smatch found;

        if (!regex_match(line, found, idpattern))
            continue;

        try
        {
            maximumid = max(maximumid, stoi(found[1].str()));
        }
        catch (const exception&)
        {
            return -1;
        }
    }

    return maximumid + 1;
}

bool appendclimatebenchmarkrunlog(
    int runid,
    const string& timestamp,
    const string& information,
    double weightedrelativeerror,
    const climatespatialmetrics& spatial)
{
    const filesystem::path logpath = getappenvironment().climateBenchmarkRunLogPath;
    string content;

    {
        ifstream logfile(logpath);

        if (logfile.is_open())
            content.assign(istreambuf_iterator<char>(logfile), istreambuf_iterator<char>());
    }

    if (content.find_first_not_of(" \t\r\n") == string::npos)
        content = "{\n  \"runs\": []\n}\n";

    const size_t runskey = content.find("\"runs\"");
    const size_t arraybegin = runskey == string::npos ? string::npos : content.find('[', runskey);
    const size_t arrayend = arraybegin == string::npos ? string::npos : content.rfind(']');

    if (arraybegin == string::npos || arrayend == string::npos || arrayend < arraybegin)
        return false;

    const bool hasentries = content.find_first_not_of(" \t\r\n", arraybegin + 1) < arrayend;
    ostringstream entry;

    if (hasentries)
        entry << ',';

    entry
        << "\n    {\n"
        << "      \"id\": " << runid << ",\n"
        << "      \"datetime\": \"" << timestamp << "\",\n"
        << "      \"information\": \"" << jsonescape(information) << "\",\n"
        << "      \"metrics\": {\n"
        << fixed << setprecision(10)
        << "        \"weighted_relative_error\": " << weightedrelativeerror << ",\n"
        << "        \"spatial_compared_cells\": " << spatial.comparedcells << ",\n"
        << "        \"spatial_exact_accuracy\": " << spatial.exactaccuracy << ",\n"
        << "        \"spatial_group_accuracy\": " << spatial.groupaccuracy << ",\n"
        << "        \"spatial_kappa\": " << spatial.kappa << "\n"
        << "      }\n"
        << "    }\n  ";

    content.insert(arrayend, entry.str());

    if (logpath.has_parent_path())
        filesystem::create_directories(logpath.parent_path());

    filesystem::path temppath = logpath;
    temppath += ".tmp";

    {
        ofstream tempfile(temppath, ios::trunc);

        if (!tempfile.is_open())
            return false;

        tempfile << content;
    }

    if (MoveFileExW(temppath.wstring().c_str(), logpath.wstring().c_str(), MOVEFILE_REPLACE_EXISTING | MOVEFILE_WRITE_THROUGH) == FALSE)
    {
        error_code ignored;
        filesystem::remove(temppath, ignored);
        return false;
    }

    return true;
}

double saferelativeerror(long long simulated, long long reference)
{
    if (reference == 0)
        return simulated == 0 ? 0.0 : 1.0;

    return fabs(static_cast<double>(simulated - reference)) / fabs(static_cast<double>(reference));
}

bool runhiddenprocessandwait(const wstring& commandline)
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

bool updateclimatebenchmarkworkbook(int runid, const vector<string>& codes, const vector<long long>& simulationcounts)
{
    const filesystem::path workbookpath = filesystem::absolute(getappenvironment().climateWorkbookPath).lexically_normal();

    if (filesystem::exists(workbookpath) == false || simulationcounts.size() != codes.size())
        return false;

    const filesystem::path temproot = filesystem::temp_directory_path();
    const filesystem::path datapath = temproot / "uw_climate_benchmark_input.txt";
    const filesystem::path scriptpath = temproot / "uw_climate_benchmark_excel.ps1";
    const filesystem::path errorpath = temproot / "uw_climate_benchmark_excel_error.txt";

    {
        ofstream datafile(datapath);

        if (!datafile.is_open())
            return false;

        datafile << "ID";

        for (const string& code : codes)
            datafile << ',' << code;

        datafile << '\n';
        datafile << runid;

        for (const long long count : simulationcounts)
            datafile << ',' << count;

        datafile << '\n';
    }

    {
        ofstream scriptfile(scriptpath);

        if (!scriptfile.is_open())
            return false;

        scriptfile << "param([string]$WorkbookPath, [string]$DataPath, [string]$ErrorPath)\n";
        scriptfile << "$ErrorActionPreference = 'Stop'\n";
        scriptfile << "Remove-Item -LiteralPath $ErrorPath -ErrorAction SilentlyContinue\n";
        scriptfile << "$excel = $null\n";
        scriptfile << "$workbook = $null\n";
        scriptfile << "$sheet = $null\n";
        scriptfile << "try {\n";
        scriptfile << "$lines = Get-Content -Path $DataPath\n";
        scriptfile << "if ($lines.Count -lt 2) { throw 'Benchmark input is incomplete' }\n";
        scriptfile << "$headers = $lines[0].Split(',')\n";
        scriptfile << "$simulation = $lines[1].Split(',')\n";
        scriptfile << "$excel = New-Object -ComObject Excel.Application\n";
        scriptfile << "$excel.Visible = $false\n";
        scriptfile << "$excel.DisplayAlerts = $false\n";
        scriptfile << "$workbook = $excel.Workbooks.Open($WorkbookPath)\n";
        scriptfile << "if ($workbook.ReadOnly) { throw 'Climate benchmark workbook is open read-only; close it in Excel and retry' }\n";
        scriptfile << "$sheet = $workbook.Worksheets.Item('RAW_PIXELS')\n";
        scriptfile << "for ($index = 0; $index -lt $headers.Count; $index++) {\n";
        scriptfile << "    $actual = [string]$sheet.Cells.Item(1, $index + 1).Text\n";
        scriptfile << "    if ($actual -ne $headers[$index]) { throw \"Workbook header mismatch at column $($index + 1): expected '$($headers[$index])', found '$actual'\" }\n";
        scriptfile << "}\n";
        scriptfile << "$runId = [int]$simulation[0]\n";
        scriptfile << "if ($runId -lt 2) { throw 'Benchmark run ID must be at least 2' }\n";
        scriptfile << "$row = $runId + 2\n";
        scriptfile << "$existingId = $sheet.Cells.Item($row, 1).Value2\n";
        scriptfile << "if ($null -ne $existingId -and -not [string]::IsNullOrWhiteSpace([string]$existingId) -and [int]$existingId -ne $runId) { throw \"Workbook row $row already belongs to run $existingId\" }\n";
        scriptfile << "for ($column = $simulation.Count + 1; $column -le 34; $column++) {\n";
        scriptfile << "    if (-not $sheet.Cells.Item($row, $column).HasFormula) { throw \"RAW_PIXELS formula missing at row $row, column $column\" }\n";
        scriptfile << "}\n";
        scriptfile << "for ($index = 0; $index -lt $simulation.Count; $index++) {\n";
        scriptfile << "    $sheet.Cells.Item($row, $index + 1).Value2 = [double]$simulation[$index]\n";
        scriptfile << "}\n";
        scriptfile << "$workbook.Save()\n";
        scriptfile << "}\n";
        scriptfile << "catch {\n";
        scriptfile << "    Set-Content -LiteralPath $ErrorPath -Value $_.Exception.Message\n";
        scriptfile << "    exit 1\n";
        scriptfile << "}\n";
        scriptfile << "finally {\n";
        scriptfile << "    if ($null -ne $sheet) { [void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($sheet) }\n";
        scriptfile << "    if ($null -ne $workbook) { try { $workbook.Close($false) } catch {}; [void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($workbook) }\n";
        scriptfile << "    if ($null -ne $excel) { try { $excel.Quit() } catch {}; [void][System.Runtime.Interopservices.Marshal]::ReleaseComObject($excel) }\n";
        scriptfile << "    [GC]::Collect()\n";
        scriptfile << "    [GC]::WaitForPendingFinalizers()\n";
        scriptfile << "}\n";
    }

    wstring commandline = L"powershell.exe -NoProfile -ExecutionPolicy Bypass -File \"";
    commandline += scriptpath.wstring();
    commandline += L"\" -WorkbookPath \"";
    commandline += workbookpath.wstring();
    commandline += L"\" -DataPath \"";
    commandline += datapath.wstring();
    commandline += L"\" -ErrorPath \"";
    commandline += errorpath.wstring();
    commandline += L"\"";

    const bool succeeded = runhiddenprocessandwait(commandline);

    if (!succeeded)
    {
        ifstream errorfile(errorpath);
        string errormessage;

        if (getline(errorfile, errormessage) && !errormessage.empty())
            cerr << "Climate benchmark workbook error: " << errormessage << '\n';
    }

    return succeeded;
}

bool exportclimatebenchmarkimages(planet& world, int runid)
{
    const AppEnvironmentConfig& appenv = getappenvironment();
    const filesystem::path outputdir = appenv.climateBenchmarkImageDirectory;
    const filesystem::path referencepath = appenv.earthKoppenImagePath;

    if (outputdir.empty())
    {
        cerr << "Climate benchmark image directory is not configured.\n";
        return false;
    }

    if (filesystem::exists(referencepath) == false)
    {
        cerr << "Climate benchmark reference image not found: " << referencepath.string() << '\n';
        return false;
    }

    error_code filesystemerror;
    filesystem::create_directories(outputdir, filesystemerror);

    if (filesystemerror)
    {
        cerr << "Failed to create climate benchmark image directory: " << filesystemerror.message() << '\n';
        return false;
    }

    sf::Image referenceimage;

    if (referenceimage.loadFromFile(referencepath.string()) == false)
    {
        cerr << "Failed to load climate benchmark reference image: " << referencepath.string() << '\n';
        return false;
    }

    const sf::Vector2u referencesize = referenceimage.getSize();

    for (unsigned int y = 0; y < referencesize.y; y++)
    {
        for (unsigned int x = 0; x < referencesize.x; x++)
        {
            int distancesquared = 0;
            const short climate = nearestbenchmarkclimate(referenceimage.getPixel(x, y), distancesquared);

            if (distancesquared <= benchmarkcolourdistancelimitsquared && climate != 31)
                referenceimage.setPixel(x, y, benchmarkclimatecolour(climate));
            else
                referenceimage.setPixel(x, y, sf::Color::Black);
        }
    }

    const filesystem::path benchmarkreferencepath = outputdir / "0.png";

    if (referenceimage.saveToFile(benchmarkreferencepath.string()) == false)
    {
        cerr << "Failed to save climate benchmark reference image: " << benchmarkreferencepath.string() << '\n';
        return false;
    }

    sf::Image simulatedimage;
    const int width = world.width();
    const int height = world.height();
    simulatedimage.create(width + 1, height + 1, sf::Color::Black);

    for (int y = 0; y <= height; y++)
    {
        for (int x = 0; x <= width; x++)
        {
            const short climate = static_cast<short>(world.climate(x, y));
            const bool water = world.sea(x, y) == 1 || world.truelake(x, y) != 0 || world.riftlakesurface(x, y) != 0;

            if (water == false && climate >= 1 && climate <= 30)
                simulatedimage.setPixel(x, y, benchmarkclimatecolour(climate));
        }
    }

    const filesystem::path simulatedpath = outputdir / (to_string(runid) + ".png");

    if (simulatedimage.saveToFile(simulatedpath.string()) == false)
    {
        cerr << "Failed to save simulated climate image: " << simulatedpath.string() << '\n';
        return false;
    }

    return true;
}

bool loadprecipitationgrid(const filesystem::path& filepath, vector<vector<double>>& grid)
{
    ifstream infile(filepath);

    if (!infile.is_open())
        return false;

    string line;

    if (!getline(infile, line))
        return false;

    while (getline(infile, line))
    {
        if (line.empty())
            continue;

        vector<double> rowvalues;
        string token;
        stringstream linestream(line);
        int column = 0;

        while (getline(linestream, token, ','))
        {
            if (column >= 2)
                rowvalues.push_back(stod(token));

            column++;
        }

        if (rowvalues.empty() == false)
            grid.push_back(rowvalues);
    }

    return grid.empty() == false;
}

comparisonmetrics compareannualprecipitation(const filesystem::path& outputdir, planet& world, const vector<zonalstats>& rows)
{
    comparisonmetrics metrics;
    const filesystem::path referencepath = getappenvironment().referencePrecipitationGridPath;
    const filesystem::path comparisonpath = outputdir / "annual_precipitation_comparison.txt";

    if (filesystem::exists(referencepath) == false)
    {
        ofstream comparisonfile(comparisonpath);

        if (comparisonfile.is_open())
            comparisonfile << "status=reference_not_found\nreference_grid_path=" << referencepath.string() << '\n';

        return metrics;
    }

    metrics.referencefound = true;

    vector<vector<double>> referencegrid;

    if (loadprecipitationgrid(referencepath, referencegrid) == false)
    {
        ofstream comparisonfile(comparisonpath);

        if (comparisonfile.is_open())
            comparisonfile << "status=reference_unreadable\nreference_grid_path=" << referencepath.string() << '\n';

        return metrics;
    }

    const int width = world.width();
    const int height = world.height();

    if (static_cast<int>(referencegrid.size()) != height + 1)
    {
        ofstream comparisonfile(comparisonpath);

        if (comparisonfile.is_open())
        {
            comparisonfile << "status=dimension_mismatch\n";
            comparisonfile << "reference_grid_path=" << referencepath.string() << '\n';
            comparisonfile << "expected_height=" << height + 1 << '\n';
            comparisonfile << "actual_height=" << referencegrid.size() << '\n';
        }

        return metrics;
    }

    for (const auto& row : referencegrid)
    {
        if (static_cast<int>(row.size()) != width + 1)
        {
            ofstream comparisonfile(comparisonpath);

            if (comparisonfile.is_open())
            {
                comparisonfile << "status=dimension_mismatch\n";
                comparisonfile << "reference_grid_path=" << referencepath.string() << '\n';
                comparisonfile << "expected_width=" << width + 1 << '\n';
                comparisonfile << "actual_width=" << row.size() << '\n';
            }

            return metrics;
        }
    }

    metrics.dimensionsmatch = true;

    double simulatedsum = 0.0;
    double referencesum = 0.0;
    double biassum = 0.0;
    double absoluteerrorsum = 0.0;
    double squarederrorsum = 0.0;
    double sumsim2 = 0.0;
    double sumref2 = 0.0;
    double sumcross = 0.0;
    double tropicalbiassum = 0.0;
    int tropicalcells = 0;

    ofstream zonalcomparisonfile(outputdir / "annual_precipitation_zonal_comparison.csv");

    if (zonalcomparisonfile.is_open())
    {
        zonalcomparisonfile << "y,latitude,simulated_mean_annual_rain,reference_mean_annual_rain,bias\n";
        zonalcomparisonfile << fixed << setprecision(4);
    }

    for (int y = 0; y <= height; y++)
    {
        double referencesumrow = 0.0;

        for (int x = 0; x <= width; x++)
        {
            const double simulated = static_cast<double>(world.averain(x, y));
            const double reference = referencegrid[y][x];
            const double diff = simulated - reference;

            metrics.comparedcells++;
            simulatedsum = simulatedsum + simulated;
            referencesum = referencesum + reference;
            biassum = biassum + diff;
            absoluteerrorsum = absoluteerrorsum + fabs(diff);
            squarederrorsum = squarederrorsum + diff * diff;
            sumsim2 = sumsim2 + simulated * simulated;
            sumref2 = sumref2 + reference * reference;
            sumcross = sumcross + simulated * reference;
            referencesumrow = referencesumrow + reference;

            if (fabs(latitudeforrow(y, height)) <= 30.0f)
            {
                tropicalbiassum = tropicalbiassum + diff;
                tropicalcells++;
            }
        }

        if (zonalcomparisonfile.is_open())
        {
            const double simulatedmeanrow = safeaverage(rows[y].annualrain, rows[y].cells);
            const double referencemeanrow = safeaverage(referencesumrow, width + 1);

            zonalcomparisonfile
                << y << ','
                << latitudeforrow(y, height) << ','
                << simulatedmeanrow << ','
                << referencemeanrow << ','
                << simulatedmeanrow - referencemeanrow << '\n';
        }
    }

    metrics.simulatedmean = safeaverage(simulatedsum, metrics.comparedcells);
    metrics.referencemean = safeaverage(referencesum, metrics.comparedcells);
    metrics.meanbias = safeaverage(biassum, metrics.comparedcells);
    metrics.meanabsoluteerror = safeaverage(absoluteerrorsum, metrics.comparedcells);
    metrics.rmse = sqrt(safeaverage(squarederrorsum, metrics.comparedcells));
    metrics.tropicalmeanbias = safeaverage(tropicalbiassum, tropicalcells);

    const double numerator = sumcross - (simulatedsum * referencesum / static_cast<double>(metrics.comparedcells));
    const double simulatedvariance = sumsim2 - (simulatedsum * simulatedsum / static_cast<double>(metrics.comparedcells));
    const double referencevariance = sumref2 - (referencesum * referencesum / static_cast<double>(metrics.comparedcells));
    const double denominator = sqrt(max(0.0, simulatedvariance) * max(0.0, referencevariance));

    if (denominator > 0.0)
        metrics.correlation = numerator / denominator;

    ofstream comparisonfile(comparisonpath);

    if (comparisonfile.is_open())
    {
        comparisonfile << "status=ok\n";
        comparisonfile << fixed << setprecision(6);
        comparisonfile << "reference_grid_path=" << referencepath.string() << '\n';
        comparisonfile << "compared_cells=" << metrics.comparedcells << '\n';
        comparisonfile << "simulated_mean=" << metrics.simulatedmean << '\n';
        comparisonfile << "reference_mean=" << metrics.referencemean << '\n';
        comparisonfile << "mean_bias=" << metrics.meanbias << '\n';
        comparisonfile << "mean_absolute_error=" << metrics.meanabsoluteerror << '\n';
        comparisonfile << "rmse=" << metrics.rmse << '\n';
        comparisonfile << "correlation=" << metrics.correlation << '\n';
        comparisonfile << "tropical_mean_bias=" << metrics.tropicalmeanbias << '\n';
    }

    return metrics;
}

void writeprecipitationgrid(const filesystem::path& filepath, planet& world, int mode)
{
    ofstream outfile(filepath);

    if (!outfile.is_open())
        return;

    const int width = world.width();
    const int height = world.height();

    outfile << "y,latitude";

    for (int x = 0; x <= width; x++)
        outfile << ",x" << x;

    outfile << '\n';
    outfile << fixed << setprecision(2);

    for (int y = 0; y <= height; y++)
    {
        outfile << y << ',' << latitudeforrow(y, height);

        for (int x = 0; x <= width; x++)
        {
            int value = 0;

            if (mode == 0)
                value = world.averain(x, y);

            if (mode == 1)
                value = world.janrain(x, y);

            if (mode == 2)
                value = world.julrain(x, y);

            if (mode == 3)
                value = world.sea(x, y) == 0 ? 1 : 0;

            outfile << ',' << value;
        }

        outfile << '\n';
    }
}
}

bool recordclimatebenchmarkrun(planet& world, const string& information, bool updateworkbook, int* runid)
{
    const int nextrunid = nextclimatebenchmarkrunid();
    const string timestamp = climatebenchmarktimestamp();
    const vector<string> codes = orderedclimatecodes();
    const vector<long long> simulationcounts = collectsimulatedclimatecounts(world);
    const vector<long long> referencecounts = referenceclimatecounts();
    const climatespatialmetrics spatial = compareclimatespatially(world);
    long long totalabsoluteerror = 0;
    long long totalreference = 0;

    for (size_t index = 0; index < simulationcounts.size() && index < referencecounts.size(); index++)
    {
        totalabsoluteerror += llabs(simulationcounts[index] - referencecounts[index]);
        totalreference += referencecounts[index];
    }

    const double weightedrelativeerror = totalreference > 0
        ? static_cast<double>(totalabsoluteerror) / static_cast<double>(totalreference)
        : 0.0;

    if (nextrunid < 2)
        return false;

    if (updateworkbook && updateclimatebenchmarkworkbook(nextrunid, codes, simulationcounts) == false)
    {
        cerr << "Failed to update climate benchmark workbook.\n";
        return false;
    }

    if (exportclimatebenchmarkimages(world, nextrunid) == false)
        return false;

    if (appendclimatebenchmarkrunlog(nextrunid, timestamp, information, weightedrelativeerror, spatial) == false)
    {
        cerr << "Failed to update climate benchmark run log.\n";
        return false;
    }

    if (runid != nullptr)
        *runid = nextrunid;

    return true;
}

void printclimaterelativeerrorreport(planet& world)
{
    const vector<string> codes = orderedclimatecodes();
    const vector<long long> references = referenceclimatecounts();
    const vector<long long> simulated = collectsimulatedclimatecounts(world);

    if (references.size() != codes.size() || simulated.size() != codes.size())
        return;

    double meanrelativeerror = 0.0;
    long long totalsimulated = 0;
    long long totalreference = 0;
    long long totalabsoluteerror = 0;
    double maxrelativeerror = -1.0;
    string maxcode;
    double adjustedmeanrelativeerror = 0.0;
    double adjustedmaxrelativeerror = -1.0;
    string adjustedmaxcode;
    int adjustedcount = 0;

    cout << "Climate relative error (raw pixels):" << '\n';
    cout << fixed << setprecision(4);

    for (size_t index = 0; index < codes.size(); index++)
    {
        const double relativeerror = saferelativeerror(simulated[index], references[index]);
        const long long absoluteerror = llabs(simulated[index] - references[index]);

        cout
            << codes[index]
            << " simulated=" << simulated[index]
            << " reference=" << references[index]
            << " relative_error=" << relativeerror
            << '\n';

        meanrelativeerror += relativeerror;
        totalsimulated += simulated[index];
        totalreference += references[index];
        totalabsoluteerror += absoluteerror;

        if (relativeerror > maxrelativeerror)
        {
            maxrelativeerror = relativeerror;
            maxcode = codes[index];
        }

        if (codes[index] == "Aw" || codes[index] == "As")
            continue;

        adjustedmeanrelativeerror += relativeerror;
        adjustedcount++;

        if (relativeerror > adjustedmaxrelativeerror)
        {
            adjustedmaxrelativeerror = relativeerror;
            adjustedmaxcode = codes[index];
        }
    }

    meanrelativeerror /= static_cast<double>(codes.size());

    const double awasrelativeerror = saferelativeerror(simulated[2] + simulated[3], references[2] + references[3]);
    adjustedmeanrelativeerror += awasrelativeerror;
    adjustedcount++;

    if (awasrelativeerror > adjustedmaxrelativeerror)
    {
        adjustedmaxrelativeerror = awasrelativeerror;
        adjustedmaxcode = "Aw+As";
    }

    adjustedmeanrelativeerror /= static_cast<double>(adjustedcount);

    const double weightedrelativeerror = totalreference > 0
        ? static_cast<double>(totalabsoluteerror) / static_cast<double>(totalreference)
        : 0.0;

    cout << "Aw+As combined"
         << " simulated=" << simulated[2] + simulated[3]
         << " reference=" << references[2] + references[3]
         << " relative_error=" << awasrelativeerror
         << '\n';
    cout << "Climate relative error summary:" << '\n';
    cout << "mean_relative_error=" << meanrelativeerror << '\n';
    cout << "mean_relative_error_adjusted=" << adjustedmeanrelativeerror << '\n';
    cout << "weighted_relative_error=" << weightedrelativeerror << '\n';
    cout << "max_relative_error=" << maxrelativeerror << " (" << maxcode << ")" << '\n';
    cout << "max_relative_error_adjusted=" << adjustedmaxrelativeerror << " (" << adjustedmaxcode << ")" << '\n';
    cout << "simulated_land_total=" << totalsimulated << '\n';
    cout << "reference_land_total=" << totalreference << '\n';

    const climatespatialmetrics spatial = compareclimatespatially(world);

    cout << "Climate spatial agreement:" << '\n';

    if (spatial.referencefound == false)
    {
        cout << "spatial_status=reference_not_found" << '\n';
        return;
    }

    if (spatial.dimensionsmatch == false)
    {
        cout << "spatial_status=dimension_mismatch" << '\n';
        return;
    }

    cout << "spatial_status=ok" << '\n';
    cout << "spatial_compared_cells=" << spatial.comparedcells << '\n';
    cout << "spatial_exact_accuracy=" << spatial.exactaccuracy << '\n';
    cout << "spatial_group_accuracy=" << spatial.groupaccuracy << '\n';
    cout << "spatial_kappa=" << spatial.kappa << '\n';
}

void exportclimatevalidationreport(planet& world)
{
    if (isworldgendebugrunactive() == false)
        return;

    const int width = world.width();
    const int height = world.height();
    const filesystem::path outputdir = climatevalidationoutputdirectory();

    vector<zonalstats> rows(height + 1);
    map<string, int> climatecounts;

    double globalannualrain = 0.0;
    double landannualrain = 0.0;
    double oceanannualrain = 0.0;
    double northernannualrain = 0.0;
    double southernannualrain = 0.0;
    int globalcells = 0;
    int landcells = 0;
    int oceancells = 0;
    int northerncells = 0;
    int southerncells = 0;

    for (int y = 0; y <= height; y++)
    {
        zonalstats& row = rows[y];

        for (int x = 0; x <= width; x++)
        {
            const bool sea = world.sea(x, y) == 1;
            const int annualrain = world.averain(x, y);
            const int januaryrain = world.janrain(x, y);
            const int julyrain = world.julrain(x, y);

            row.cells++;
            row.annualrain = row.annualrain + annualrain;
            row.januaryrain = row.januaryrain + januaryrain;
            row.julyrain = row.julyrain + julyrain;
            row.januarypressure = row.januarypressure + world.seasonalpressure(seasonjanuary, x, y);
            row.julypressure = row.julypressure + world.seasonalpressure(seasonjuly, x, y);
            row.januaryuwind = row.januaryuwind + world.seasonaluwind(seasonjanuary, x, y);
            row.januaryvwind = row.januaryvwind + world.seasonalvwind(seasonjanuary, x, y);
            row.julyuwind = row.julyuwind + world.seasonaluwind(seasonjuly, x, y);
            row.julyvwind = row.julyvwind + world.seasonalvwind(seasonjuly, x, y);
            row.januaryevaporation = row.januaryevaporation + world.seasonalevaporation(seasonjanuary, x, y);
            row.julyevaporation = row.julyevaporation + world.seasonalevaporation(seasonjuly, x, y);
            row.januarymoisture = row.januarymoisture + world.seasonalmoisture(seasonjanuary, x, y);
            row.julymoisture = row.julymoisture + world.seasonalmoisture(seasonjuly, x, y);

            globalcells++;
            globalannualrain = globalannualrain + annualrain;

            if (y < height / 2)
            {
                northernannualrain = northernannualrain + annualrain;
                northerncells++;
            }
            else if (y > height / 2)
            {
                southernannualrain = southernannualrain + annualrain;
                southerncells++;
            }

            if (sea)
            {
                row.oceancells++;
                row.oceanannualrain = row.oceanannualrain + annualrain;
                row.januarysst = row.januarysst + world.seasonalsst(seasonjanuary, x, y);
                row.julysst = row.julysst + world.seasonalsst(seasonjuly, x, y);
                row.januarycurrentu = row.januarycurrentu + world.seasonalcurrentu(seasonjanuary, x, y);
                row.januarycurrentv = row.januarycurrentv + world.seasonalcurrentv(seasonjanuary, x, y);
                row.julycurrentu = row.julycurrentu + world.seasonalcurrentu(seasonjuly, x, y);
                row.julycurrentv = row.julycurrentv + world.seasonalcurrentv(seasonjuly, x, y);

                oceancells++;
                oceanannualrain = oceanannualrain + annualrain;
            }
            else
            {
                row.landcells++;
                row.landannualrain = row.landannualrain + annualrain;

                landcells++;
                landannualrain = landannualrain + annualrain;

                const short climate = static_cast<short>(world.climate(x, y));

                if (climate > 0)
                    climatecounts[getclimatecode(climate)]++;
            }
        }
    }

    ofstream zonalfile(outputdir / "precipitation_zonal.csv");

    if (zonalfile.is_open())
    {
        zonalfile << "y,latitude,cells,land_cells,ocean_cells,mean_annual_rain,mean_jan_rain,mean_jul_rain,land_mean_annual_rain,ocean_mean_annual_rain,mean_jan_pressure,mean_jul_pressure,mean_jan_u_wind,mean_jan_v_wind,mean_jul_u_wind,mean_jul_v_wind,ocean_mean_jan_sst,ocean_mean_jul_sst,ocean_mean_jan_current_u,ocean_mean_jan_current_v,ocean_mean_jul_current_u,ocean_mean_jul_current_v,mean_jan_evaporation,mean_jul_evaporation,mean_jan_moisture,mean_jul_moisture\n";
        zonalfile << fixed << setprecision(4);

        for (int y = 0; y <= height; y++)
        {
            const zonalstats& row = rows[y];

            zonalfile
                << y << ','
                << latitudeforrow(y, height) << ','
                << row.cells << ','
                << row.landcells << ','
                << row.oceancells << ','
                << safeaverage(row.annualrain, row.cells) << ','
                << safeaverage(row.januaryrain, row.cells) << ','
                << safeaverage(row.julyrain, row.cells) << ','
                << safeaverage(row.landannualrain, row.landcells) << ','
                << safeaverage(row.oceanannualrain, row.oceancells) << ','
                << safeaverage(row.januarypressure, row.cells) << ','
                << safeaverage(row.julypressure, row.cells) << ','
                << safeaverage(row.januaryuwind, row.cells) << ','
                << safeaverage(row.januaryvwind, row.cells) << ','
                << safeaverage(row.julyuwind, row.cells) << ','
                << safeaverage(row.julyvwind, row.cells) << ','
                << safeaverage(row.januarysst, row.oceancells) << ','
                << safeaverage(row.julysst, row.oceancells) << ','
                << safeaverage(row.januarycurrentu, row.oceancells) << ','
                << safeaverage(row.januarycurrentv, row.oceancells) << ','
                << safeaverage(row.julycurrentu, row.oceancells) << ','
                << safeaverage(row.julycurrentv, row.oceancells) << ','
                << safeaverage(row.januaryevaporation, row.cells) << ','
                << safeaverage(row.julyevaporation, row.cells) << ','
                << safeaverage(row.januarymoisture, row.cells) << ','
                << safeaverage(row.julymoisture, row.cells) << '\n';
        }
    }

    int annualitczrow = height / 2;
    int januaryitczrow = height / 2;
    int julyitczrow = height / 2;
    int northdryrow = height / 4;
    int southdryrow = (height * 3) / 4;
    double annualitczrain = -1.0;
    double januaryitczrain = -1.0;
    double julyitczrain = -1.0;
    double northdryrain = 1e30;
    double southdryrain = 1e30;

    for (int y = 0; y <= height; y++)
    {
        const float latitude = latitudeforrow(y, height);
        const zonalstats& row = rows[y];
        const double meanannualrain = safeaverage(row.annualrain, row.cells);
        const double meanjanuaryrain = safeaverage(row.januaryrain, row.cells);
        const double meanjulyrain = safeaverage(row.julyrain, row.cells);

        if (fabs(latitude) <= 30.0f)
        {
            if (meanannualrain > annualitczrain)
            {
                annualitczrain = meanannualrain;
                annualitczrow = y;
            }

            if (meanjanuaryrain > januaryitczrain)
            {
                januaryitczrain = meanjanuaryrain;
                januaryitczrow = y;
            }

            if (meanjulyrain > julyitczrain)
            {
                julyitczrain = meanjulyrain;
                julyitczrow = y;
            }
        }

        if (latitude >= 10.0f && latitude <= 45.0f && meanannualrain < northdryrain)
        {
            northdryrain = meanannualrain;
            northdryrow = y;
        }

        if (latitude <= -10.0f && latitude >= -45.0f && meanannualrain < southdryrain)
        {
            southdryrain = meanannualrain;
            southdryrow = y;
        }
    }

    const comparisonmetrics comparison = compareannualprecipitation(outputdir, world, rows);
    ofstream summaryfile(outputdir / "precipitation_summary.txt");

    if (summaryfile.is_open())
    {
        summaryfile << fixed << setprecision(4);
        summaryfile << "seed=" << worldgenerationdebugseed() << '\n';
        summaryfile << "width=" << width << '\n';
        summaryfile << "height=" << height << '\n';
        summaryfile << "global_mean_annual_precip=" << safeaverage(globalannualrain, globalcells) << '\n';
        summaryfile << "land_mean_annual_precip=" << safeaverage(landannualrain, landcells) << '\n';
        summaryfile << "ocean_mean_annual_precip=" << safeaverage(oceanannualrain, oceancells) << '\n';
        summaryfile << "northern_mean_annual_precip=" << safeaverage(northernannualrain, northerncells) << '\n';
        summaryfile << "southern_mean_annual_precip=" << safeaverage(southernannualrain, southerncells) << '\n';
        summaryfile << "annual_itcz_latitude=" << latitudeforrow(annualitczrow, height) << '\n';
        summaryfile << "january_itcz_latitude=" << latitudeforrow(januaryitczrow, height) << '\n';
        summaryfile << "july_itcz_latitude=" << latitudeforrow(julyitczrow, height) << '\n';
        summaryfile << "north_subtropical_dry_latitude=" << latitudeforrow(northdryrow, height) << '\n';
        summaryfile << "south_subtropical_dry_latitude=" << latitudeforrow(southdryrow, height) << '\n';
        summaryfile << "annual_itcz_zonal_precip=" << annualitczrain << '\n';
        summaryfile << "north_subtropical_dry_zonal_precip=" << northdryrain << '\n';
        summaryfile << "south_subtropical_dry_zonal_precip=" << southdryrain << '\n';
        summaryfile << "grid_files=annual_precipitation_grid.csv,january_precipitation_grid.csv,july_precipitation_grid.csv,land_mask_grid.csv\n";
        summaryfile << "reference_grid_path=" << getappenvironment().referencePrecipitationGridPath.string() << '\n';
        summaryfile << "reference_found=" << (comparison.referencefound ? 1 : 0) << '\n';
        summaryfile << "reference_dimensions_match=" << (comparison.dimensionsmatch ? 1 : 0) << '\n';

        if (comparison.dimensionsmatch)
        {
            summaryfile << "reference_compared_cells=" << comparison.comparedcells << '\n';
            summaryfile << "reference_mean_bias=" << comparison.meanbias << '\n';
            summaryfile << "reference_mae=" << comparison.meanabsoluteerror << '\n';
            summaryfile << "reference_rmse=" << comparison.rmse << '\n';
            summaryfile << "reference_correlation=" << comparison.correlation << '\n';
            summaryfile << "reference_tropical_mean_bias=" << comparison.tropicalmeanbias << '\n';
        }
    }

    ofstream climatefile(outputdir / "climate_counts.csv");

    if (climatefile.is_open())
    {
        climatefile << "code,name,cells\n";

        for (const auto& entry : climatecounts)
        {
            const short climate = climatefromcode(entry.first);

            climatefile
                << entry.first << ','
                << csvescape(getclimatename(climate)) << ','
                << entry.second << '\n';
        }
    }

    writeprecipitationgrid(outputdir / "annual_precipitation_grid.csv", world, 0);
    writeprecipitationgrid(outputdir / "january_precipitation_grid.csv", world, 1);
    writeprecipitationgrid(outputdir / "july_precipitation_grid.csv", world, 2);
    writeprecipitationgrid(outputdir / "land_mask_grid.csv", world, 3);
}
