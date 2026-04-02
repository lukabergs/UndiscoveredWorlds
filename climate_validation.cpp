#include "app_environment.hpp"
#include "functions.hpp"

#include <windows.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
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

bool updateclimatebenchmarkworkbook(const vector<string>& codes, const vector<long long>& simulationcounts)
{
    const filesystem::path workbookpath = getappenvironment().climateWorkbookPath;

    if (filesystem::exists(workbookpath) == false || simulationcounts.size() != codes.size())
        return false;

    const filesystem::path temproot = filesystem::temp_directory_path();
    const filesystem::path datapath = temproot / "uw_climate_benchmark_input.txt";
    const filesystem::path scriptpath = temproot / "uw_climate_benchmark_excel.ps1";

    {
        ofstream datafile(datapath);

        if (!datafile.is_open())
            return false;

        for (size_t index = 0; index < codes.size(); index++)
        {
            if (index > 0)
                datafile << ',';

            datafile << codes[index];
        }

        datafile << '\n';

        for (size_t index = 0; index < simulationcounts.size(); index++)
        {
            if (index > 0)
                datafile << ',';

            datafile << simulationcounts[index];
        }

        datafile << '\n';
    }

    {
        ofstream scriptfile(scriptpath);

        if (!scriptfile.is_open())
            return false;

        scriptfile << "param([string]$WorkbookPath, [string]$DataPath)\n";
        scriptfile << "$ErrorActionPreference = 'Stop'\n";
        scriptfile << "$lines = Get-Content -Path $DataPath\n";
        scriptfile << "if ($lines.Count -lt 2) { exit 1 }\n";
        scriptfile << "$headers = $lines[0].Split(',')\n";
        scriptfile << "$simulation = $lines[1].Split(',')\n";
        scriptfile << "$excel = New-Object -ComObject Excel.Application\n";
        scriptfile << "$excel.Visible = $false\n";
        scriptfile << "$excel.DisplayAlerts = $false\n";
        scriptfile << "$workbook = $excel.Workbooks.Open($WorkbookPath)\n";
        scriptfile << "$sheet = $workbook.Worksheets.Item('RAW_PIXELS')\n";
        scriptfile << "for ($index = 0; $index -lt $headers.Count; $index++) {\n";
        scriptfile << "    $actual = [string]$sheet.Cells.Item(1, $index + 1).Text\n";
        scriptfile << "    if ($actual -ne $headers[$index]) { throw \"Workbook header mismatch at column $($index + 1): expected '$($headers[$index])', found '$actual'\" }\n";
        scriptfile << "}\n";
        scriptfile << "$row = 3\n";
        scriptfile << "while ($true) {\n";
        scriptfile << "    $value = $sheet.Cells.Item($row, 1).Value2\n";
        scriptfile << "    if ($null -eq $value -or [string]::IsNullOrWhiteSpace([string]$value)) { break }\n";
        scriptfile << "    $row++\n";
        scriptfile << "}\n";
        scriptfile << "for ($index = 0; $index -lt $simulation.Count; $index++) {\n";
        scriptfile << "    $sheet.Cells.Item($row, $index + 1).Value2 = [double]$simulation[$index]\n";
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

    return runhiddenprocessandwait(commandline);
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

bool appendclimatebenchmarkworkbook(planet& world)
{
    const vector<string> codes = orderedclimatecodes();
    const vector<long long> simulationcounts = collectsimulatedclimatecounts(world);

    return updateclimatebenchmarkworkbook(codes, simulationcounts);
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
