#include "climate_benchmark_outputs.hpp"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <stdexcept>

// Keep behavior checks active in Release, too.
#define assert(condition) do { if (!(condition)) { std::cerr << "FAIL line " << __LINE__ << ": " #condition << '\n'; std::exit(1); } } while (false)

int main()
{
    const climatebenchmarkmapselection defaults =
        climatebenchmarkmapselection::defaultselection();
    assert(defaults.requests().size() == 20);
    assert(defaults.windmapcolumns() == 2048);
    assert(defaults.includes(climatebenchmarkmapkind::oceancurrents, 0));
    assert(defaults.includes(climatebenchmarkmapkind::oceancurrents, 2));
    assert(defaults.includes(climatebenchmarkmapkind::koppen));
    assert(defaults.includes(climatebenchmarkmapkind::surfacewindlic, 0));
    assert(defaults.includes(climatebenchmarkmapkind::surfacewindparticles, 0));
    assert(defaults.includes(climatebenchmarkmapkind::precipitation));
    assert(defaults.includes(climatebenchmarkmapkind::precipitationtiff));
    assert(defaults.includes(climatebenchmarkmapkind::surfacewindlic, 2));
    assert(defaults.includes(climatebenchmarkmapkind::surfacewindparticles, 2));
    assert(defaults.includes(climatebenchmarkmapkind::temperature));
    assert(defaults.includes(climatebenchmarkmapkind::seasurfacetemperature, 0));
    assert(defaults.includes(climatebenchmarkmapkind::seasurfacetemperature, 2));
    for (const int season : { 0, 2 })
    {
        assert(defaults.includes(climatebenchmarkmapkind::surfacewindspeed, season));
        assert(defaults.includes(climatebenchmarkmapkind::upperwindspeed, season));
        assert(defaults.includes(climatebenchmarkmapkind::upperwindlic, season));
        assert(defaults.includes(climatebenchmarkmapkind::upperwindparticles, season));
    }

    climatebenchmarkmapselection selected;
    selected.setwindmapcolumns(512);
    selected.clear();
    assert(selected.windmapcolumns() == 512);
    for (const int width : { 0, 63, 65, 2049, 4096 })
    {
        bool rejected = false;
        try { selected.setwindmapcolumns(width); } catch (const std::invalid_argument&) { rejected = true; }
        assert(rejected && selected.windmapcolumns() == 512);
    }
    std::string failure;
    assert(addclimatebenchmarkmapargument(
        selected,
        "koppen,oct-s-wind-lic,jan-era5-s-wind-part,precip-tif",
        &failure));
    assert(selected.requests().size() == 4);
    assert(selected.includes(climatebenchmarkmapkind::surfacewindlic, 3));
    assert(selected.includes(climatebenchmarkmapkind::era5surfacewindparticles, 0));
    assert(selected.requiresreferencewinds());
    assert(climatebenchmarkmapfilename(
        { climatebenchmarkmapkind::surfacewindparticles, 0 }, 138) ==
        "138_jan_s_wind_part.png");
    assert(climatebenchmarkmapfilename(
        { climatebenchmarkmapkind::precipitationtiff, -1 }, 138) ==
        "138_precip.tif");

    failure.clear();
    assert(!addclimatebenchmarkmapargument(selected, "january-s-wind-lic", &failure));
    assert(!failure.empty());

    selected.clear();
    assert(addclimatebenchmarkmapargument(selected, "all", &failure));
    assert(selected.requests().size() == 96);
    assert(selected.includes(climatebenchmarkmapkind::upperwindspeed, 2));
    assert(climatebenchmarkmaprelativepath({ climatebenchmarkmapkind::upperwindspeed, 2 }, 273) ==
        "climate/wind/speed/3/u/273.png");
    assert(selected.includes(climatebenchmarkmapkind::columnmoistureflux, 0));
    assert(selected.includes(climatebenchmarkmapkind::surfacewindconsistency, 0));
    assert(selected.includes(climatebenchmarkmapkind::era5upperwindparticles, 3));
    assert(climatebenchmarkmaprelativepath(
        { climatebenchmarkmapkind::surfacewindlic, 2 }, 160) ==
        "climate/wind/lic/3/s/160.png");
    assert(climatebenchmarkmaprelativepath(
        { climatebenchmarkmapkind::precipitationtiff, -1 }, 160) ==
        "climate/rain/annual/s/160.tif");
    std::set<std::string> paths;
    for (const auto& request : selected.requests())
    {
        if (climatebenchmarkmapisreference(request.kind))
            continue;
        const auto path = climatebenchmarkmaprelativepath(request, 160);
        assert(path.rfind("climate/", 0) == 0);
        assert(path.substr(path.size() - 8) == "/160.png" || path.substr(path.size() - 8) == "/160.tif");
        assert(path.find("/0/") == std::string::npos);
        assert(paths.insert(path).second);
    }

    const std::vector<double> referencecounts = { 16.0, 32.0, 48.0 };
    const std::vector<long long> exactquartercounts = { 1, 2, 3 };
    assert(climatebenchmarkweightedrelativeerror(
        exactquartercounts, referencecounts, 512, 2048) == 0.0);
    const std::vector<long long> biasedcounts = { 2, 2, 2 };
    assert(climatebenchmarkweightedrelativeerror(
        biasedcounts, referencecounts, 512, 2048) == 1.0 / 3.0);

    const std::vector<long long> run145counts = {
        1939, 392, 1198, 93, 2175, 1037, 1535, 1135, 1539, 851, 23,
        324, 74, 1, 1594, 1475, 100, 54, 80, 2691, 0, 0, 2, 982, 0,
        278, 439, 6271, 0, 3598, 14862
    };
    const std::vector<double> workbookreferencecounts = {
        62064.0, 42738.0, 75295.5, 75295.5, 200752.0, 75185.0,
        74499.0, 90464.0, 15833.0, 9967.0, 18.0, 36456.0, 13728.0,
        20.0, 59648.0, 37240.0, 51.0, 2341.0, 6761.0, 28161.0,
        200.0, 11837.0, 19978.0, 37404.0, 4615.0, 19815.0,
        124260.0, 266833.0, 2338.0, 149192.0, 666267.0
    };
    assert(std::abs(climatebenchmarkweightedrelativeerror(
        run145counts, workbookreferencecounts, 512) -
        0.40125772177493896) < 1e-12);
    return 0;
}
