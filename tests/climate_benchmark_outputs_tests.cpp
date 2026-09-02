#include "climate_benchmark_outputs.hpp"

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

// Keep behavior checks active in Release, too.
#define assert(condition) do { if (!(condition)) { std::cerr << "FAIL line " << __LINE__ << ": " #condition << '\n'; std::exit(1); } } while (false)

int main()
{
    const climatebenchmarkmapselection defaults =
        climatebenchmarkmapselection::defaultselection();
    assert(defaults.requests().size() == 5);
    assert(defaults.includes(climatebenchmarkmapkind::koppen));
    assert(defaults.includes(climatebenchmarkmapkind::surfacewindlic, 0));
    assert(defaults.includes(climatebenchmarkmapkind::surfacewindparticles, 0));
    assert(defaults.includes(climatebenchmarkmapkind::precipitation));
    assert(defaults.includes(climatebenchmarkmapkind::precipitationtiff));
    assert(!defaults.includes(climatebenchmarkmapkind::surfacewindlic, 2));

    climatebenchmarkmapselection selected;
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
    assert(selected.requests().size() == 92);
    assert(selected.includes(climatebenchmarkmapkind::columnmoistureflux, 0));
    assert(selected.includes(climatebenchmarkmapkind::surfacewindconsistency, 0));
    assert(selected.includes(climatebenchmarkmapkind::era5upperwindparticles, 3));

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
