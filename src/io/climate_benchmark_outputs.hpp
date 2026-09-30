#pragma once

#include <string>
#include <vector>

enum class climatebenchmarkmapkind
{
    koppen,
    temperature,
    precipitation,
    precipitationtiff,
    surfacewindlic,
    surfacewindparticles,
    upperwindlic,
    upperwindparticles,
    surfacewindvectorerror,
    upperwindvectorerror,
    surfacewindspeed,
    upperwindspeed,
    surfacedivergence,
    moisturefluxconvergence,
    boundarymoistureflux,
    freemoistureflux,
    columnmoistureflux,
    verticalascent,
    surfacewindconsistency,
    upperwindconsistency,
    seasurfacetemperature,
    oceancurrents,
    columnheating,
    era5surfacewindlic,
    era5surfacewindparticles,
    era5upperwindlic,
    era5upperwindparticles
};

struct climatebenchmarkmaprequest
{
    climatebenchmarkmapkind kind = climatebenchmarkmapkind::koppen;
    int season = -1;
};

class climatebenchmarkmapselection
{
public:
    static climatebenchmarkmapselection defaultselection();
    static climatebenchmarkmapselection allselection();

    void clear();
    void add(climatebenchmarkmapkind kind, int season = -1);
    bool includes(climatebenchmarkmapkind kind, int season = -1) const;
    bool requirescirculation() const;
    bool requiresreferencewinds() const;
    // Rendering only: does not change world or simulation grids.
    int windmapcolumns() const { return windmapcolumns_; }
    void setwindmapcolumns(int columns);
    const std::vector<climatebenchmarkmaprequest>& requests() const;

private:
    std::vector<climatebenchmarkmaprequest> requests_;
    int windmapcolumns_ = 2048;
};

bool addclimatebenchmarkmapargument(
    climatebenchmarkmapselection& selection,
    const std::string& argument,
    std::string* failuremessage = nullptr);
const char* climatebenchmarkseasonname(int season);
std::string climatebenchmarkmapid(const climatebenchmarkmaprequest& request);
std::string climatebenchmarkmapfilename(
    const climatebenchmarkmaprequest& request,
    int runid);
// Relative to the maps root (fields root for the precipitation GeoTIFF).
std::string climatebenchmarkmaprelativepath(const climatebenchmarkmaprequest& request, int runid);
bool climatebenchmarkmapisreference(climatebenchmarkmapkind kind);
double climatebenchmarkweightedrelativeerror(
    const std::vector<long long>& simulationcounts,
    const std::vector<double>& referencecounts,
    int horizontalresolution,
    int referencehorizontalresolution = 3600);
