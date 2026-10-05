#pragma once

#include "io/climate_benchmark_outputs.hpp"
#include "io/raster_image.hpp"
#include "simulations/climate/climate_weather.hpp"
#include "simulations/climate/climate_hydrology.hpp"
#include "simulations/climate/climate_ocean_dynamics.hpp"

#include <array>
#include <filesystem>
#include <vector>

class climateworld;

namespace climatevalidation
{
struct circulationreferencewindfields
{
    int columns = 0;
    int rows = 0;
    std::array<std::vector<float>, 4> surfaceu;
    std::array<std::vector<float>, 4> surfacev;
    std::array<std::vector<float>, 4> upperu;
    std::array<std::vector<float>, 4> upperv;
};

void capturecirculationwindfields(
    climateworld& world,
    int season,
    const std::vector<std::vector<float>>& surfaceu,
    const std::vector<std::vector<float>>& surfacev,
    const std::vector<std::vector<float>>& upperu,
    const std::vector<std::vector<float>>& upperv);
bool capturedcirculationwind(
    const climateworld& world,
    int season,
    int x,
    int y,
    bool upper,
    float& u,
    float& v);
// Every exported PNG map has the world grid's size (width x width/2).
// Coarser rasters are enlarged without blending, so each pixel keeps its
// cell's colour.
rasterimage::Image expandclimatemapimage(const rasterimage::Image& cells, int mapcolumns);
bool saveclimatemapimage(const rasterimage::Image& cells, const std::filesystem::path& path, int mapcolumns);
bool exportcirculationdiagnostics(
    const std::filesystem::path& outputdirectory,
    climateworld& world,
    const climatebenchmarkmapselection& selection);
bool circulationflowvisualizationenabled(const climateworld& world);
int circulationflowvisualizationmaxcolumns();
int circulationflowvisualizationmaxcells();
bool exportcirculationreferencecomparisons(
    const std::filesystem::path& diagnosticoutputdirectory,
    const std::filesystem::path& referencecacheoutputdirectory,
    climateworld& world,
    const circulationreferencewindfields& referencefields,
    const climatebenchmarkmapselection& selection);
}
