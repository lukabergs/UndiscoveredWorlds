#pragma once

#include "climate_benchmark_outputs.hpp"
#include "climate_weather.hpp"

#include <array>
#include <filesystem>
#include <vector>

class planet;

namespace climatevalidation
{
void captureweatherstatistics(const planet& world, int season, int columns, int rows,
    const climateweather::WeatherStatistics& statistics);
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
    planet& world,
    int season,
    const std::vector<std::vector<float>>& surfaceu,
    const std::vector<std::vector<float>>& surfacev,
    const std::vector<std::vector<float>>& upperu,
    const std::vector<std::vector<float>>& upperv);
bool capturedcirculationwind(
    const planet& world,
    int season,
    int x,
    int y,
    bool upper,
    float& u,
    float& v);
bool exportcirculationdiagnostics(
    const std::filesystem::path& outputdirectory,
    planet& world,
    const climatebenchmarkmapselection& selection);
bool circulationflowvisualizationenabled(const planet& world);
int circulationflowvisualizationmaxcolumns();
int circulationflowvisualizationmaxcells();
bool exportcirculationreferencecomparisons(
    const std::filesystem::path& diagnosticoutputdirectory,
    const std::filesystem::path& referencecacheoutputdirectory,
    planet& world,
    const circulationreferencewindfields& referencefields,
    const climatebenchmarkmapselection& selection);
}
