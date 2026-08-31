#pragma once

#include <filesystem>
#include <vector>

class planet;

namespace climatevalidation
{
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
    planet& world);
}
