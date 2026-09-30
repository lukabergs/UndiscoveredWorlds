#pragma once
#include <vector>
class planet;

// Overrides of configs/tectonics-baseline.json; zero selects its value.
struct WorldTectonicsOptions {
    double durationMyr = 0;
    int plateCount = 0;
    double deltaTimeMyr = 0;
    bool mineralHistory = true;
};
void generateworldtectonics(planet& world, std::vector<std::vector<bool>>& shelves,
                            const WorldTectonicsOptions& options = {});
