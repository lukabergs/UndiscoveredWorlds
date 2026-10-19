#pragma once

#include "simulations/geology/mineral_systems/gold_formation.hpp"
#include "simulations/geology/mineral_systems/mineral_deposit.hpp"

#include <cstdint>
#include <vector>

// Mass-conserving transport of liberated gold along a drainage network. Each
// cell traps a fraction 1 - exp(-step / L) of the gold passing through it,
// where L is the local travel distance of detrital gold. Gold reaching the sea
// is exported; lakes and closed depressions retain all remaining gold.
namespace mineralsystems
{
struct DrainageNetwork
{
    // A width x height grid, row-major (a lat-lon world), or any cell grid (a cube) when
    // centreU and centreV are given: per cell the centre's fractions u eastward from
    // longitude -180 and v southward from the north pole, as Deposit::u and v. A placer
    // lies at a uniform point of its row-major cell, or at its cell's centre.
    int width = 0;
    int height = 0;
    std::vector<double> centreU;
    std::vector<double> centreV;
    std::vector<std::int32_t> receiver;   // Downstream cell, or -1.
    std::vector<std::uint8_t> sea;        // Outlet: gold is exported.
    std::vector<std::uint8_t> sink;       // Lake or closed depression.
    std::vector<float> stepLengthKm;      // Distance to the receiver.
    std::vector<float> transportLengthKm; // e-folding travel distance of gold.

    std::size_t cellCount() const
    {
        return centreU.empty() ? static_cast<std::size_t>(width) * static_cast<std::size_t>(height) : centreU.size();
    }
};

struct PlacerSource
{
    std::uint32_t depositId = 0;
    std::int32_t cell = -1;
    double goldT = 0.0;
};

struct PlacerContribution
{
    std::int32_t cell = -1;
    std::uint32_t depositId = 0;
    double goldT = 0.0;
};

struct PlacerRouting
{
    std::vector<double> riverGoldT; // Trapped in channels and floodplains, per cell.
    std::vector<PlacerContribution> contributions;
    double injectedT = 0.0;
    double riverT = 0.0;
    double lakeT = 0.0;
    double exportedT = 0.0;
};

PlacerRouting routeplacergold(const DrainageNetwork& network, const std::vector<PlacerSource>& sources);

struct PlacerParameters
{
    double minimumDepositGoldT = 0.5; // Smaller accumulations count as dispersed gold.
    LognormalRange gravelGradeGt{ 0.25, 0.6 };
    UniformRange overburdenM{ 0.0, 20.0 };
    UniformRange payThicknessM{ 1.0, 5.0 };
    std::size_t maximumSources = 6;
};

// Placer deposits for cells whose trapped gold reaches the threshold. Updates
// the placer entries of the budget; ids start at firstId.
std::vector<Deposit> makeplacerdeposits(const PlacerRouting& routing, const DrainageNetwork& network,
    const PlacerParameters& parameters, double timeMyr, std::uint32_t firstId, std::uint64_t seed, MineralBudget& budget);
}
