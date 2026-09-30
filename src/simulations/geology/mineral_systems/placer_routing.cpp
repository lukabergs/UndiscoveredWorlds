#include "placer_routing.hpp"

#include "deterministic_random.hpp"
#include "mineral_random.hpp"

#include <algorithm>
#include <cmath>
#include <map>

namespace mineralsystems
{
namespace
{
// Remaining flux below this fraction of a source settles in the current cell.
constexpr double kNegligibleFlux = 1.0e-9;
constexpr double kMinimumTransportKm = 1.0e-3;
// Sources below this share of a placer are summarised in the total only.
constexpr double kMinimumSourceShare = 0.01;
}

PlacerRouting routeplacergold(const DrainageNetwork& network, const std::vector<PlacerSource>& sources)
{
    PlacerRouting routing;
    const std::size_t count = network.cellCount();
    routing.riverGoldT.assign(count, 0.0);
    const auto trap = [&](std::int32_t cell, std::uint32_t depositid, double gold)
    {
        routing.riverGoldT[static_cast<std::size_t>(cell)] += gold;
        routing.riverT += gold;
        routing.contributions.push_back({ cell, depositid, gold });
    };

    for (const PlacerSource& source : sources)
    {
        if (source.goldT <= 0.0 || source.cell < 0 || static_cast<std::size_t>(source.cell) >= count)
            continue;
        routing.injectedT += source.goldT;
        double flux = source.goldT;
        std::int32_t cell = source.cell;

        for (std::size_t visited = 0;; ++visited)
        {
            const auto index = static_cast<std::size_t>(cell);
            if (network.sea[index] != 0)
            {
                routing.exportedT += flux;
                break;
            }
            if (network.sink[index] != 0)
            {
                routing.lakeT += flux;
                break;
            }
            const std::int32_t next = network.receiver[index];
            // Network cycles and dead ends retain what is left.
            if (next < 0 || static_cast<std::size_t>(next) >= count || visited >= count
                || flux <= kNegligibleFlux * source.goldT)
            {
                trap(cell, source.depositId, flux);
                break;
            }
            const double length = std::max(kMinimumTransportKm, static_cast<double>(network.transportLengthKm[index]));
            const double trapped = flux * (1.0 - std::exp(-static_cast<double>(network.stepLengthKm[index]) / length));
            if (trapped > 0.0)
                trap(cell, source.depositId, trapped);
            flux -= trapped;
            cell = next;
        }
    }
    return routing;
}

std::vector<Deposit> makeplacerdeposits(const PlacerRouting& routing, const DrainageNetwork& network,
    const PlacerParameters& parameters, double timeMyr, std::uint32_t firstId, std::uint64_t seed, MineralBudget& budget)
{
    std::map<std::int32_t, std::map<std::uint32_t, double>> sharesbycell;
    for (const PlacerContribution& contribution : routing.contributions)
    {
        if (routing.riverGoldT[static_cast<std::size_t>(contribution.cell)] >= parameters.minimumDepositGoldT)
            sharesbycell[contribution.cell][contribution.depositId] += contribution.goldT;
    }

    budget.placerInjectedGoldT = routing.injectedT;
    budget.placerLakeGoldT = routing.lakeT;
    budget.placerExportedGoldT = routing.exportedT;
    budget.placerDepositGoldT = 0.0;
    budget.placerDispersedGoldT = 0.0;

    std::vector<Deposit> deposits;
    std::uint32_t id = firstId;
    for (std::size_t cell = 0; cell < routing.riverGoldT.size(); ++cell)
    {
        const double gold = routing.riverGoldT[cell];
        if (gold < parameters.minimumDepositGoldT)
        {
            budget.placerDispersedGoldT += gold;
            continue;
        }
        budget.placerDepositGoldT += gold;

        detail::Random random(deterministiccombine(seed, static_cast<std::int64_t>(cell)));
        Deposit deposit;
        deposit.id = id++;
        deposit.system = DepositSystem::placer_gold;
        deposit.origin = DepositOrigin::redistributed;
        const int x = static_cast<int>(cell % static_cast<std::size_t>(network.width));
        const int y = static_cast<int>(cell / static_cast<std::size_t>(network.width));
        deposit.u = (x + random.uniform()) / network.width;
        deposit.v = (y + random.uniform()) / network.height;
        const double grade = std::max(1.0e-6, random.lognormal(parameters.gravelGradeGt.median, parameters.gravelGradeGt.sigmaLn));
        deposit.goldGradeGt = grade;
        deposit.oreMt = gold / grade;
        deposit.formationDepthM = random.uniform(parameters.overburdenM.minimum, parameters.overburdenM.maximum);
        deposit.verticalExtentM = random.uniform(parameters.payThicknessM.minimum, parameters.payThicknessM.maximum);
        deposit.formation.timeMyr = timeMyr;

        std::vector<SourceShare> shares;
        for (const auto& [depositid, tonnes] : sharesbycell[static_cast<std::int32_t>(cell)])
            shares.push_back({ depositid, tonnes });
        std::sort(shares.begin(), shares.end(), [](const SourceShare& a, const SourceShare& b)
        {
            return a.goldTonnes != b.goldTonnes ? a.goldTonnes > b.goldTonnes : a.depositId < b.depositId;
        });
        for (const SourceShare& share : shares)
        {
            if (deposit.sources.size() >= parameters.maximumSources || share.goldTonnes < kMinimumSourceShare * gold)
                break;
            deposit.sources.push_back(share);
        }
        deposits.push_back(std::move(deposit));
    }
    return deposits;
}
}
