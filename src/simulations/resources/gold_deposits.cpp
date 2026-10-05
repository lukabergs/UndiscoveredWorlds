#include "simulations/resources/gold_deposits.hpp"

#include "core/deterministic_random.hpp"
#include "simulations/geology/mineral_systems/placer_routing.hpp"
#include "simulations/resources/physical_layers_internal.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace physical_layers
{
using namespace detail;

namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kCircumferenceKm = 2.0 * kPi * 6371.0; // As assumed by the climate solvers.
// e-folding travel distance of detrital gold in an ordinary channel; gold is
// dense, so most settles within tens of km of its source.
constexpr double kChannelTransportKm = 25.0;
constexpr double kDepositionTrapping = 0.75; // Share of travel lost on fully depositional reaches.
constexpr double kDeltaTransportFactor = 0.5;
constexpr int kPlacerSeedSalt = 0x504c43;

int worldcoordinate(double fraction, int cells)
{
    return std::clamp(static_cast<int>(std::floor(fraction * cells)), 0, cells - 1);
}

float localrelief(const planet& world, int x, int y)
{
    int lowest = world.map(x, y);
    int highest = lowest;
    for (int dy = -1; dy <= 1; ++dy)
    {
        const int ny = y + dy;
        if (ny < 0 || ny > world.height())
            continue;
        for (int dx = -1; dx <= 1; ++dx)
        {
            const int elevation = world.map(wrap(x + dx, world.width()), ny);
            lowest = std::min(lowest, elevation);
            highest = std::max(highest, elevation);
        }
    }
    return static_cast<float>(highest - lowest);
}

void evaluatesurface(const planet& world, mineralsystems::Deposit& deposit)
{
    mineralsystems::SurfaceState& surface = deposit.surface;
    const int x = worldcoordinate(deposit.u, world.width() + 1);
    const int y = worldcoordinate(deposit.v, world.height() + 1);
    surface.evaluated = true;
    surface.worldX = x;
    surface.worldY = y;
    surface.submerged = world.sea(x, y) != 0;
    surface.elevationM = static_cast<float>(world.map(x, y) - world.sealevel());
    surface.localReliefM = localrelief(world, x, y);
    surface.annualRainfallMm = 12.0f * world.averainfloat(x, y);
    surface.iceCovered = !surface.submerged && world.maxtemp(x, y) < 0;
}

// Present drainage on the final world; distances treat rows as latitudes.
mineralsystems::DrainageNetwork drainagenetwork(const planet& world)
{
    mineralsystems::DrainageNetwork network;
    network.width = world.width() + 1;
    network.height = world.height() + 1;
    network.receiver = drainagereceivers(world);
    const size_t count = network.cellCount();
    network.sea.assign(count, 0);
    network.sink.assign(count, 0);
    network.stepLengthKm.assign(count, 0.0f);
    network.transportLengthKm.assign(count, 0.0f);

    const double cellkm = kCircumferenceKm / network.width;
    for (int y = 0; y < network.height; ++y)
    {
        const double latitude = kPi * (0.5 - (y + 0.5) / network.height);
        const double eastwest = cellkm * std::cos(latitude);
        for (int x = 0; x < network.width; ++x)
        {
            const size_t index = size_t(y) * size_t(network.width) + size_t(x);
            const bool delta = isdeltacell(world, x, y);
            network.sea[index] = world.sea(x, y) != 0 && !delta ? 1 : 0;
            network.sink[index] = islakecell(world, x, y) ? 1 : 0;
            double step = cellkm;
            const int32_t next = network.receiver[index];
            if (next >= 0)
            {
                int dx = static_cast<int>(next % network.width) - x;
                if (dx > network.width / 2) dx -= network.width;
                if (dx < -network.width / 2) dx += network.width;
                const int dy = static_cast<int>(next / network.width) - y;
                step = std::hypot(dx * eastwest, dy * cellkm);
            }
            network.stepLengthKm[index] = static_cast<float>(step);
            const double deposition = world.depositionpotential(x, y) / 100.0;
            network.transportLengthKm[index] = static_cast<float>(kChannelTransportKm
                * (1.0 - kDepositionTrapping * deposition) * (delta ? kDeltaTransportFactor : 1.0));
        }
    }
    return network;
}
}

void generate_gold_deposits(planet& world)
{
    using namespace mineralsystems;
    std::vector<Deposit>& deposits = world.mineraldeposits();
    deposits.erase(std::remove_if(deposits.begin(), deposits.end(),
        [](const Deposit& deposit) { return deposit.system == DepositSystem::placer_gold; }), deposits.end());

    MineralBudget budget = world.mineralbudget();
    std::vector<PlacerSource> sources;
    std::uint32_t lastid = 0;
    for (Deposit& deposit : deposits)
    {
        evaluatesurface(world, deposit);
        lastid = std::max(lastid, deposit.id);
        if (deposit.recentLiberatedGoldT > 0.0)
            sources.push_back({ deposit.id, deposit.surface.worldY * (world.width() + 1) + deposit.surface.worldX,
                deposit.recentLiberatedGoldT });
    }

    const DrainageNetwork network = drainagenetwork(world);
    std::vector<Deposit> placers = makeplacerdeposits(routeplacergold(network, sources), network, PlacerParameters{},
        world.tectonictimemyr(), lastid + 1, deterministiccontextseed(world.seed(), kPlacerSeedSalt), budget);
    for (Deposit& placer : placers)
    {
        evaluatesurface(world, placer);
        deposits.push_back(std::move(placer));
    }
    world.setmineralbudget(budget);
}
}
