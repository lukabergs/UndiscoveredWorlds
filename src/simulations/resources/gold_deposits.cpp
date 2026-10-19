#include "simulations/resources/gold_deposits.hpp"

#include "core/deterministic_random.hpp"
#include "simulations/geology/mineral_systems/mineral_history.hpp"
#include "simulations/geology/mineral_systems/placer_routing.hpp"
#include "simulations/resources/physical_layers_internal.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace physical_layers
{
using namespace detail;
using namespace cellqueries;

namespace
{
constexpr double kPi = 3.14159265358979323846;
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

// The lat-lon (export) columns of the world: the grid's on lat-lon, 4 N on a cube.
int exportcolumns(const mesh::CellGrid& grid)
{
    return grid.kind() == mesh::GridKind::cube ? 4 * grid.facecells() : grid.columns();
}

// The world cell of a deposit's fractions (u eastward from longitude -180, v southward
// from the north pole) and its export pixel: on lat-lon the pixel's cell, on a cube the
// cell containing the point.
CellId depositcell(const mesh::CellGrid& grid, double u, double v, int& x, int& y)
{
    const int columns = exportcolumns(grid);
    x = worldcoordinate(u, columns);
    y = worldcoordinate(v, columns / 2);
    if (grid.kind() == mesh::GridKind::latlon)
        return grid.cellof(x, y);
    const double latitude = kPi * (0.5 - v);
    const double longitude = kPi * (2.0 * u - 1.0);
    return grid.locate(mesh::unitvector(latitude, longitude));
}

// The relief of the 3 x 3 box around a cell.
float localrelief(const planet& world, CellId c)
{
    float lowest = world.mapfloat(c);
    float highest = lowest;
    anyinbox(world, c, 1, [&](CellId n)
    {
        const float elevation = world.mapfloat(n);
        lowest = std::min(lowest, elevation);
        highest = std::max(highest, elevation);
        return false;
    });
    return highest - lowest;
}

void evaluatesurface(const planet& world, mineralsystems::Deposit& deposit)
{
    mineralsystems::SurfaceState& surface = deposit.surface;
    int x = 0;
    int y = 0;
    const CellId c = depositcell(world.grid(), deposit.u, deposit.v, x, y);
    surface.evaluated = true;
    surface.cell = c;
    surface.worldX = x;
    surface.worldY = y;
    surface.submerged = world.sea(c) != 0;
    surface.elevationM = static_cast<float>(world.mapfloat(c) - world.sealevel());
    surface.localReliefM = localrelief(world, c);
    surface.annualRainfallMm = 12.0f * world.averainfloat(c);
    surface.iceCovered = !surface.submerged && world.maxtempfloat(c) < 0;
}

// The erosion run's placer residence weights (erosion plan B4, review): over a run
// of D Myr the gold liberated before it decays by exp(-D / tau), and gold released
// at a constant rate during it remains at tau / D x (1 - exp(-D / tau)) of the total,
// as the tectonic exhumation's per-step decay gives (mineral_history.cpp).
struct ResidenceWeights
{
    double carried = 1.0;   // for the gold liberated before the erosion run
    double released = 1.0;  // for the gold the erosion run releases
};

ResidenceWeights residenceweights(double durationMyr)
{
    const double residence = mineralsystems::ErosionParameters{}.placerResidenceMyr;
    ResidenceWeights weights;
    if (durationMyr > 0.0 && residence > 0.0)
    {
        weights.carried = std::exp(-durationMyr / residence);
        weights.released = residence / durationMyr * (1.0 - weights.carried);
    }
    return weights;
}

// Post-erosion exhumation of a deposit by the erosion run's net lowering at its
// cell (denudation minus aggradation, never negative: a buried deposit is not
// exhumed), booked as the tectonic exhumation books it (mineral_history.cpp): the
// gold released joins the budget's eroded gold in full and the recent liberated gold
// at its residence weight, the earlier liberated gold decays over the run, and a
// fully eroded ore body takes the eroded fate at the end of the run. Applied once
// per world: layers.finalize must not run twice on an eroded world (plan section 6).
void exhumeaftererosion(mineralsystems::Deposit& deposit, double lowering_m, double timeMyr,
    const ResidenceWeights& weights, mineralsystems::MineralBudget& budget)
{
    using namespace mineralsystems;
    deposit.recentLiberatedGoldT *= weights.carried;
    if (deposit.fate != DepositFate::preserved)
        return;
    const double released = exhumedeposit(deposit, std::max(0.0, lowering_m));
    deposit.recentLiberatedGoldT += released * weights.released;
    budget.erodedGoldT += released;
    if (deposit.erodedFraction >= 1.0)
    {
        deposit.fate = DepositFate::eroded;
        deposit.fateTimeMyr = timeMyr;
        ++budget.erodedCount;
    }
}

// Present drainage on the final world, by cell id; a step is the neighbour table's distance
// to the receiver on the world's sphere (on lat-lon SphereMetric's: dy, the row's dx, or
// the hypotenuse with the two rows' mean dx for a diagonal).
mineralsystems::DrainageNetwork drainagenetwork(const planet& world)
{
    const mesh::CellGrid& grid = world.grid();
    const bool cube = grid.kind() == mesh::GridKind::cube;
    mineralsystems::DrainageNetwork network;
    std::vector<std::int8_t> slots;
    network.receiver = drainagereceivers(world, &slots);
    const size_t count = grid.cellcount();
    if (cube)
    {
        network.centreU.resize(count);
        network.centreV.resize(count);
        for (CellId c = 0; c < count; ++c)
        {
            network.centreU[c] = std::clamp((grid.longitude_rad(c) + kPi) / (2.0 * kPi), 0.0, 1.0);
            network.centreV[c] = std::clamp(0.5 - grid.latitude_rad(c) / kPi, 0.0, 1.0);
        }
    }
    else
    {
        network.width = grid.columns();
        network.height = grid.rows();
    }
    network.sea.assign(count, 0);
    network.sink.assign(count, 0);
    network.stepLengthKm.assign(count, 0.0f);
    network.transportLengthKm.assign(count, 0.0f);

    for (CellId c = 0; c < count; ++c)
    {
        const bool delta = isdeltacell(world, c);
        network.sea[c] = world.sea(c) != 0 && !delta ? 1 : 0;
        network.sink[c] = islakecell(world, c) ? 1 : 0;
        // The routing stops at a cell without a receiver, so its step is never read.
        const double step = network.receiver[c] >= 0 ? grid.distance_m(c, slots[c]) / 1000.0 : 0.0;
        network.stepLengthKm[c] = static_cast<float>(step);
        const double deposition = world.depositionpotential(c) / 100.0;
        network.transportLengthKm[c] = static_cast<float>(kChannelTransportKm
            * (1.0 - kDepositionTrapping * deposition) * (delta ? kDeltaTransportFactor : 1.0));
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
    // An eroded world (erosion_duration > 0) exhumes the deposits once more by
    // the erosion run's net rock removal (denudation minus aggradation) before the
    // liberated gold is routed. That is rock removed, not surface lowering: with the
    // isostatic rebound (landscape.isostasy.*) the surface lowers by only part of it.
    const planet& reader = world;
    const bool eroded = reader.erosionduration() > 0.0f;
    const worldgrid<float>* denudation = eroded ? &reader.denudation() : nullptr;
    const worldgrid<float>* aggradation = eroded ? &reader.aggradation() : nullptr;
    const ResidenceWeights weights = residenceweights(eroded ? static_cast<double>(reader.erosionduration()) : 0.0);
    const double erodedtimeMyr = static_cast<double>(world.tectonictimemyr()) +
        (eroded ? static_cast<double>(reader.erosionduration()) : 0.0);
    std::vector<PlacerSource> sources;
    std::uint32_t lastid = 0;
    for (Deposit& deposit : deposits)
    {
        evaluatesurface(world, deposit);
        lastid = std::max(lastid, deposit.id);
        const CellId c = deposit.surface.cell;
        if (denudation != nullptr)
            exhumeaftererosion(deposit, static_cast<double>((*denudation)[c]) - static_cast<double>((*aggradation)[c]),
                erodedtimeMyr, weights, budget);
        if (deposit.recentLiberatedGoldT > 0.0)
            sources.push_back({ deposit.id, static_cast<std::int32_t>(c), deposit.recentLiberatedGoldT });
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
