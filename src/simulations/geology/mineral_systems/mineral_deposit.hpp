#pragma once

#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

#include "core/mesh/cell_grid.hpp"

// Discrete mineral deposits with formation provenance. Deposits are classified
// by the process that formed them; commodities are attributes of a deposit.
// Quantities are geological endowment assumptions (ore tonnage and grade), not
// reserves: discovery, extraction technology and value belong to society.
namespace mineralsystems
{
enum class DepositSystem : std::uint8_t
{
    orogenic_gold = 0,        // Metamorphic-fluid lode gold in deformed crust.
    epithermal_gold = 1,      // Shallow gold-silver veins above arc magmatism.
    porphyry_copper_gold = 2, // Arc intrusion-hosted copper with by-product gold.
    placer_gold = 3,          // Detrital gold concentrated by rivers.
};
inline constexpr int kDepositSystemCount = 4;

enum class DepositOrigin : std::uint8_t
{
    simulated = 0,     // Formed during the tectonic run from recorded conditions.
    inherited = 1,     // Seeded at initialization for pre-simulation history.
    redistributed = 2, // Concentrated from liberated gold by surface transport.
};

enum class DepositFate : std::uint8_t
{
    preserved = 0, // Ore remains in place, possibly partly eroded.
    eroded = 1,    // Ore body completely removed by exhumation.
    subducted = 2, // Host crust consumed at a subduction zone.
};

// Depth classes for display; the depth itself remains the numerical output.
enum class Exposure : std::uint8_t
{
    outcropping = 0, // Ore body intersects the surface.
    shallow = 1,     // Top within 100 m of the surface.
    buried = 2,      // Top within 1 km.
    deep = 3,
};

// Conditions recorded when a deposit formed. Codes use the tectonic contract
// enums (tectonics::contract::GeologicRegime, CrustClass, BoundaryType).
struct FormationRecord
{
    double timeMyr = 0.0; // Simulation time; negative for inherited deposits.
    std::uint8_t geologicRegime = 0;
    std::uint8_t boundaryType = 0;
    std::uint8_t hostCrustClass = 0;
    std::uint32_t boundarySegmentId = 0;
    std::uint32_t provinceId = 0; // Inherited province; 0 for simulated deposits.
    float boundaryDistanceKm = 0.0f;
    float boundaryAgeMyr = 0.0f;
    float settingIntensity = 0.0f; // 0..1 rate multiplier that drove formation.
};

struct SourceShare
{
    std::uint32_t depositId = 0;
    double goldTonnes = 0.0;
};

// Present-day surface conditions at the deposit, evaluated on the final world.
struct SurfaceState
{
    bool evaluated = false;
    // The world cell (docs/architecture/sphere-mesh.md section 4) and its lat-lon (export)
    // pixel: the cell's column and row on a lat-lon world, the deposit's pixel of the
    // 4 N x 2 N export on a cube.
    mesh::CellId cell = mesh::nocell;
    int worldX = -1;
    int worldY = -1;
    bool submerged = false;
    bool iceCovered = false;
    float elevationM = 0.0f;
    float localReliefM = 0.0f;
    float annualRainfallMm = 0.0f;
};

struct Deposit
{
    std::uint32_t id = 0;
    DepositSystem system = DepositSystem::orogenic_gold;
    DepositOrigin origin = DepositOrigin::simulated;
    DepositFate fate = DepositFate::preserved;
    double fateTimeMyr = 0.0;

    // Fractions of the native tectonic grid: u eastward, v southward, in [0, 1).
    double u = 0.0;
    double v = 0.0;

    FormationRecord formation;

    double oreMt = 0.0;          // Ore (or pay gravel) tonnage, million tonnes.
    double goldGradeGt = 0.0;    // g/t; Mt * g/t = tonnes of gold.
    double copperGradePct = 0.0; // Percent copper by mass.
    double formationDepthM = 0.0; // Depth of the ore-body top when it formed.
    double verticalExtentM = 0.0;
    double exhumationM = 0.0;     // Rock removed above the deposit since formation.
    double erodedFraction = 0.0;  // Fraction of the ore body removed, 0..1.
    // Gold released by erosion, discounted by the placer residence time: the
    // part still expected in the surface sediment system at the end of the run.
    double recentLiberatedGoldT = 0.0;

    std::vector<SourceShare> sources; // Placers: contributing primary deposits.
    SurfaceState surface;

    double containedGoldT() const { return oreMt * goldGradeGt; }
    double remainingGoldT() const { return containedGoldT() * (1.0 - erodedFraction); }
    double containedCopperT() const { return oreMt * 1.0e6 * copperGradePct / 100.0; }
    double remainingCopperT() const { return containedCopperT() * (1.0 - erodedFraction); }
    // Negative once erosion has cut into the ore body.
    double topDepthM() const { return formationDepthM - exhumationM; }
    bool isPrimary() const { return system != DepositSystem::placer_gold; }
};

// Totals in tonnes of gold, by deposit system where indexed.
struct MineralBudget
{
    double formedGoldT[kDepositSystemCount] = {};
    double inheritedGoldT = 0.0;
    double erodedGoldT = 0.0;    // Removed from ore bodies during the tectonic run and, in an eroded world, the erosion run.
    double subductedGoldT = 0.0; // Remaining gold lost with consumed host crust.
    int formedCount[kDepositSystemCount] = {};
    int inheritedCount = 0;
    int erodedCount = 0;
    int subductedCount = 0;

    // Present-day placer routing on the final drainage network.
    double placerInjectedGoldT = 0.0;
    double placerDepositGoldT = 0.0;   // Trapped in cells that form placer deposits.
    double placerDispersedGoldT = 0.0; // Trapped below the deposit threshold.
    double placerLakeGoldT = 0.0;      // Trapped in lakes and closed depressions.
    double placerExportedGoldT = 0.0;  // Delivered to the sea.
};

const char* depositsystemlabel(DepositSystem system);
const char* depositoriginlabel(DepositOrigin origin);
const char* depositfatelabel(DepositFate fate);
const char* exposurelabel(Exposure exposure);
Exposure exposureclass(double topDepthM);

// Multi-line explanation of why the deposit exists and what it contains.
// Placer sources are named when the full deposit list is supplied.
std::string describedeposit(const Deposit& deposit, const std::vector<Deposit>* deposits = nullptr);

// One row per deposit: enum codes, positions as world fractions, provenance,
// quantities and depths. Placer sources are omitted. cells adds the world cell of
// each deposit (a cube bundle; sphere-mesh.md section 4), empty where it has none.
void writedepositscsv(std::ostream& out, const std::vector<Deposit>& deposits, bool cells = false);
std::string mineralbudgetjson(const MineralBudget& budget);
}
