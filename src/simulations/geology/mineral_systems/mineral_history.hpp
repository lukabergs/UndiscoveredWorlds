#pragma once

#include "gold_formation.hpp"
#include "mineral_deposit.hpp"
#include "plate_markers.hpp"
#include "tectonic_frame.hpp"

#include <cstdint>
#include <utility>
#include <vector>

// Records mineral deposits through a tectonic run. Deposits form from the
// setting of each update, travel with their host crust, lose ore to erosion as
// they are exhumed, and are removed with subducted crust.
namespace mineralsystems
{
// Subaerial denudation from native elevation under a reference climate. Erosion
// rises steeply with relief, and relief with elevation, so the rate follows a
// power of elevation: about 1-5 m/Myr on cratons, tens of m/Myr on uplands and
// hundreds in high orogens. It is an explicit proxy until the terrain model
// supplies per-cell erosion; submerged crust does not erode.
struct ErosionParameters
{
    double referenceRateMMyr = 5.0;
    double referenceElevationM = 1000.0;
    double exponent = 3.0;
    double maximumRateMMyr = 1000.0;
    // Liberated gold remains in the fluvial system for about this long, so a
    // placer inventory reflects erosion over the last few residence times.
    double placerResidenceMyr = 5.0;
};

// History stream for a tectonic seed; the app and world_pipeline share it,
// so the same seed and recipe give the same deposits.
std::uint64_t mineralseed(long tectonicSeed);

double denudationrate(double elevationM, const ErosionParameters& parameters);
// Removes rock above and then within the ore body. Returns the gold released;
// the ore body is eroded uniformly through its vertical extent.
double exhumedeposit(Deposit& deposit, double denudationM);

struct MineralParameters
{
    FormationParameters formation;
    ErosionParameters erosion;
    double birthToleranceMyr = 3.0;
    bool inheritedProvinces = true;
};

struct HistoryStepDiagnostics
{
    MarkerStepSummary markers;
    double dtMyr = 0.0;
    double orogenicAreaKm2 = 0.0;
    double arcAreaKm2 = 0.0;
    int formed = 0;
};

class MineralHistory
{
public:
    MineralHistory(const MineralParameters& parameters, std::uint64_t seed);

    // Seeds inherited provinces in the initial crust.
    void begin(const TectonicFrame& frame);
    // Call after every physical tectonic update.
    void step(const TectonicFrame& frame);
    // Writes fractional positions from the markers.
    void finish(const TectonicFrame& frame);

    const std::vector<Deposit>& deposits() const { return deposits_; }
    std::vector<Deposit> takedeposits() { return std::move(deposits_); }
    const std::vector<Marker>& markers() const { return markers_; }
    const MineralBudget& budget() const { return budget_; }
    const HistoryStepDiagnostics& laststep() const { return last_; }
    std::uint32_t provincecount() const { return provinces_; }

private:
    bool add(NewDeposit item, const TectonicFrame& frame);
    void exhume(const TectonicFrame& frame, double dtMyr);

    MineralParameters parameters_;
    std::uint64_t seed_;
    MarkerTracker tracker_;
    std::vector<Deposit> deposits_;
    std::vector<Marker> markers_; // Parallel to deposits_.
    MineralBudget budget_;
    HistoryStepDiagnostics last_;
    double timeMyr_ = 0.0;
    std::uint64_t updates_ = 0;
    std::uint32_t nextId_ = 1;
    std::uint32_t provinces_ = 0;
};
}
