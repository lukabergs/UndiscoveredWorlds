#include "mineral_history.hpp"

#include "deterministic_random.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace mineralsystems
{
namespace
{
enum class Stream : int
{
    inherited = 1,
    formation = 2,
};

double wrapfraction(double value, std::uint32_t period)
{
    const double size = static_cast<double>(period);
    double wrapped = std::fmod(value, size);
    if (wrapped < 0.0)
        wrapped += size;
    const double fraction = wrapped / size;
    return fraction >= 1.0 ? 0.0 : fraction;
}
}

std::uint64_t mineralseed(long tectonicSeed)
{
    constexpr int kMineralSeedSalt = 0x4d494e;
    return deterministiccontextseed(tectonicSeed, kMineralSeedSalt);
}

double denudationrate(double elevationM, const ErosionParameters& parameters)
{
    if (elevationM <= 0.0 || parameters.referenceElevationM <= 0.0)
        return 0.0;
    const double rate = parameters.referenceRateMMyr * std::pow(elevationM / parameters.referenceElevationM, parameters.exponent);
    return std::min(rate, parameters.maximumRateMMyr);
}

double exhumedeposit(Deposit& deposit, double denudationM)
{
    deposit.exhumationM += std::max(0.0, denudationM);
    const double cut = deposit.exhumationM - deposit.formationDepthM;
    const double fraction = deposit.verticalExtentM > 0.0
        ? std::clamp(cut / deposit.verticalExtentM, 0.0, 1.0)
        : (cut >= 0.0 ? 1.0 : 0.0);
    if (fraction <= deposit.erodedFraction)
        return 0.0;
    const double released = deposit.containedGoldT() * (fraction - deposit.erodedFraction);
    deposit.erodedFraction = fraction;
    return released;
}

MineralHistory::MineralHistory(const MineralParameters& parameters, std::uint64_t seed)
    : parameters_(parameters), seed_(seed), tracker_(parameters.birthToleranceMyr)
{
}

bool MineralHistory::add(NewDeposit item, const TectonicFrame& frame)
{
    Marker marker;
    marker.x = item.x;
    marker.y = item.y;
    if (!tracker_.attach(marker, frame))
        return false;

    Deposit& deposit = item.deposit;
    deposit.id = nextId_++;
    const int system = static_cast<int>(deposit.system);
    if (deposit.origin == DepositOrigin::inherited)
    {
        budget_.inheritedGoldT += deposit.containedGoldT();
        ++budget_.inheritedCount;
    }
    else
    {
        budget_.formedGoldT[system] += deposit.containedGoldT();
        ++budget_.formedCount[system];
    }
    deposits_.push_back(std::move(deposit));
    markers_.push_back(marker);
    return true;
}

void MineralHistory::begin(const TectonicFrame& frame)
{
    tracker_.begin(frame);
    timeMyr_ = frame.timeMyr;
    if (!parameters_.inheritedProvinces)
        return;
    const auto inherited = sampleinheritedprovinces(frame, parameters_.formation,
        deterministiccombine(seed_, static_cast<int>(Stream::inherited)), provinces_);
    for (const NewDeposit& item : inherited)
        add(item, frame);
}

void MineralHistory::exhume(const TectonicFrame& frame, double dtMyr)
{
    const double residence = std::max(1.0e-9, parameters_.erosion.placerResidenceMyr);
    const double decay = std::exp(-dtMyr / residence);
    for (std::size_t i = 0; i < deposits_.size(); ++i)
    {
        Deposit& deposit = deposits_[i];
        deposit.recentLiberatedGoldT *= decay;
        if (deposit.fate != DepositFate::preserved)
            continue;

        const double elevation = frame.elevationM != nullptr ? frame.elevationM[markercell(markers_[i], frame)] : 0.0;
        const double denudation = denudationrate(elevation, parameters_.erosion) * dtMyr;
        const double released = exhumedeposit(deposit, denudation);
        deposit.recentLiberatedGoldT += released;
        budget_.erodedGoldT += released;
        if (deposit.erodedFraction >= 1.0)
        {
            deposit.fate = DepositFate::eroded;
            deposit.fateTimeMyr = frame.timeMyr;
            ++budget_.erodedCount;
        }
    }
}

void MineralHistory::step(const TectonicFrame& frame)
{
    last_ = {};
    if (frame.timeMyr <= timeMyr_) return;
    last_.markers = tracker_.advance(markers_, frame);
    for (std::size_t i = 0; i < deposits_.size(); ++i)
    {
        Deposit& deposit = deposits_[i];
        if (markers_[i].state != MarkerState::subducted || deposit.fate != DepositFate::preserved)
            continue;
        deposit.fate = DepositFate::subducted;
        deposit.fateTimeMyr = markers_[i].stateTimeMyr;
        budget_.subductedGoldT += deposit.remainingGoldT();
        ++budget_.subductedCount;
    }

    // Re-reading a frame does not advance geological time.
    const double dt = std::max(0.0, frame.timeMyr - timeMyr_);
    timeMyr_ = std::max(timeMyr_, frame.timeMyr);
    last_.dtMyr = dt;
    ++updates_;
    if (dt <= 0.0)
        return;

    exhume(frame, dt);
    FormationSample sample = sampleformation(frame, parameters_.formation, dt,
        deterministiccombine(seed_, static_cast<int>(Stream::formation), static_cast<std::int64_t>(updates_)));
    last_.orogenicAreaKm2 = sample.orogenicAreaKm2;
    last_.arcAreaKm2 = sample.arcAreaKm2;
    for (NewDeposit& item : sample.deposits)
        last_.formed += add(std::move(item), frame) ? 1 : 0;
}

void MineralHistory::finish(const TectonicFrame& frame)
{
    for (std::size_t i = 0; i < deposits_.size(); ++i)
    {
        deposits_[i].u = wrapfraction(markers_[i].x, frame.width);
        deposits_[i].v = wrapfraction(markers_[i].y, frame.height);
    }
}
}
