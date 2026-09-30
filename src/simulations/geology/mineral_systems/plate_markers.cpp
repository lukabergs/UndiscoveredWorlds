#include "plate_markers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace mineralsystems
{
namespace
{
// Snapped markers stay this far inside the target cell, in cells.
constexpr double kSnapMargin = 0.25;

double wrapcoordinate(double value, std::uint32_t period)
{
    const double size = static_cast<double>(period);
    double wrapped = std::fmod(value, size);
    if (wrapped < 0.0)
        wrapped += size;
    return wrapped >= size ? 0.0 : wrapped;
}

double wrappedoffset(double point, double centre, std::uint32_t period)
{
    const double size = static_cast<double>(period);
    double delta = std::fmod(point - centre, size);
    if (delta > 0.5 * size)
        delta -= size;
    else if (delta < -0.5 * size)
        delta += size;
    return delta;
}

std::vector<std::pair<std::uint64_t, std::uint32_t>> uidindex(const TectonicFrame& frame)
{
    std::vector<std::pair<std::uint64_t, std::uint32_t>> index;
    index.reserve(frame.plates.size());
    for (std::uint32_t plate = 0; plate < frame.plates.size(); ++plate)
        index.emplace_back(frame.plates[plate].uid, plate);
    std::sort(index.begin(), index.end());
    return index;
}

std::uint32_t findplate(const std::vector<std::pair<std::uint64_t, std::uint32_t>>& index, std::uint64_t uid)
{
    const auto found = std::lower_bound(index.begin(), index.end(), std::make_pair(uid, std::uint32_t{0}));
    return found != index.end() && found->first == uid ? found->second : kNoPlate;
}

bool continentalcell(const TectonicFrame& frame, std::size_t cell)
{
    const std::uint8_t crust = frame.crustClass[cell];
    return crust == crustcode::continental || crust == crustcode::transitional;
}

// Moves a continental-hosted marker to the nearest neighbouring cell where its
// host owns continental or transitional crust.
bool snaptohost(Marker& marker, const TectonicFrame& frame, std::uint32_t host)
{
    const double cellx = std::floor(marker.x);
    const double celly = std::floor(marker.y);
    double bestdistance = std::numeric_limits<double>::max();
    double bestx = 0.0;
    double besty = 0.0;

    for (int dy = -1; dy <= 1; ++dy)
    {
        for (int dx = -1; dx <= 1; ++dx)
        {
            if (dx == 0 && dy == 0)
                continue;
            const double left = cellx + dx;
            const double top = celly + dy;
            const auto nx = static_cast<std::uint32_t>(wrapcoordinate(left, frame.width));
            const auto ny = static_cast<std::uint32_t>(wrapcoordinate(top, frame.height));
            const std::size_t cell = frame.cellIndex(nx, ny);
            if (frame.plateIndex[cell] != host || !continentalcell(frame, cell))
                continue;
            const double px = std::clamp(marker.x, left + kSnapMargin, left + 1.0 - kSnapMargin);
            const double py = std::clamp(marker.y, top + kSnapMargin, top + 1.0 - kSnapMargin);
            const double distance = (px - marker.x) * (px - marker.x) + (py - marker.y) * (py - marker.y);
            if (distance < bestdistance)
            {
                bestdistance = distance;
                bestx = px;
                besty = py;
            }
        }
    }

    if (bestdistance == std::numeric_limits<double>::max())
        return false;
    marker.x = wrapcoordinate(bestx, frame.width);
    marker.y = wrapcoordinate(besty, frame.height);
    return true;
}

std::vector<std::uint64_t> sorteduids(const TectonicFrame& frame)
{
    std::vector<std::uint64_t> uids;
    uids.reserve(frame.plates.size());
    for (const PlateIdentity& plate : frame.plates)
        uids.push_back(plate.uid);
    std::sort(uids.begin(), uids.end());
    return uids;
}
}

std::size_t markercell(const Marker& marker, const TectonicFrame& frame)
{
    const auto x = std::min(frame.width - 1, static_cast<std::uint32_t>(wrapcoordinate(marker.x, frame.width)));
    const auto y = std::min(frame.height - 1, static_cast<std::uint32_t>(wrapcoordinate(marker.y, frame.height)));
    return frame.cellIndex(x, y);
}

void applymaterialshift(Marker& marker, const TectonicFrame& frame)
{
    // Approximate point tracers of the column displacement. Conservative
    // subcell mixing and root/sediment redistribution are not marker transport.
    constexpr double kMaximumSubstepCells = 0.25;
    const CellShift start = frame.materialShift[markercell(marker, frame)];
    const double largest = std::max(std::abs(start.x), std::abs(start.y));
    const int steps = std::max(1, static_cast<int>(std::ceil(largest / kMaximumSubstepCells)));
    for (int step = 0; step < steps; ++step)
    {
        const CellShift shift = frame.materialShift[markercell(marker, frame)];
        marker.x = wrapcoordinate(marker.x + shift.x / steps, frame.width);
        marker.y = wrapcoordinate(marker.y + shift.y / steps, frame.height);
    }
}

void MarkerTracker::begin(const TectonicFrame& frame)
{
    knownUids_ = sorteduids(frame);
}

bool MarkerTracker::attach(Marker& marker, const TectonicFrame& frame) const
{
    const std::size_t cell = markercell(marker, frame);
    const std::uint32_t owner = frame.plateIndex[cell];
    if (owner >= frame.plates.size())
        return false;
    marker.hostUid = frame.plates[owner].uid;
    marker.continentalHost = frame.crustClass == nullptr || continentalcell(frame, cell);
    marker.hostCrustBirthMyr = frame.crustAgeMyr != nullptr ? frame.timeMyr - frame.crustAgeMyr[cell] : frame.timeMyr;
    marker.state = MarkerState::attached;
    return true;
}

MarkerStepSummary MarkerTracker::advance(std::vector<Marker>& markers, const TectonicFrame& frame)
{
    MarkerStepSummary summary;
    const auto index = uidindex(frame);
    const auto platecount = static_cast<std::uint32_t>(frame.plates.size());

    for (Marker& marker : markers)
    {
        if (marker.state != MarkerState::attached)
            continue;

        const std::uint32_t host = findplate(index, marker.hostUid);
        if (frame.materialShift != nullptr)
        {
            applymaterialshift(marker, frame);
            ++summary.moved;
        }
        else
            ++summary.orphaned;

        const std::size_t cell = markercell(marker, frame);
        const std::uint32_t owner = frame.plateIndex[cell];
        if (owner >= platecount)
            continue; // Uncovered gap: keep the host.
        const std::uint64_t owneruid = frame.plates[owner].uid;

        if (marker.continentalHost)
        {
            // Keep approximate tracers on nearby host crust after cell mixing.
            if (frame.crustClass != nullptr && !continentalcell(frame, cell) && host < platecount
                && snaptohost(marker, frame, host))
            {
                ++summary.snapped;
                continue;
            }
            if (owneruid != marker.hostUid)
            {
                marker.hostUid = owneruid;
                ++summary.rehomed;
            }
            continue;
        }

        if (owneruid == marker.hostUid)
            continue;
        // Crust overridden by an existing plate is consumed unless the owner
        // carries crust of the same age, i.e. the host crust was relabelled.
        const bool ownerexisted = std::binary_search(knownUids_.begin(), knownUids_.end(), owneruid);
        const double birth = frame.crustAgeMyr != nullptr ? frame.timeMyr - frame.crustAgeMyr[cell] : marker.hostCrustBirthMyr;
        if (ownerexisted && std::abs(birth - marker.hostCrustBirthMyr) > birthToleranceMyr_)
        {
            marker.state = MarkerState::subducted;
            marker.stateTimeMyr = frame.timeMyr;
            ++summary.subducted;
            continue;
        }
        marker.hostUid = owneruid;
        ++summary.rehomed;
    }

    knownUids_ = sorteduids(frame);
    return summary;
}
}
