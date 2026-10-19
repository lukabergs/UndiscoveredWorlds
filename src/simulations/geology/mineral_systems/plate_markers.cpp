#include "simulations/geology/mineral_systems/plate_markers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace mineralsystems
{
namespace
{
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

// Moves a continental-hosted marker to the centre of the nearest neighbouring cell (of the
// eight slots) where its host owns continental or transitional crust.
bool snaptohost(Marker& marker, const TectonicFrame& frame, std::uint32_t host)
{
    const mesh::CellGrid& grid = *frame.grid;
    const mesh::CellId here = grid.locate(marker.position);
    double nearest = -std::numeric_limits<double>::infinity();
    mesh::CellId best = mesh::nocell;
    for (int k = 1; k <= mesh::slotcount; ++k)
    {
        const std::int32_t neighbour = grid.neighbour(here, k);
        if (neighbour < 0)
            continue;
        const auto cell = static_cast<mesh::CellId>(neighbour);
        if (frame.plateIndex[cell] != host || !continentalcell(frame, cell))
            continue;
        // The largest cosine is the nearest centre (ties to the lower slot).
        const double along = mesh::dot(grid.centre(cell), marker.position);
        if (along > nearest)
        {
            nearest = along;
            best = cell;
        }
    }
    if (best == mesh::nocell)
        return false;
    marker.position = grid.centre(best);
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
    return frame.grid->locate(marker.position);
}

mesh::Vec3 movealong(const mesh::Vec3& position, const mesh::Vec3& displacement, double radius_m)
{
    // The displacement's part in the tangent plane at the position, then the great circle.
    const mesh::Vec3 tangent = displacement - mesh::dot(displacement, position) * position;
    const double distance = mesh::length(tangent);
    if (!(distance > 0.0))
        return position;
    const double angle = distance / radius_m;
    return mesh::normalise(std::cos(angle) * position + (std::sin(angle) / distance) * tangent);
}

void applymaterialshift(Marker& marker, const TectonicFrame& frame)
{
    // Approximate point tracers of the column displacement. Conservative
    // subcell mixing and root/sediment redistribution are not marker transport.
    constexpr double kMaximumSubstepSpacings = 0.25;
    const double radius_m = frame.grid->radius_m();
    const double spacing_m = frame.spacingKm * 1000.0;
    const double largest = mesh::length(frame.materialShift[markercell(marker, frame)]);
    const int steps = std::max(1, static_cast<int>(std::ceil(largest / (kMaximumSubstepSpacings * spacing_m))));
    for (int step = 0; step < steps; ++step)
    {
        const mesh::Vec3 shift = frame.materialShift[markercell(marker, frame)];
        marker.position = movealong(marker.position, (1.0 / steps) * shift, radius_m);
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
