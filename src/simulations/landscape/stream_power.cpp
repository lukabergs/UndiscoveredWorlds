#include "simulations/landscape/stream_power.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <string>
#include <vector>

// The scheme of step A2 (docs/architecture/erosion-hydrology-plan.md, section 3,
// Kernels). Per land cell i, with a_i its own area, A_i its drainage area,
// L_i its receiver distance (an outlet's to its sea cell), Q_i its discharge in
// km^3/yr and dt the step; an outlet's receiver elevation z_r is baselevel_m:
//
//   k_i = K_i Q_i^m / L_i                                 (1/Myr)
//   E_i = k_i max(0, z_i,new - z_r)                       fluvial incision (m/Myr)
//   D_i = G (S_up,i + E_i a_i) / (A_i + G a_i)            Yuan deposition (m/Myr)
//   (z_i,new - z_i,old) / dt = U_i - E_i + D_i
//
// where z_r is the receiver's elevation already updated in this sweep and S_up,i
// the sediment flux (m^3/Myr) arriving from upstream, taken from the previous
// sweep (zero in the first). D_i follows from Q_s,i = S_up,i + (E_i - D_i) a_i
// and D_i = G Q_s,i / A_i. Writing c_i = G a_i / (A_i + G a_i) and
// s_i = G S_up,i / (A_i + G a_i), D_i = s_i + c_i E_i, so
//   z_i,new - z_i,old = (U_i + s_i) dt - (1 - c_i) k_i dt (z_i,new - z_r)
// which is linear in z_i,new: with F'_i = k_i dt (1 - c_i) = k_i dt A_i / (A_i + G a_i)
// and b_i = z_i,old + (U_i + s_i) dt,
//   z_i,new = (b_i + F'_i z_r) / (1 + F'_i)   when b_i > z_r (then z_i,new > z_r),
//   z_i,new = b_i, E_i = 0                     otherwise (the cell is not above its
//                                              receiver, so it is not incised).
// With G = 0, c = s = 0 and this is Braun & Willett 2013:
// z_new = (z_old + U dt + F z_r) / (1 + F) with F = k dt. A land root (its own
// receiver, no sea cell) has k = 0 and keeps all the sediment reaching it.
//
// The glacial term (erosion plan Phase G, GlacialInput): on a land cell with ice
// share gamma = min(1, Q_ice / Q) the erodibility in k_i is
// K_eff = gamma^m K_g + (1 - gamma)^m K with K_g = erodibilityratio K, and z is the
// ice surface z + h (h the cell's thickness) throughout the solve. A cell with
// h > 0 is not incised below the surface baselevel_m - flotationratio h + h (its bed
// at the depth where the ice floats): where the solution would end below it, z_new
// is that surface and E_i = (b_i - z_new) / ((1 - c_i) dt), the incision the mass
// balance above gives for it, so the sediment budget stays exact; a cell whose b_i
// is already below it keeps z_new = b_i and is not incised. After the solve the bed
// is z_new - h. Without ice (gamma = 0, h = 0) no glacial arithmetic is done, so the
// step is bit-identical to the one without the input.

namespace landscape
{
namespace
{
constexpr double cubicmetrespercubickilometre = 1.0e9;
constexpr double squaremetrespersquarekilometre = 1.0e6;
constexpr double metresperkilometre = 1000.0;
// Stability limiter of the explicit hillslope step, carried over from the former
// climate_erosion.cpp (plan section 3): a face moves at most this fraction of its
// drop times the smaller of the two cell areas. Numerical, not a tuned value.
constexpr double hillslopefacelimit = 0.05;
constexpr int northdirection = 1;
constexpr int eastdirection = 3;
constexpr int southdirection = 5;
constexpr int westdirection = 7;

// The hillslope step visits every face shared by two cells once, from the side
// whose cardinal slot ranks first: east and south before north and west, then the
// lower cell id (a cube's frames rotate across face edges, so a face may be east or
// north on both sides). On latitude-longitude every face is a cell's east or south
// face, as before.
bool ownsface(int direction, int reversedirection, std::size_t cell, std::size_t neighbour)
{
    const auto rank = [](int slot) { return slot == eastdirection || slot == southdirection ? 0 : 1; };
    const int own = rank(direction), other = rank(reversedirection);
    return own != other ? own < other : cell < neighbour;
}

void checksize(std::size_t actual, std::size_t expected, const char* name)
{
    if (actual != expected)
        throw std::invalid_argument(std::string("stepstreampower: ") + name + " does not have metric.size() elements");
}

std::size_t checkedindex(int32_t value, std::size_t cells, const char* name)
{
    if (value < 0 || static_cast<std::size_t>(value) >= cells)
        throw std::invalid_argument(std::string("stepstreampower: ") + name + " is not a cell index");
    return static_cast<std::size_t>(value);
}
}

StreamPowerBudget stepstreampower(const SphereMetric& metric, const FlowGraph& graph,
                                  const std::vector<float>& erodibility, const std::vector<float>& uplift_m_myr,
                                  const StreamPowerParameters& parameters, std::vector<float>& elevation,
                                  std::vector<float>& erosion_m, std::vector<float>& deposition_m,
                                  std::vector<double>* sedimentflux_m3_myr, const GlacialInput* glacial)
{
    if (parameters.slopeexponent != 1.0)
        throw std::invalid_argument("stepstreampower: only a slope exponent of 1 is supported");
    const std::size_t cells = metric.size();
    checksize(graph.receiver.size(), cells, "graph.receiver");
    checksize(graph.direction.size(), cells, "graph.direction");
    checksize(graph.distance_m.size(), cells, "graph.distance_m");
    checksize(graph.outletsea.size(), cells, "graph.outletsea");
    checksize(graph.land.size(), cells, "graph.land");
    checksize(graph.area_m2.size(), cells, "graph.area_m2");
    checksize(graph.discharge_m3_yr.size(), cells, "graph.discharge_m3_yr");
    checksize(erodibility.size(), cells, "erodibility");
    checksize(uplift_m_myr.size(), cells, "uplift_m_myr");
    checksize(elevation.size(), cells, "elevation");
    checksize(erosion_m.size(), cells, "erosion_m");
    checksize(deposition_m.size(), cells, "deposition_m");
    if (glacial != nullptr)
    {
        if (glacial->iceflux_m3_yr == nullptr || glacial->thickness_m == nullptr)
            throw std::invalid_argument("stepstreampower: the glacial input lacks its ice flux or thickness");
        checksize(glacial->iceflux_m3_yr->size(), cells, "glacial iceflux_m3_yr");
        checksize(glacial->thickness_m->size(), cells, "glacial thickness_m");
    }

    // Every cell of `order` is land and drains to a land receiver, or is its own
    // receiver: an outlet into its outletsea cell (a sea cell), or a land root
    // without a sea cell (outletsea < 0: a region the flood never reached).
    for (int32_t entry : graph.order)
    {
        const std::size_t cell = checkedindex(entry, cells, "graph.order entry");
        if (graph.land[cell] == 0)
            throw std::invalid_argument("stepstreampower: graph.order lists a sea cell");
        const std::size_t receiver = checkedindex(graph.receiver[cell], cells, "graph.receiver");
        if (receiver == cell)
        {
            if (graph.outletsea[cell] >= 0
                && graph.land[checkedindex(graph.outletsea[cell], cells, "graph.outletsea of an outlet")] != 0)
                throw std::invalid_argument("stepstreampower: an outlet's outletsea is a land cell");
        }
        else if (graph.land[receiver] == 0)
            throw std::invalid_argument("stepstreampower: a land cell other than an outlet drains to a sea cell");
    }

    const double dt = parameters.step_myr;
    const double deposition = parameters.depositioncoefficient;  // G
    const double exponent = parameters.areaexponent;            // m

    // Double working copies: z is updated in place, zold keeps the start of the step.
    // With ice, z of a cell of the order is its ice surface until the solve ends.
    std::vector<double> z(elevation.begin(), elevation.end());
    std::vector<double> ice;  // h per cell (m); empty without glacial input
    if (glacial != nullptr)
    {
        ice.assign(cells, 0.0);
        for (int32_t entry : graph.order)
        {
            const std::size_t i = static_cast<std::size_t>(entry);
            const double thickness = static_cast<double>((*glacial->thickness_m)[i]);
            if (thickness > 0.0)
            {
                ice[i] = thickness;
                z[i] += thickness;
            }
        }
    }
    const std::vector<double> zold = z;

    // k_i = K_i Q_i^m / L_i in 1/Myr, Q in km^3/yr. An outlet's distance is to its
    // sea cell. A land root has no base level: its rate stays 0, so E = 0 and D = s.
    auto isroot = [&](std::size_t i) {
        return static_cast<std::size_t>(graph.receiver[i]) == i && graph.outletsea[i] < 0;
    };
    // With ice, K is K_eff and glacialshare holds gamma^m K_g / K_eff.
    std::vector<double> rate(cells, 0.0);
    std::vector<double> glacialshare;
    if (glacial != nullptr)
        glacialshare.assign(cells, 0.0);
    for (int32_t entry : graph.order)
    {
        const std::size_t i = static_cast<std::size_t>(entry);
        if (isroot(i))
            continue;
        const double length = static_cast<double>(graph.distance_m[i]);
        if (length <= 0.0)
            continue;
        double erodibilityhere = static_cast<double>(erodibility[i]);
        if (glacial != nullptr)
        {
            const double water = graph.discharge_m3_yr[i];
            const double gamma =
                water > 0.0 ? std::max(0.0, std::min(1.0, (*glacial->iceflux_m3_yr)[i] / water)) : 0.0;
            if (gamma > 0.0)
            {
                const double glacialpart = std::pow(gamma, exponent) * glacial->erodibilityratio * erodibilityhere;
                erodibilityhere = glacialpart + std::pow(1.0 - gamma, exponent) * erodibilityhere;
                if (erodibilityhere > 0.0)
                    glacialshare[i] = glacialpart / erodibilityhere;
            }
        }
        const double discharge = std::max(0.0, graph.discharge_m3_yr[i]) / cubicmetrespercubickilometre;
        rate[i] = erodibilityhere * std::pow(discharge, exponent) / length;
    }

    std::vector<double> incision(cells, 0.0);    // E_i, m/Myr
    std::vector<double> settling(cells, 0.0);    // D_i, m/Myr
    std::vector<double> upstream(cells, 0.0);    // S_up,i, m^3/Myr
    std::vector<double> exported(cells, 0.0);    // Q_s of an outlet, m^3/Myr

    // One Gauss-Seidel sweep, outlets first so that a receiver is updated before its donors.
    auto sweep = [&]() {
        for (int32_t entry : graph.order)
        {
            const std::size_t i = static_cast<std::size_t>(entry);
            const std::size_t receiver = static_cast<std::size_t>(graph.receiver[i]);
            // An outlet's receiver elevation is the sea surface (baselevel_m), not
            // the seabed of its outletsea cell; its sediment still enters outletsea.
            const double zr = receiver == i ? parameters.baselevel_m : z[receiver];
            const double area = metric.area_m2(i);
            const double drainage = graph.area_m2[i];
            double captured = 0.0;   // c_i
            double supplied = 0.0;   // s_i
            double retained = 1.0;   // 1 - c_i
            const double denominator = drainage + deposition * area;
            if (deposition > 0.0 && denominator > 0.0)
            {
                captured = deposition * area / denominator;
                supplied = deposition * upstream[i] / denominator;
                retained = drainage / denominator;
            }
            const double factor = rate[i] * dt * retained;  // F'_i
            const double base = zold[i] + (static_cast<double>(uplift_m_myr[i]) + supplied) * dt;
            double znew = base;
            double eroding = 0.0;
            if (base > zr)
            {
                znew = (base + factor * zr) / (1.0 + factor);
                eroding = rate[i] * std::max(0.0, znew - zr);  // rounding can leave znew a hair below zr
                if (!ice.empty() && ice[i] > 0.0)
                {
                    // The flotation floor of the bed, as an ice surface.
                    const double floor = (parameters.baselevel_m - glacial->flotationratio * ice[i]) + ice[i];
                    if (znew < floor && base > floor)
                    {
                        znew = floor;  // znew < base here, so factor, retained and dt are positive
                        eroding = (base - floor) / (retained * dt);
                    }
                    else if (znew < floor)
                    {
                        znew = base;  // the bed starts below the floor: not incised, not raised
                        eroding = 0.0;
                    }
                }
            }
            z[i] = znew;
            incision[i] = eroding;
            settling[i] = supplied + captured * eroding;
        }
    };

    // Sediment flux from the current E and D, upstream cells first: S_up of every
    // land cell and the flux Q_s = S_up + (E - D) a an outlet hands to the sea.
    auto accumulate = [&]() {
        std::fill(upstream.begin(), upstream.end(), 0.0);
        for (auto it = graph.order.rbegin(); it != graph.order.rend(); ++it)
        {
            const std::size_t i = static_cast<std::size_t>(*it);
            const std::size_t receiver = static_cast<std::size_t>(graph.receiver[i]);
            const double flux = upstream[i] + (incision[i] - settling[i]) * metric.area_m2(i);
            if (receiver == i)
                exported[i] = flux;
            else
                upstream[receiver] += flux;
        }
    };

    StreamPowerBudget budget;
    if (deposition > 0.0)
    {
        const int maxsweeps = std::max(1, parameters.depositionmaxsweeps);
        std::vector<double> previous(cells, 0.0);
        for (;;)
        {
            sweep();
            budget.depositionsweeps++;
            double change = std::numeric_limits<double>::infinity();
            if (budget.depositionsweeps > 1)
            {
                change = 0.0;
                for (int32_t entry : graph.order)
                {
                    const std::size_t i = static_cast<std::size_t>(entry);
                    change = std::max(change, std::abs(z[i] - previous[i]));
                }
            }
            if (change < parameters.depositiontolerance_m || budget.depositionsweeps >= maxsweeps)
                break;
            for (int32_t entry : graph.order)
                previous[static_cast<std::size_t>(entry)] = z[static_cast<std::size_t>(entry)];
            accumulate();
        }
    }
    else
        sweep();
    // The exported flux of the final E and D, so that the budget closes exactly.
    accumulate();
    // From here on z is the bed again: the solved surface less the ice.
    if (!ice.empty())
        for (int32_t entry : graph.order)
        {
            const std::size_t i = static_cast<std::size_t>(entry);
            if (ice[i] > 0.0)
                z[i] -= ice[i];
        }
    if (sedimentflux_m3_myr != nullptr)
    {
        std::vector<double>& flux = *sedimentflux_m3_myr;
        flux.assign(cells, 0.0);
        for (int32_t entry : graph.order)
        {
            const std::size_t i = static_cast<std::size_t>(entry);
            if (isroot(i))
                continue;
            flux[i] = static_cast<std::size_t>(graph.receiver[i]) == i
                ? exported[i]
                : upstream[i] + (incision[i] - settling[i]) * metric.area_m2(i);
        }
    }

    double erodedvolume = 0.0;
    double landdepositedvolume = 0.0;
    for (int32_t entry : graph.order)
    {
        const std::size_t i = static_cast<std::size_t>(entry);
        erodedvolume += incision[i] * dt * metric.area_m2(i);
        landdepositedvolume += settling[i] * dt * metric.area_m2(i);
    }

    // A land root keeps all the sediment that reaches it (it has no sea cell), so
    // the budget closes by construction.
    std::vector<double> trapped(cells, 0.0);  // m
    for (int32_t entry : graph.order)
    {
        const std::size_t i = static_cast<std::size_t>(entry);
        if (!isroot(i) || exported[i] == 0.0)
            continue;
        const double volume = exported[i] * dt;
        trapped[i] = volume / metric.area_m2(i);
        z[i] += trapped[i];
        landdepositedvolume += volume;
    }

    // Hillslope diffusion between land cells across every cardinal face, each once
    // (ownsface). All transfers use the elevations after the fluvial step and
    // are applied together, so the result does not depend on the face order.
    std::vector<double> lowered(cells, 0.0);  // m
    std::vector<double> raised(cells, 0.0);   // m
    const double diffusivity = parameters.hillslopediffusivity_km2_myr * squaremetrespersquarekilometre;  // m^2/Myr
    if (diffusivity > 0.0)
    {
        struct Transfer
        {
            std::size_t from;
            std::size_t to;
            double volume;
        };
        std::vector<Transfer> transfers;
        for (std::size_t i = 0; i < cells; i++)
        {
            if (graph.land[i] == 0)
                continue;
            for (int direction : { eastdirection, southdirection, northdirection, westdirection })
            {
                const std::int32_t neighbour = metric.neighbourid(i, direction);
                if (neighbour < 0)
                    continue;
                const std::size_t j = static_cast<std::size_t>(neighbour);
                if (graph.land[j] == 0 || !ownsface(direction, metric.reverse(i, direction), i, j))
                    continue;
                if (z[i] == z[j])
                    continue;
                const std::size_t high = z[i] > z[j] ? i : j;
                const std::size_t low = high == i ? j : i;
                const double drop = z[high] - z[low];
                double gradient = drop / metric.distance_m(i, direction);
                if (parameters.criticalslope > 0.0)
                    gradient = std::max(0.0, gradient - parameters.criticalslope);
                if (gradient <= 0.0)
                    continue;
                const double limit = hillslopefacelimit * drop * std::min(metric.area_m2(high), metric.area_m2(low));
                const double flux = diffusivity * dt * gradient * metric.facelength_m(i, direction);
                transfers.push_back({ high, low, std::min(limit, flux) });
            }
        }
        for (const Transfer& transfer : transfers)
        {
            lowered[transfer.from] += transfer.volume / metric.area_m2(transfer.from);
            raised[transfer.to] += transfer.volume / metric.area_m2(transfer.to);
            erodedvolume += transfer.volume;
            landdepositedvolume += transfer.volume;
        }
        for (std::size_t i = 0; i < cells; i++)
            z[i] += raised[i] - lowered[i];
    }

    // Marine deposition: each outlet's load enters its sea cell and walks down the
    // seabed's receivers, leaving 1 - exp(-d / travel) of what remains at every cell
    // it leaves over a distance d and the rest at the walk's end: a seabed pit, the
    // last sea cell before land, or the cell where the cap of one step per cell
    // (against a cycle) is reached. Progradation (erosion plan C5b), fill and spill:
    // no sea cell is raised above the delta top (baselevel_m + deltatop_m). A cell
    // keeps at most the room left below it (less what earlier outlets of this step
    // laid there); the rest of its share passes on with the load. At the walk's end
    // the cell keeps what fits and the remainder spills over the sea, lowest current
    // surface first (a priority walk from the end cell over its 8-connected sea
    // neighbours, through full cells), filling each cell it reaches up to the delta
    // top. Only when no sea cell with room is reachable does the remainder stay on
    // the end cell (counted as stranded), so the budget closes exactly.
    // Sea cells rise (or sink) by U dt when parameters.upliftsea is set: no fluvial
    // erosion there, and no budget entry (uplift moves rock in from below, it is not
    // eroded or deposited). It comes before the marine deposition so that the delta
    // top's room is measured on the seabed where this step leaves it: a deposit never
    // builds a sea cell above the delta top, as without uplift; the uplift itself may.
    // Outlets drain to baselevel_m, not to their sea cell, so the fluvial part of
    // the step does not depend on it.
    if (parameters.upliftsea)
        for (std::size_t i = 0; i < cells; i++)
            if (graph.land[i] == 0)
                z[i] += static_cast<double>(uplift_m_myr[i]) * dt;

    std::vector<double> marine(cells, 0.0);  // m^3
    const double travel = parameters.sedimenttravel_km * metresperkilometre;
    const double deltatop = parameters.baselevel_m + parameters.deltatop_m;
    auto room = [&](std::size_t cell) {
        if (!std::isfinite(deltatop))
            return std::numeric_limits<double>::infinity();
        return std::max(0.0, (deltatop - z[cell]) * metric.area_m2(cell) - marine[cell]);
    };
    std::vector<uint32_t> reached;  // spill stamps, allocated at the first spill
    uint32_t stamp = 0;
    using Candidate = std::pair<double, std::size_t>;  // current surface, cell
    auto spill = [&](std::size_t start, double load) {
        if (reached.empty())
            reached.assign(cells, 0);
        stamp++;
        std::priority_queue<Candidate, std::vector<Candidate>, std::greater<Candidate>> queue;
        auto reach = [&](std::size_t from) {
            for (int direction = 1; direction <= directioncount; direction++)
            {
                const std::int32_t neighbour = metric.neighbourid(from, direction);
                if (neighbour < 0)
                    continue;
                const std::size_t next = static_cast<std::size_t>(neighbour);
                if (graph.land[next] != 0 || reached[next] == stamp)
                    continue;
                reached[next] = stamp;
                queue.push({ z[next] + marine[next] / metric.area_m2(next), next });
            }
        };
        reached[start] = stamp;
        reach(start);
        while (load > 0.0 && !queue.empty())
        {
            const std::size_t cell = queue.top().second;
            queue.pop();
            const double placed = std::min(load, room(cell));
            marine[cell] += placed;
            load -= placed;
            reach(cell);
        }
        return load;
    };
    for (int32_t entry : graph.order)
    {
        const std::size_t outlet = static_cast<std::size_t>(entry);
        if (static_cast<std::size_t>(graph.receiver[outlet]) != outlet || isroot(outlet))
            continue;
        double load = exported[outlet] * dt;
        if (load == 0.0)
            continue;
        std::size_t cell = static_cast<std::size_t>(graph.outletsea[outlet]);
        for (std::size_t step = 0;; step++)
        {
            const std::size_t next = checkedindex(graph.receiver[cell], cells, "graph.receiver of a sea cell");
            if (next == cell || graph.land[next] != 0 || step >= cells)
            {
                const double kept = std::min(load, room(cell));
                marine[cell] += kept;
                load -= kept;
                if (load > 0.0)
                    load = spill(cell, load);
                if (load > 0.0)
                {
                    marine[cell] += load;
                    budget.strandedloads++;
                    budget.stranded_km3 += load / cubicmetrespercubickilometre;
                }
                break;
            }
            const double distance = static_cast<double>(graph.distance_m[cell]);
            const double fraction = travel > 0.0 ? 1.0 - std::exp(-distance / travel) : 1.0;
            const double settled = std::min(fraction * load, room(cell));
            marine[cell] += settled;
            load -= settled;
            cell = next;
        }
    }
    double marinevolume = 0.0;
    for (std::size_t i = 0; i < cells; i++)
    {
        if (marine[i] == 0.0)
            continue;
        z[i] += marine[i] / metric.area_m2(i);
        marinevolume += marine[i];
    }

    for (std::size_t i = 0; i < cells; i++)
    {
        elevation[i] = static_cast<float>(z[i]);
        if (graph.land[i] != 0)
        {
            erosion_m[i] = static_cast<float>(incision[i] * dt + lowered[i]);
            deposition_m[i] = static_cast<float>(settling[i] * dt + trapped[i] + raised[i]);
        }
        else
        {
            erosion_m[i] = 0.0f;
            deposition_m[i] = static_cast<float>(marine[i] / metric.area_m2(i));
        }
    }
    if (glacial != nullptr && glacial->glacialerosion_m != nullptr)
    {
        // The glacial part of each cell's fluvial-glacial incision, never more than its erosion.
        std::vector<float>& glacialerosion = *glacial->glacialerosion_m;
        glacialerosion.assign(cells, 0.0f);
        for (int32_t entry : graph.order)
        {
            const std::size_t i = static_cast<std::size_t>(entry);
            if (glacialshare[i] <= 0.0)
                continue;
            const float part = static_cast<float>(glacialshare[i] * incision[i] * dt);
            glacialerosion[i] = std::min(std::max(part, 0.0f), erosion_m[i]);
        }
    }

    budget.eroded_km3 = erodedvolume / cubicmetrespercubickilometre;
    budget.landdeposited_km3 = landdepositedvolume / cubicmetrespercubickilometre;
    budget.marinedeposited_km3 = marinevolume / cubicmetrespercubickilometre;
    return budget;
}
}
