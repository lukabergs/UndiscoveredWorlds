#ifndef TECTONIC_TOPOLOGY_HPP
#define TECTONIC_TOPOLOGY_HPP
#include "simulations/geology/plate_tectonics/tectonic_transport.hpp"
#include "simulations/geology/plate_tectonics/tectonic_dynamics.hpp"
#include "simulations/geology/plate_tectonics/tectonic_grid.hpp"
#include <memory>
#include <cstdint>
#include <map>
#include <vector>

namespace tectonics::topology {
enum class EventKind { Split, Fusion };
struct Event {
    EventKind kind = EventKind::Split;
    double time_myr = 0;
    uint32_t parent_a = 0, parent_b = 0, result = 0; // Persistent identities, not array indices.
    uint32_t cells = 0;
};
struct Change { EventKind kind; uint32_t a, b; }; // Live array indices, before compaction.
// Effective failure/welding rules, not a resolved damage rheology or Earth fit.
inline constexpr double rift_yield_n_m = 2e12;
inline constexpr double rift_damage_time_myr = 30;
inline constexpr double event_cooldown_myr = 20;
inline constexpr double weld_time_myr = 30;
inline constexpr double weld_speed_m_s = .2 / 3.15576e9;
inline constexpr double minimum_plate_area_m2 = 1e12;
struct Parameters {
    bool progressive_rifting = false;
    double yield_n_m = rift_yield_n_m;
    double damage_time_myr = rift_damage_time_myr;
    double strain_rate_myr = .002; // Maximum distributed extension per Myr at fully damaged rifts.
};
struct Diagnostics {
    bool enabled = false;
    uint64_t reassigned_cells = 0; // Cumulative connectivity repairs, not crust loss.
    uint32_t max_detached_components = 0;
    std::vector<uint32_t> stable_ids;
    std::vector<Event> events;
    double max_rift_damage = 0;
    double max_weld_time_myr = 0;
    double max_tension_n_m = 0;
    uint32_t damaged_rift_cells = 0;
};
// Keeps subcell plate coverage between updates. Ownership is separate from
// material columns: resolving a boundary must never erase or flatten terrain.
class Model {
public:
    // Neighbourhoods, connectivity and contacts on the engine's cell table (W2-E); the rift
    // geometry and the coverage transport are still planar (E2).
    void initialize(std::shared_ptr<const grid::TectonicGrid> grid, uint32_t plates,
                    const std::vector<uint32_t>& owners);
    // The same on a planar scaffold of width x height (unit tests).
    void initialize(uint32_t width, uint32_t height, uint32_t plates, const std::vector<uint32_t>& owners);
    std::vector<uint32_t> advance(const std::vector<transport::Displacement>& plate_displacement);
    // Changes ownership only. No height, rock, sediment or age input can be modified.
    std::vector<Change> evolve(std::vector<uint32_t>& owners,
        const std::vector<dynamics::Vector>& mantle_m_s,
        const std::vector<dynamics::Vector>& velocity_m_s,
        const std::vector<double>& continental_fraction,
        double cell_m, double dt_myr, double time_myr,
        const Parameters& parameters = {}, const std::vector<double>& strength = {});
    // Incipient rifting deforms material before a new rigid plate identity exists.
    std::vector<transport::Displacement> riftDisplacement(const std::vector<uint32_t>& owners,
        double dt_myr, const Parameters& parameters) const;
    std::vector<double> damageField(const std::vector<uint32_t>& owners) const;
    void removePlate(uint32_t index); // Same swap-with-last convention as lithosphere.
    const Diagnostics& diagnostics() const { return diagnostics_; }
private:
    void connect(std::vector<uint32_t>& owners);
    std::vector<uint32_t> partition(uint32_t plate, const std::vector<uint32_t>& owners,
        const std::vector<dynamics::Vector>& mantle, const std::vector<double>& strength = {}) const;
    struct Rift { std::vector<double> side, cracks; double damage = 0, duration = 0, cooldown = 0; };
    std::shared_ptr<const grid::TectonicGrid> grid_;
    uint32_t width_ = 0, height_ = 0;
    std::vector<std::vector<double>> coverage_;
    std::vector<Rift> rifts_;
    std::map<std::pair<uint32_t, uint32_t>, double> welding_;
    uint32_t next_id_ = 0;
    Diagnostics diagnostics_;
};
}
#endif
