#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <queue>
#include <tuple>

namespace tectonics::topology {
namespace {
constexpr uint32_t missing = std::numeric_limits<uint32_t>::max();
bool finite(const mesh::Vec3& v) { return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z); }
// The unit tangent of v at the unit vector r (zero for a radial v).
mesh::Vec3 tangent(const mesh::Vec3& v, const mesh::Vec3& r) {
    const mesh::Vec3 t = v - mesh::dot(v, r) * r;
    const double length = mesh::length(t);
    return length > 0 ? (1.0 / length) * t : mesh::Vec3{};
}
}
void Model::initialize(std::shared_ptr<const grid::TectonicGrid> grid, uint32_t plates,
        const std::vector<uint32_t>& owners) {
    if (!grid || grid->facecells() < 2 || plates == 0 || owners.size() != grid->cellcount())
        throw std::runtime_error("Invalid plate topology grid");
    grid_ = std::move(grid);
    lattice_ = sphere_transport::Lattice::make(grid_->facecells(), 1, grid_->radius_m());
    diagnostics_ = {}; diagnostics_.enabled = true;
    coverage_.assign(plates, std::vector<double>(owners.size(), 0));
    rifts_.assign(plates, {}); welding_.clear(); next_id_ = plates;
    diagnostics_.stable_ids.resize(plates);
    std::iota(diagnostics_.stable_ids.begin(), diagnostics_.stable_ids.end(), 0U);
    for (size_t i = 0; i < owners.size(); ++i) {
        if (owners[i] >= plates) throw std::runtime_error("Invalid plate topology owner");
        coverage_[owners[i]][i] = 1;
    }
}
void Model::removePlate(uint32_t index) {
    if (index >= coverage_.size()) throw std::runtime_error("Invalid removed topology plate");
    coverage_[index] = std::move(coverage_.back());
    coverage_.pop_back();
    rifts_[index] = std::move(rifts_.back()); rifts_.pop_back();
    diagnostics_.stable_ids[index] = diagnostics_.stable_ids.back();
    diagnostics_.stable_ids.pop_back();
}
std::vector<uint32_t> Model::advance(const std::vector<mesh::Vec3>& rotation) {
    if (!diagnostics_.enabled || rotation.size() != coverage_.size())
        throw std::runtime_error("Invalid plate topology motion");
    for (const auto& w : rotation) if (!finite(w)) throw std::runtime_error("Invalid plate topology rotation");
    const size_t size = grid_->cellcount();
    for (size_t p = 0; p < coverage_.size(); ++p) {
        // The plate's coverage (and its rift fields) moved by its rotation alone, over the
        // cells it covers: the per-plate kernel keeps a rigid interior at exactly 1.
        sphere_transport::Motion motion;
        motion.rotation = {rotation[p]};
        motion.plate.assign(size, sphere_transport::noplate);
        for (size_t i = 0; i < size; ++i) if (coverage_[p][i] > 0) motion.plate[i] = 0;
        const sphere_transport::Stencil stencil(lattice_, motion);
        std::vector<std::vector<double>*> fields{&coverage_[p]};
        if (!rifts_[p].side.empty()) fields.push_back(&rifts_[p].side);
        if (!rifts_[p].cracks.empty()) fields.push_back(&rifts_[p].cracks);
        std::vector<const std::vector<double>*> inputs(fields.begin(), fields.end());
        const auto layers = stencil.transport(inputs);
        for (auto* field : fields) std::fill(field->begin(), field->end(), 0.0);
        for (const auto& layer : layers)
            for (size_t k = 0; k < layer.cells.size(); ++k)
                for (size_t f = 0; f < fields.size(); ++f)
                    (*fields[f])[layer.cells[k]] += layer.values[f][k] * layer.coverage[k];
    }
    std::vector<uint32_t> owners(size, 0);
    for (size_t i = 0; i < size; ++i) {
        double sum = 0, best = -1;
        for (uint32_t p = 0; p < coverage_.size(); ++p) {
            sum += coverage_[p][i];
            if (coverage_[p][i] > best) { best = coverage_[p][i]; owners[i] = p; }
        }
        if (!(sum > 0)) throw std::runtime_error("Uncovered plate topology cell");
        // Normalize competing coverage, preserving subcell boundary position.
        // New ocean crust and continental stretching are handled by columns.
        for (size_t p = 0; p < coverage_.size(); ++p) {
            if (!rifts_[p].side.empty()) rifts_[p].side[i] = std::min(rifts_[p].side[i], coverage_[p][i]) / sum;
            if (!rifts_[p].cracks.empty()) rifts_[p].cracks[i] /= sum;
            coverage_[p][i] /= sum;
        }
    }
    connect(owners);
    return owners;
}

std::vector<uint32_t> Model::partition(uint32_t plate, const std::vector<uint32_t>& owners,
        const std::vector<mesh::Vec3>& mantle, const std::vector<double>& strength) const {
    // Two connected traction domains, not random seeds or a new global Voronoi
    // map. The strongest varying mantle-flow axis selects opposite anchors; the flow is read
    // on the tangent plane at the plate's centroid (east, north there).
    mesh::Vec3 centre;
    uint32_t count = 0;
    for (size_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        centre = centre + grid_->area_m2(static_cast<grid::CellId>(i)) * grid_->position(static_cast<grid::CellId>(i));
        ++count;
    }
    if (count == 0 || !(mesh::length(centre) > 0)) return {};
    centre = mesh::normalise(centre);
    const double latitude = mesh::latitudeof(centre), longitude = mesh::longitudeof(centre);
    const mesh::Vec3 east = mesh::eastof(latitude, longitude), north = mesh::northof(latitude, longitude);
    const auto planar = [&](size_t i) { return std::array<double, 2>{mesh::dot(mantle[i], east), mesh::dot(mantle[i], north)}; };
    double mean_x = 0, mean_y = 0;
    for (size_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        const auto v = planar(i); mean_x += v[0]; mean_y += v[1];
    }
    mean_x /= count; mean_y /= count;
    double xx = 0, xy = 0, yy = 0;
    for (size_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        const auto v = planar(i);
        const double x = v[0] - mean_x, y = v[1] - mean_y;
        xx += x*x; xy += x*y; yy += y*y;
    }
    if (xx + yy < 1e-30) return {};
    const double angle = .5 * std::atan2(2*xy, xx-yy);
    const double nx = std::cos(angle), ny = std::sin(angle);
    uint32_t low = missing, high = missing;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (uint32_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        const auto v = planar(i);
        const double projection = v[0] * nx + v[1] * ny;
        if (projection < lo) { lo = projection; low = i; }
        if (projection > hi) { hi = projection; high = i; }
    }
    if (low == high) return {};
    std::vector<uint32_t> side(owners.size(), missing);
    if (!strength.empty()) {
        // Strength-weighted watershed: strong interiors join the traction
        // domains first, leaving the interface preferentially in weak rock.
        // This is a coarse fault-path approximation, not an elastic stress solve.
        using Entry=std::tuple<double,uint32_t,uint32_t>;
        std::priority_queue<Entry,std::vector<Entry>,std::greater<Entry>> queue;
        queue.emplace(0,low,0);queue.emplace(0,high,1);
        while (!queue.empty()) {
            const auto [distance,i,label]=queue.top();queue.pop();
            if (side[i]!=missing) continue;
            side[i]=label;
            for (auto n:grid_->neighbours(i))
                if (owners[n]==plate && side[n]==missing)
                    queue.emplace(distance+2/(strength[i]+strength[n]),n,label);
        }
        return side;
    }
    // Each cell to the nearer anchor by great-circle distance: the cut is the anchors' bisecting
    // great circle wherever the plate lies (a breadth-first split counts grid steps, whose
    // bisector bends and lengthens across a cube corner; E3).
    const mesh::Vec3 a = grid_->position(low), b = grid_->position(high);
    for (uint32_t i = 0; i < owners.size(); ++i)
        if (owners[i] == plate) side[i] = mesh::dot(grid_->position(i), b) > mesh::dot(grid_->position(i), a) ? 1 : 0;
    return side;
}

std::vector<Change> Model::evolve(std::vector<uint32_t>& owners,
        const std::vector<mesh::Vec3>& mantle, const std::vector<mesh::Vec3>& rotation,
        const std::vector<double>& continent, double dt, double time,
        const Parameters& parameters, const std::vector<double>& strength) {
    const size_t size = grid_ ? grid_->cellcount() : 0, plates = coverage_.size();
    if (!std::isfinite(parameters.yield_n_m) || parameters.yield_n_m<=0 ||
        !std::isfinite(parameters.damage_time_myr) || parameters.damage_time_myr<=0 ||
        !std::isfinite(parameters.strain_rate_myr) || parameters.strain_rate_myr<0 ||
        (!strength.empty() && strength.size()!=size)) throw std::runtime_error("Invalid rift parameters");
    for (double v:strength) if (!std::isfinite(v) || v<=0) throw std::runtime_error("Invalid rock strength");
    if (!diagnostics_.enabled || owners.size() != size || mantle.size() != size || continent.size() != size ||
        rotation.size() != plates || !std::isfinite(dt) || dt <= 0 || !std::isfinite(time) || time < 0)
        throw std::runtime_error("Invalid topology evolution inputs");
    for (size_t i = 0; i < size; ++i)
        if (owners[i] >= plates || !finite(mantle[i]) ||
            !std::isfinite(continent[i]) || continent[i] < 0 || continent[i] > 1 + 1e-8)
            throw std::runtime_error("Invalid topology forcing cell");
    for (const auto& w : rotation) if (!finite(w)) throw std::runtime_error("Invalid topology velocity");
    const double radius = grid_->radius_m();
    const double cell_area = grid_->totalarea_m2() / static_cast<double>(size);
    const auto neighbors = [&](uint32_t i) { return grid_->neighbours(i); };
    const auto& edges = grid_->edges();
    for (size_t p=0;p<rifts_.size();++p) {
        auto& r=rifts_[p];r.cooldown=std::max(0.0,r.cooldown-dt);
        // Material damage survives a rejected geometric candidate. Heal even
        // while no connected candidate exists, rather than erasing the fault.
        if (parameters.progressive_rifting) for (size_t i=0;i<r.cracks.size();++i)
            r.cracks[i]=std::min(r.cracks[i],coverage_[p][i])*std::exp(-dt/200);
    }
    const auto reset_candidate=[&](Rift& r) {
        if (parameters.progressive_rifting) {r.side.clear();r.damage=0;r.duration=0;}
        else r={};
    };
    diagnostics_.max_rift_damage = 0; diagnostics_.max_weld_time_myr = 0;
    diagnostics_.max_tension_n_m=0;diagnostics_.damaged_rift_cells=0;
    const auto count_damage=[&]() {
        diagnostics_.damaged_rift_cells=0;
        for (size_t i=0;i<size;++i) {
            const auto p=owners[i];const auto& cracks=rifts_[p].cracks;
            diagnostics_.damaged_rift_cells+=!cracks.empty() &&
                cracks[i]>.1*std::max(coverage_[p][i],1e-12);
        }
    };
    // Each contact's edges: closing ((w_a - w_b) x R x) . n at the edge's midpoint, the
    // relative speed there and the shared length.
    struct Contact { uint32_t faces=0, continental=0; double compression=0, normal_motion=0, speed=0, length_m=0; };
    std::map<std::pair<uint32_t,uint32_t>,Contact> contacts;
    for (const grid::Edge& edge : edges) {
        const auto i = edge.a, j = edge.b;
        auto a = owners[i], b = owners[j];
        if (a == b) continue;
        // From the lower plate index toward the higher (the contact's orientation).
        mesh::Vec3 normal = edge.normal;
        if (a > b) { std::swap(a, b); normal = -1.0 * normal; }
        auto& c = contacts[{a, b}]; ++c.faces;
        c.continental += continent[i] >= .5 && continent[j] >= .5;
        const mesh::Vec3 relative = mesh::cross(rotation[a] - rotation[b], radius * edge.midpoint);
        const double closing = mesh::dot(relative, normal);
        c.compression += std::max(0.0, closing);
        c.normal_motion += std::abs(closing);
        c.speed += mesh::length(relative);
        c.length_m += edge.length_m;
    }
    std::map<std::pair<uint32_t,uint32_t>,double> next_welding;
    for (const auto& entry : contacts) {
        const auto [a,b] = entry.first;
        const auto& c = entry.second;
        const auto key = std::minmax(diagnostics_.stable_ids[a], diagnostics_.stable_ids[b]);
        const double speed = c.speed / c.faces;
        const double previous = welding_.count(key) ? welding_.at(key) : 0;
        // Require broad continental contact, compatible full velocity (including
        // transform slip), and evidence of compression. Mere proximity is not welding.
        const bool locked = c.faces >= 4 && c.length_m >= 500000 &&
            c.continental >= .8*c.faces && speed <= weld_speed_m_s &&
            (c.compression > c.faces*1e-13 || (previous > 0 && speed < 1e-13)) &&
            c.normal_motion >= .25*c.faces*speed && rifts_[a].cooldown == 0 && rifts_[b].cooldown == 0;
        const double duration = locked ? previous + dt : 0;
        next_welding[key] = duration;
        diagnostics_.max_weld_time_myr = std::max(diagnostics_.max_weld_time_myr, duration);
        if (duration < weld_time_myr) continue;
        if (next_id_ == missing) throw std::runtime_error("Plate identity exhausted");
        uint32_t transferred = 0;
        for (size_t i = 0; i < size; ++i) {
            coverage_[a][i] += coverage_[b][i]; coverage_[b][i] = 0;
            if (owners[i] == b) { owners[i] = a; ++transferred; }
        }
        diagnostics_.events.push_back({EventKind::Fusion, time, diagnostics_.stable_ids[a],
            diagnostics_.stable_ids[b], next_id_, transferred});
        diagnostics_.stable_ids[a] = next_id_++;
        rifts_[a] = {}; rifts_[a].cooldown = event_cooldown_myr; rifts_[b] = {};
        // Retired identities cannot carry welding history into a different pair.
        welding_.clear();
        count_damage();
        return {{EventKind::Fusion,a,b}};
    }
    welding_ = std::move(next_welding);
    const auto min_cells = static_cast<uint32_t>(std::max(16.0, std::ceil(minimum_plate_area_m2/cell_area)));
    for (uint32_t p = 0; p < plates; ++p) {
        auto& r = rifts_[p];
        if (r.cooldown > 0) continue;
        const auto count = static_cast<uint32_t>(std::count(owners.begin(), owners.end(), p));
        if (count < 2*min_cells) continue;
        if (r.side.empty()) {
            auto inherited_strength=parameters.progressive_rifting?strength:std::vector<double>{};
            if (parameters.progressive_rifting && !r.cracks.empty()) {
                if (inherited_strength.empty()) inherited_strength.assign(size,1);
                for (size_t i=0;i<size;++i) inherited_strength[i]*=1-.6*
                    std::clamp(r.cracks[i]/std::max(coverage_[p][i],1e-12),0.0,1.0);
            }
            const auto sides = partition(p, owners, mantle,inherited_strength);
            if (sides.empty()) continue;
            r.side.resize(size, 0);
            for (size_t i = 0; i < size; ++i) if (sides[i] == 1) r.side[i] = coverage_[p][i];
            if (parameters.progressive_rifting && r.cracks.empty()) r.cracks.assign(size,0);
        }
        std::vector<uint32_t> side(size, missing);
        uint32_t counts[2]{}, anchors[2]{missing,missing};
        double areas[2]{};
        mesh::Vec3 flow[2]{};
        for (uint32_t i = 0; i < size; ++i) if (owners[i] == p) {
            const uint32_t s = r.side[i] > .5*coverage_[p][i] ? 1U : 0U;
            side[i] = s; ++counts[s]; anchors[s] = i; areas[s] += grid_->area_m2(i);
            flow[s] = flow[s] + mantle[i];
        }
        if (std::min(counts[0],counts[1]) < min_cells || std::min(counts[0],counts[1]) < .15*count) {
            reset_candidate(r); continue;
        }
        // A throughgoing cut must leave exactly two connected resolved domains.
        // Advected candidates that develop islands are discarded, not promoted.
        std::vector<uint8_t> visited(size,0);
        bool connected = true;
        mesh::Vec3 center[2]{};
        for (uint32_t s = 0; s < 2; ++s) {
            std::deque<uint32_t> queue{anchors[s]}; visited[anchors[s]] = 1;
            uint32_t found = 0;
            while (!queue.empty()) {
                const auto i = queue.front(); queue.pop_front(); ++found;
                center[s] = center[s] + grid_->area_m2(i) * grid_->position(i);
                for (auto n : neighbors(i)) if (!visited[n] && side[n] == s) { visited[n]=1; queue.push_back(n); }
            }
            connected &= found == counts[s];
            flow[s] = (1.0 / counts[s]) * flow[s];
        }
        if (!connected || !(mesh::length(center[0]) > 0) || !(mesh::length(center[1]) > 0)) { reset_candidate(r); continue; }
        // The direction from side 0's centroid toward side 1's, at their midpoint; their
        // great-circle distance in cell spacings.
        const mesh::Vec3 c0 = mesh::normalise(center[0]), c1 = mesh::normalise(center[1]);
        const double distance = mesh::greatcircle_rad(c0, c1) * radius / grid_->cellspacing_m();
        if (distance < 1) { reset_candidate(r); continue; }
        const mesh::Vec3 middle = mesh::normalise(c0 + c1);
        const mesh::Vec3 n = tangent(c1 - c0, middle);
        // The fault's length: each cut edge projected on the cut's normal there (the direction
        // between the sides' centroids), so a cut across the grid's axes (a staircase of
        // edges) measures its own length, not the staircase's (E3).
        double interface_m = 0;
        for (const grid::Edge& edge : edges)
            if (side[edge.a] != missing && side[edge.b] != missing && side[edge.a] != side[edge.b])
                interface_m += edge.length_m * std::abs(mesh::dot(edge.normal, tangent(c1 - c0, edge.midpoint)));
        const double opening = mesh::dot(flow[1] - flow[0], n);
        const double reduced_area = areas[0]*areas[1]/(areas[0]+areas[1]);
        const double tension = interface_m > 0 ? dynamics::kBasalDrag*reduced_area*opening/interface_m : 0;
        diagnostics_.max_tension_n_m=std::max(diagnostics_.max_tension_n_m,tension);
        if (parameters.progressive_rifting) {
            double total=0;uint32_t interface_cells=0;
            auto next=r.cracks;
            for (uint32_t i=0;i<size;++i) if (owners[i]==p) {
                bool interface=false;double neighbor_damage=0;
                for (auto nb:neighbors(i)) if (owners[nb]==p) {
                    interface|=side[nb]!=side[i];
                    neighbor_damage=std::max(neighbor_damage,r.cracks[nb]/std::max(coverage_[p][nb],1e-12));
                }
                double damage=r.cracks[i]/std::max(coverage_[p][i],1e-12);
                // Inherited thin/warm rock and nearby damage reduce strength.
                // Damage grows only under tension, heals slowly under unloading,
                // and is carried with the plate rather than reseeded each update.
                const double rock=strength.empty()?1:strength[i];
                const double yield=parameters.yield_n_m*rock*(1-.6*damage);
                const double growth=interface?std::max(0.0,tension/yield-1):0;
                damage=std::clamp(damage+
                    growth*(.2+.8*neighbor_damage)*dt/parameters.damage_time_myr,0.0,1.0);
                next[i]=damage*coverage_[p][i];
                if (interface) {total+=damage;++interface_cells;}
            }
            r.cracks=std::move(next);
            r.damage=interface_cells?total/interface_cells:0;
            r.duration=opening>0?r.duration+dt:0;
        } else if (tension > parameters.yield_n_m) {
            r.damage += (tension/parameters.yield_n_m-1)*dt/parameters.damage_time_myr;
            r.duration += dt;
        } else {
            r.damage *= std::exp(-dt/20); r.duration = 0;
            if (r.damage < 1e-6) { r = {}; continue; }
        }
        diagnostics_.max_rift_damage = std::max(diagnostics_.max_rift_damage, r.damage);
        if (r.damage < (parameters.progressive_rifting?.9:1) || r.duration < event_cooldown_myr) continue;
        if (next_id_ == missing) throw std::runtime_error("Plate identity exhausted");
        const auto daughter = static_cast<uint32_t>(plates);
        coverage_.emplace_back(size,0);
        for (size_t i = 0; i < size; ++i) {
            // Split the parent's fractional coverage; other plates are untouched.
            const double fraction = r.side[i] > .5*coverage_[p][i] ? coverage_[p][i] : 0;
            coverage_[daughter][i] = fraction; coverage_[p][i] -= fraction;
            if (side[i] == 1) owners[i] = daughter;
        }
        diagnostics_.events.push_back({EventKind::Split,time,diagnostics_.stable_ids[p],
            diagnostics_.stable_ids[p],next_id_,counts[1]});
        diagnostics_.stable_ids.push_back(next_id_++);
        r = {}; r.cooldown = event_cooldown_myr;
        rifts_.push_back({}); rifts_.back().cooldown = event_cooldown_myr;
        welding_.clear();
        count_damage();
        return {{EventKind::Split,p,daughter}};
    }
    count_damage();
    return {};
}
std::vector<double> Model::damageField(const std::vector<uint32_t>& owners) const {
    std::vector<double> result(owners.size());
    for (size_t i=0;i<owners.size();++i) {
        const auto p=owners[i];
        if (!rifts_[p].cracks.empty()) result[i]=std::clamp(rifts_[p].cracks[i]/std::max(coverage_[p][i],1e-12),0.0,1.0);
    }
    return result;
}

std::vector<mesh::Vec3> Model::riftDisplacement(const std::vector<uint32_t>& owners,
        double dt,const Parameters& parameters) const {
    std::vector<mesh::Vec3> result(owners.size());
    if (!parameters.progressive_rifting || parameters.strain_rate_myr==0) return result;
    const auto damage=damageField(owners);
    // Opening motion points away from an incipient fault. Conservative column
    // transport converts this distributed extension to thinning and subsidence;
    // no rift depth is painted onto the terrain and no rock is deleted. Each damaged face
    // across the fault moves the cell away from its neighbour by the opening strain times
    // their centre distance, along the line between the centres.
    for (uint32_t i=0;i<owners.size();++i) {
        const auto p=owners[i];
        if (rifts_[p].side.empty() || damage[i]==0) continue;
        const bool side=rifts_[p].side[i]>.5*coverage_[p][i];
        const mesh::Vec3 r=grid_->position(i);
        for (std::uint8_t k=0;k<grid::directioncount;++k) {
            const uint32_t j=grid_->neighbour(i,k);
            if (owners[j]!=p || (rifts_[p].side[j]>.5*coverage_[p][j])==side) continue;
            const double opening=.5*parameters.strain_rate_myr*dt*std::min(damage[i],damage[j]);
            const double distance=grid_->edges()[grid_->edge(i,k)].distance_m;
            result[i]=result[i]-(opening*distance)*tangent(grid_->position(j)-r,r);
        }
    }
    return result;
}

void Model::connect(std::vector<uint32_t>& owners) {
    const auto neighbors = [&](uint32_t i) { return grid_->neighbours(i); };
    std::vector<uint32_t> component(owners.size(), missing), largest(coverage_.size(), missing);
    std::vector<uint32_t> sizes, plate_ids;
    std::deque<uint32_t> queue;
    for (uint32_t i = 0; i < owners.size(); ++i) {
        if (component[i] != missing) continue;
        const auto id = static_cast<uint32_t>(sizes.size());
        uint32_t count = 0;
        queue.push_back(i); component[i] = id;
        while (!queue.empty()) {
            const auto cell = queue.front(); queue.pop_front(); ++count;
            for (auto n : neighbors(cell)) if (component[n] == missing && owners[n] == owners[i]) {
                component[n] = id; queue.push_back(n);
            }
        }
        sizes.push_back(count); plate_ids.push_back(owners[i]);
        auto& core = largest[owners[i]];
        if (core == missing || count > sizes[core]) core = id;
    }
    uint32_t detached = 0;
    for (uint32_t c = 0; c < sizes.size(); ++c) detached += largest[plate_ids[c]] != c;
    diagnostics_.max_detached_components = std::max(diagnostics_.max_detached_components, detached);
    if (detached == 0) return;
    // Retain every largest component. Simultaneous flood from these connected
    // cores assigns detached fragments to adjacent plates without creating new
    // fragments (unlike sequential orphan-to-orphan reassignment). This is a
    // numerical grid repair, not a claim to simulate physical rifting or fusion.
    std::vector<uint32_t> connected(owners.size(), missing);
    for (uint32_t i = 0; i < owners.size(); ++i)
        if (component[i] == largest[owners[i]]) { connected[i] = owners[i]; queue.push_back(i); }
    while (!queue.empty()) {
        const auto cell = queue.front(); queue.pop_front();
        for (auto n : neighbors(cell)) if (connected[n] == missing) {
            connected[n] = connected[cell]; queue.push_back(n);
        }
    }
    for (size_t i = 0; i < owners.size(); ++i) if (owners[i] != connected[i]) {
        ++diagnostics_.reassigned_cells;
        for (auto& plate : coverage_) plate[i] = 0;
        for (auto& r : rifts_) if (!r.side.empty()) r.side[i] = 0;
        for (auto& r : rifts_) if (!r.cracks.empty()) r.cracks[i] = 0;
        coverage_[connected[i]][i] = 1;
    }
    owners = std::move(connected);
}
}
