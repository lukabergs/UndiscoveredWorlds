#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"
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
namespace { constexpr uint32_t missing = std::numeric_limits<uint32_t>::max(); }
void Model::initialize(uint32_t width, uint32_t height, uint32_t plates, const std::vector<uint32_t>& owners) {
    if (width < 2 || height < 2 || plates == 0 || owners.size() != static_cast<size_t>(width) * height)
        throw std::runtime_error("Invalid plate topology grid");
    width_ = width; height_ = height;
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
std::vector<uint32_t> Model::advance(const std::vector<transport::Displacement>& displacement) {
    if (!diagnostics_.enabled || displacement.size() != coverage_.size())
        throw std::runtime_error("Invalid plate topology motion");
    const size_t size = static_cast<size_t>(width_) * height_;
    for (size_t p = 0; p < coverage_.size(); ++p) {
        const std::vector<transport::Displacement> motion(size, displacement[p]);
        if (!rifts_[p].side.empty())
            rifts_[p].side = transport::remap(width_, height_, rifts_[p].side, motion);
        if (!rifts_[p].cracks.empty())
            rifts_[p].cracks = transport::remap(width_, height_, rifts_[p].cracks, motion);
        coverage_[p] = transport::remap(width_, height_, coverage_[p],
            motion);
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
        const std::vector<dynamics::Vector>& mantle, const std::vector<double>& strength) const {
    // Two connected traction domains, not random seeds or a new global Voronoi
    // map. The strongest varying mantle-flow axis selects opposite anchors.
    dynamics::Vector mean;
    uint32_t count = 0;
    for (size_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        mean.x += mantle[i].x; mean.y += mantle[i].y; ++count;
    }
    if (count == 0) return {};
    mean.x /= count; mean.y /= count;
    double xx = 0, xy = 0, yy = 0;
    for (size_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        const double x = mantle[i].x - mean.x, y = mantle[i].y - mean.y;
        xx += x*x; xy += x*y; yy += y*y;
    }
    if (xx + yy < 1e-30) return {};
    const double angle = .5 * std::atan2(2*xy, xx-yy);
    const double nx = std::cos(angle), ny = std::sin(angle);
    uint32_t low = missing, high = missing;
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (uint32_t i = 0; i < owners.size(); ++i) if (owners[i] == plate) {
        const double projection = mantle[i].x * nx + mantle[i].y * ny;
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
            const auto x=i%width_,y=i/width_;
            for (auto n:{y*width_+(x+width_-1)%width_,y*width_+(x+1)%width_,
                    ((y+height_-1)%height_)*width_+x,((y+1)%height_)*width_+x})
                if (owners[n]==plate && side[n]==missing)
                    queue.emplace(distance+2/(strength[i]+strength[n]),n,label);
        }
        return side;
    }
    std::deque<uint32_t> queue{low, high}; side[low] = 0; side[high] = 1;
    while (!queue.empty()) {
        const uint32_t i = queue.front(); queue.pop_front();
        const uint32_t x = i % width_, y = i / width_;
        for (auto n : {y*width_+(x+width_-1)%width_, y*width_+(x+1)%width_,
                ((y+height_-1)%height_)*width_+x, ((y+1)%height_)*width_+x})
            if (owners[n] == plate && side[n] == missing) { side[n] = side[i]; queue.push_back(n); }
    }
    return side;
}

std::vector<Change> Model::evolve(std::vector<uint32_t>& owners,
        const std::vector<dynamics::Vector>& mantle, const std::vector<dynamics::Vector>& velocity,
        const std::vector<double>& continent, double cell_m, double dt, double time,
        const Parameters& parameters, const std::vector<double>& strength) {
    const size_t size = static_cast<size_t>(width_) * height_, plates = coverage_.size();
    if (!std::isfinite(parameters.yield_n_m) || parameters.yield_n_m<=0 ||
        !std::isfinite(parameters.damage_time_myr) || parameters.damage_time_myr<=0 ||
        !std::isfinite(parameters.strain_rate_myr) || parameters.strain_rate_myr<0 ||
        (!strength.empty() && strength.size()!=size)) throw std::runtime_error("Invalid rift parameters");
    for (double v:strength) if (!std::isfinite(v) || v<=0) throw std::runtime_error("Invalid rock strength");
    if (!diagnostics_.enabled || owners.size() != size || mantle.size() != size || continent.size() != size ||
        velocity.size() != plates || !std::isfinite(cell_m) || cell_m <= 0 ||
        !std::isfinite(dt) || dt <= 0 || !std::isfinite(time) || time < 0)
        throw std::runtime_error("Invalid topology evolution inputs");
    for (size_t i = 0; i < size; ++i)
        if (owners[i] >= plates || !std::isfinite(mantle[i].x) || !std::isfinite(mantle[i].y) ||
            !std::isfinite(continent[i]) || continent[i] < 0 || continent[i] > 1 + 1e-8)
            throw std::runtime_error("Invalid topology forcing cell");
    for (auto v : velocity) if (!std::isfinite(v.x) || !std::isfinite(v.y))
        throw std::runtime_error("Invalid topology velocity");
    const auto neighbors = [&](uint32_t i) {
        const uint32_t x = i % width_, y = i / width_;
        return std::array<uint32_t,4>{y*width_+(x+width_-1)%width_, y*width_+(x+1)%width_,
            ((y+height_-1)%height_)*width_+x, ((y+1)%height_)*width_+x};
    };
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
    struct Contact { uint32_t faces=0, continental=0; double compression=0, normal_motion=0; };
    std::map<std::pair<uint32_t,uint32_t>,Contact> contacts;
    for (uint32_t i = 0; i < size; ++i) {
        const auto ns = neighbors(i);
        for (size_t axis = 0; axis < 2; ++axis) {
            const auto j = ns[axis == 0 ? 1 : 3], a = owners[i], b = owners[j];
            if (a == b) continue;
            auto& c = contacts[std::minmax(a,b)]; ++c.faces;
            c.continental += continent[i] >= .5 && continent[j] >= .5;
            const double closing = axis == 0 ? velocity[a].x-velocity[b].x : velocity[a].y-velocity[b].y;
            c.compression += std::max(0.0, closing);
            c.normal_motion += std::abs(closing);
        }
    }
    std::map<std::pair<uint32_t,uint32_t>,double> next_welding;
    for (const auto& entry : contacts) {
        const auto [a,b] = entry.first;
        const auto& c = entry.second;
        const auto key = std::minmax(diagnostics_.stable_ids[a], diagnostics_.stable_ids[b]);
        const double speed = std::hypot(velocity[a].x-velocity[b].x, velocity[a].y-velocity[b].y);
        const double previous = welding_.count(key) ? welding_.at(key) : 0;
        // Require broad continental contact, compatible full velocity (including
        // transform slip), and evidence of compression. Mere proximity is not welding.
        const bool locked = c.faces >= 4 && c.faces*cell_m >= 500000 &&
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
    const auto min_cells = static_cast<uint32_t>(std::max(16.0, std::ceil(minimum_plate_area_m2/(cell_m*cell_m))));
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
        dynamics::Vector flow[2]{};
        for (uint32_t i = 0; i < size; ++i) if (owners[i] == p) {
            const uint32_t s = r.side[i] > .5*coverage_[p][i] ? 1U : 0U;
            side[i] = s; ++counts[s]; anchors[s] = i;
            flow[s].x += mantle[i].x; flow[s].y += mantle[i].y;
        }
        if (std::min(counts[0],counts[1]) < min_cells || std::min(counts[0],counts[1]) < .15*count) {
            reset_candidate(r); continue;
        }
        // A throughgoing cut must leave exactly two connected resolved domains.
        // Advected candidates that develop islands are discarded, not promoted.
        std::vector<uint8_t> visited(size,0);
        bool connected = true;
        double center[2][2]{};
        for (uint32_t s = 0; s < 2; ++s) {
            std::deque<uint32_t> queue{anchors[s]}; visited[anchors[s]] = 1;
            uint32_t found = 0;
            while (!queue.empty()) {
                const auto i = queue.front(); queue.pop_front(); ++found;
                // Coordinates relative to a common anchor, unwrapped on the torus.
                center[s][0] += std::remainder(static_cast<double>(i%width_)-anchors[0]%width_, static_cast<double>(width_));
                center[s][1] += std::remainder(static_cast<double>(i/width_)-anchors[0]/width_, static_cast<double>(height_));
                for (auto n : neighbors(i)) if (!visited[n] && side[n] == s) { visited[n]=1; queue.push_back(n); }
            }
            connected &= found == counts[s];
            flow[s].x /= counts[s]; flow[s].y /= counts[s];
            center[s][0] /= counts[s]; center[s][1] /= counts[s];
        }
        if (!connected) { reset_candidate(r); continue; }
        double nx = center[1][0]-center[0][0], ny = center[1][1]-center[0][1];
        const double distance = std::hypot(nx,ny);
        if (distance < 1) { reset_candidate(r); continue; }
        nx /= distance; ny /= distance;
        uint32_t faces = 0;
        for (uint32_t i = 0; i < size; ++i) if (side[i] == 0)
            for (auto n : neighbors(i)) faces += side[n] == 1;
        const double opening = (flow[1].x-flow[0].x)*nx + (flow[1].y-flow[0].y)*ny;
        const double reduced_area = cell_m*cell_m*counts[0]*counts[1]/count;
        const double tension = faces > 0 ? dynamics::kBasalDrag*reduced_area*opening/(faces*cell_m) : 0;
        diagnostics_.max_tension_n_m=std::max(diagnostics_.max_tension_n_m,tension);
        if (parameters.progressive_rifting) {
            double total=0;uint32_t interface_cells=0;
            auto next=r.cracks;
            for (uint32_t i=0;i<size;++i) if (owners[i]==p) {
                bool interface=false;double neighbor_damage=0;
                for (auto n:neighbors(i)) if (owners[n]==p) {
                    interface|=side[n]!=side[i];
                    neighbor_damage=std::max(neighbor_damage,r.cracks[n]/std::max(coverage_[p][n],1e-12));
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

std::vector<transport::Displacement> Model::riftDisplacement(const std::vector<uint32_t>& owners,
        double dt,const Parameters& parameters) const {
    std::vector<transport::Displacement> result(owners.size());
    if (!parameters.progressive_rifting || parameters.strain_rate_myr==0) return result;
    const auto damage=damageField(owners);
    // Opening motion points away from an incipient fault. Conservative column
    // transport converts this distributed extension to thinning and subsidence;
    // no rift depth is painted onto the terrain and no rock is deleted.
    for (uint32_t i=0;i<owners.size();++i) {
        const auto p=owners[i],x=i%width_,y=i/width_;
        if (rifts_[p].side.empty() || damage[i]==0) continue;
        const bool side=rifts_[p].side[i]>.5*coverage_[p][i];
        for (const auto [j,dx,dy]:{std::tuple<uint32_t,int,int>{y*width_+(x+width_-1)%width_,-1,0},
                {y*width_+(x+1)%width_,1,0},{((y+height_-1)%height_)*width_+x,0,-1},
                {((y+1)%height_)*width_+x,0,1}}) {
            if (owners[j]!=p || (rifts_[p].side[j]>.5*coverage_[p][j])==side) continue;
            const double opening=.5*parameters.strain_rate_myr*dt*std::min(damage[i],damage[j]);
            result[i].x-=dx*opening;result[i].y-=dy*opening;
        }
    }
    return result;
}

void Model::connect(std::vector<uint32_t>& owners) {
    const auto neighbors = [&](uint32_t i) {
        const uint32_t x = i % width_, y = i / width_;
        return std::array<uint32_t, 4>{y * width_ + (x + width_ - 1) % width_,
            y * width_ + (x + 1) % width_, ((y + height_ - 1) % height_) * width_ + x,
            ((y + 1) % height_) * width_ + x};
    };
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
