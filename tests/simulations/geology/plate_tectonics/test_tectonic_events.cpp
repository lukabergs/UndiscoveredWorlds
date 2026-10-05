#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_topology.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include <algorithm>
#include <deque>
#include <set>

namespace {
using namespace tectonics::topology;
constexpr uint32_t w = 64, h = 32;
std::vector<uint32_t> island() {
    std::vector<uint32_t> owners(w*h,1);
    for (uint32_t y=4; y<28; ++y) for (uint32_t x=8; x<40; ++x) owners[y*w+x]=0;
    return owners;
}
void check_connected(const std::vector<uint32_t>& owners) {
    std::vector<bool> visited(owners.size()); std::set<uint32_t> seen;
    for (uint32_t i=0; i<owners.size(); ++i) if (!visited[i]) {
        EXPECT_TRUE(seen.insert(owners[i]).second) << "detached owner " << owners[i];
        std::deque<uint32_t> queue{i}; visited[i]=true;
        while (!queue.empty()) {
            const auto a=queue.front(); queue.pop_front(); const auto x=a%w,y=a/w;
            for (auto b : {y*w+(x+w-1)%w,y*w+(x+1)%w,((y+h-1)%h)*w+x,((y+1)%h)*w+x})
                if (!visited[b] && owners[b]==owners[a]) { visited[b]=true; queue.push_back(b); }
        }
    }
}
}

TEST(TectonicEvents, SustainedExtensionSplitsConnectedDomainsAndRecordsAncestry) {
    auto owners=island(), initial=owners;
    Model model; model.initialize(w,h,2,owners);
    std::vector<tectonics::dynamics::Vector> mantle(w*h);
    for (uint32_t i=0; i<w*h; ++i) if (owners[i]==0) mantle[i].x=(static_cast<double>(i%w)-23.5)*2e-10;
    std::vector<Change> changes;
    for (int t=1; t<=100 && changes.empty(); ++t) {
        owners=model.advance({{}, {}});
        changes=model.evolve(owners,mantle,{{},{}},std::vector<double>(w*h,1),250000,1,t);
        if (t<20) EXPECT_TRUE(changes.empty());
    }
    ASSERT_EQ(changes.size(),1U); EXPECT_EQ(changes[0].kind,EventKind::Split);
    ASSERT_EQ(model.diagnostics().events.size(),1U);
    const auto& event=model.diagnostics().events[0];
    EXPECT_EQ(event.parent_a,0U); EXPECT_EQ(event.result,2U); EXPECT_GE(event.time_myr,20);
    EXPECT_EQ(model.diagnostics().stable_ids,(std::vector<uint32_t>{0,1,2}));
    for (size_t i=0;i<owners.size();++i) if (initial[i]==1) EXPECT_EQ(owners[i],1U);
    check_connected(owners);
    // Once born, daughters are independent bodies with continuous coverage.
    owners=model.advance({{-.1,0},{},{.1,0}}); check_connected(owners);
    EXPECT_EQ(model.diagnostics().events.size(),1U);
}

TEST(TectonicEvents, ProgressiveDamagePrecedesConnectedSplitAndCreatesExtension) {
    auto owners=island();Model model;model.initialize(w,h,2,owners);
    std::vector<tectonics::dynamics::Vector> mantle(w*h);
    for (uint32_t i=0;i<w*h;++i) if (owners[i]==0)
        mantle[i].x=(static_cast<double>(i%w)-23.5)*2e-10;
    Parameters p;p.progressive_rifting=true;p.yield_n_m=5e11;p.damage_time_myr=30;
    std::vector<double> strength(w*h,1);
    for (uint32_t i=0;i<w*h;++i) if (i%w>=22 && i%w<=25) strength[i]=.5;
    bool damaged=false,extending=false;std::vector<Change> changes;
    for (int t=1;t<=200 && changes.empty();++t) {
        owners=model.advance({{}, {}});
        changes=model.evolve(owners,mantle,{{},{}},std::vector<double>(w*h,1),250000,1,t,p,strength);
        const auto damage=model.damageField(owners);
        damaged|=*std::max_element(damage.begin(),damage.end())>0;
        for (auto v:model.riftDisplacement(owners,1,p)) extending|=std::hypot(v.x,v.y)>0;
    }
    EXPECT_TRUE(damaged);EXPECT_TRUE(extending);
    ASSERT_EQ(changes.size(),1U);EXPECT_EQ(changes[0].kind,EventKind::Split);
    check_connected(owners);
}

TEST(TectonicEvents, CompressionTransformShearAndUniformFlowDoNotSplit) {
    for (int kind=0; kind<3; ++kind) {
        auto owners=island(); Model model; model.initialize(w,h,2,owners);
        std::vector<tectonics::dynamics::Vector> mantle(w*h);
        for (uint32_t i=0; i<w*h; ++i) if (owners[i]==0) {
            const double delta=(static_cast<double>(i%w)-23.5)*2e-10;
            mantle[i]=kind==0 ? tectonics::dynamics::Vector{-delta,0}
                : kind==1 ? tectonics::dynamics::Vector{0,delta} : tectonics::dynamics::Vector{1e-9,1e-9};
        }
        for (int t=1;t<=100;++t)
            EXPECT_TRUE(model.evolve(owners,mantle,{{},{}},std::vector<double>(w*h,1),250000,1,t).empty());
        EXPECT_TRUE(model.diagnostics().events.empty());
    }
}

TEST(TectonicEvents, FusionRequiresPersistentContinentalLockingAndPreservesIdentityAfterCompaction) {
    std::vector<uint32_t> owners(w*h);
    for (uint32_t i=0;i<w*h;++i) owners[i]=i%w<w/2?0:1;
    Model model; model.initialize(w,h,2,owners);
    const std::vector<tectonics::dynamics::Vector> velocity{{1e-11,0},{-1e-11,0}}, mantle(w*h);
    const std::vector<double> continent(w*h,1);
    for (int t=1;t<30;++t) EXPECT_TRUE(model.evolve(owners,mantle,velocity,continent,250000,1,t).empty());
    const auto change=model.evolve(owners,mantle,velocity,continent,250000,1,30);
    ASSERT_EQ(change.size(),1U); EXPECT_EQ(change[0].kind,EventKind::Fusion);
    check_connected(owners);
    EXPECT_EQ(*std::max_element(owners.begin(),owners.end()),0U);
    ASSERT_EQ(model.diagnostics().events.size(),1U);
    EXPECT_EQ(model.diagnostics().events[0].parent_a,0U);
    EXPECT_EQ(model.diagnostics().events[0].parent_b,1U);
    EXPECT_EQ(model.diagnostics().events[0].result,2U);
    model.removePlate(1);
    EXPECT_EQ(model.diagnostics().stable_ids,(std::vector<uint32_t>{2}));
    EXPECT_EQ(model.advance({{.1,-.1}}),owners);
}

TEST(TectonicEvents, OceanContactsFastSlipAndStationaryNeighborsDoNotFuse) {
    for (int kind=0;kind<4;++kind) {
        std::vector<uint32_t> owners(w*h);
        for (uint32_t i=0;i<w*h;++i) owners[i]=i%w<w/2?0:1;
        Model model; model.initialize(w,h,2,owners);
        const tectonics::dynamics::Vector v=kind==0 ? tectonics::dynamics::Vector{1e-11,0}
            : kind==1 ? tectonics::dynamics::Vector{1e-9,0}
            : kind==2 ? tectonics::dynamics::Vector{0,1e-11} : tectonics::dynamics::Vector{};
        for (int t=1;t<=100;++t)
            EXPECT_TRUE(model.evolve(owners,std::vector<tectonics::dynamics::Vector>(w*h),{v,{}},
                std::vector<double>(w*h,kind==0?0:1),250000,1,t).empty());
    }
}

TEST(TectonicEvents, ShortLockingEpisodeDoesNotCarryOverAcrossRenewedSlip) {
    std::vector<uint32_t> owners(w*h);
    for (uint32_t i=0;i<w*h;++i) owners[i]=i%w<w/2?0:1;
    Model model; model.initialize(w,h,2,owners);
    for (int t=1;t<=50;++t) {
        const double speed=t==25?1e-9:1e-11;
        EXPECT_TRUE(model.evolve(owners,std::vector<tectonics::dynamics::Vector>(w*h),{{speed,0},{}},
            std::vector<double>(w*h,1),250000,1,t).empty());
    }
}

TEST(TectonicEvents, TimestepRefinementPreservesSustainedFailureTime) {
    double times[2]{};
    for (int k=0;k<2;++k) {
        auto owners=island(); Model model; model.initialize(w,h,2,owners);
        std::vector<tectonics::dynamics::Vector> mantle(w*h);
        for (uint32_t i=0;i<w*h;++i) if (owners[i]==0) mantle[i].x=(static_cast<double>(i%w)-23.5)*2e-10;
        const double dt=k==0?1:.25;
        for (double t=dt;t<=100;t+=dt)
            if (!model.evolve(owners,mantle,{{},{}},std::vector<double>(w*h,1),250000,dt,t).empty()) { times[k]=t; break; }
        EXPECT_GT(times[k],0);
    }
    EXPECT_NEAR(times[0],times[1],1);
}

TEST(TectonicEvents, RiftFailureIsUnchangedWhenTheSamePlateCrossesBothPeriodicSeams) {
    double times[2]{}; uint32_t areas[2]{};
    const auto initial=island();
    for (int k=0;k<2;++k) {
        std::vector<uint32_t> owners(w*h);
        std::vector<tectonics::dynamics::Vector> mantle(w*h);
        for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x) {
            const uint32_t target=((y+14*k)%h)*w+(x+40*k)%w;
            owners[target]=initial[y*w+x];
            if (owners[target]==0) mantle[target].x=(static_cast<double>(x)-23.5)*2e-10;
        }
        Model model; model.initialize(w,h,2,owners);
        for (int t=1;t<=100;++t)
            if (!model.evolve(owners,mantle,{{},{}},std::vector<double>(w*h,1),250000,1,t).empty()) {
                times[k]=t; areas[k]=model.diagnostics().events[0].cells; break;
            }
        EXPECT_GT(times[k],0); check_connected(owners);
    }
    EXPECT_DOUBLE_EQ(times[0],times[1]); EXPECT_EQ(areas[0],areas[1]);
}
