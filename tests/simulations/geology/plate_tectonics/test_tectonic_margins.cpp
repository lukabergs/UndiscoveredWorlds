#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_initialization.hpp"
#include "simulations/geology/plate_tectonics/tectonic_variation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
struct MarginFixture {
    static constexpr uint32_t w=64,h=16;
    std::vector<uint32_t> owners=std::vector<uint32_t>(w*h);
    std::vector<tectonics::columns::Cell> cells=std::vector<tectonics::columns::Cell>(w*h);
    std::vector<tectonics::transport::Displacement> velocity=std::vector<tectonics::transport::Displacement>(w*h);
    MarginFixture() {
        for (uint32_t i=0;i<w*h;++i) {
            const bool left=i%w<w/2;
            owners[i]=left?0:1;
            cells[i]=tectonics::columns::seed(-4000,false,left?100:20);
            velocity[i]={left?.2:-.2,0};
        }
    }
};
}
TEST(TectonicMargins, SubductionPolarityPlacesTrenchAndArcOnOppositePlates) {
    MarginFixture f;
    tectonics::margins::Parameters p; p.enabled=true;
    auto a=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,f.cells,f.velocity,p);
    EXPECT_GT(a.trench_m[31],0);
    EXPECT_DOUBLE_EQ(a.trench_m[32],0);
    EXPECT_EQ(f.owners[a.arc_receiver[31]],1U);
    EXPECT_GT(a.arc_receiver[31]%f.w,32U);
    EXPECT_LT(a.thermal_age_myr[0],100);
    EXPECT_DOUBLE_EQ(a.thermal_age_myr[31],-1);
    const auto b=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,f.cells,f.velocity,p);
    EXPECT_EQ(a.trench_m,b.trench_m);EXPECT_EQ(a.arc_receiver,b.arc_receiver);
    tectonics::margins::Polarity polarity;
    tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,f.cells,f.velocity,p,&polarity);
    for (size_t i=0;i<f.cells.size();++i) f.cells[i].oceanic_age_volume=f.cells[i].ocean*(f.owners[i]==0?10:200);
    const auto persistent=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,f.cells,f.velocity,p,&polarity);
    EXPECT_EQ(f.owners[persistent.arc_receiver[31]],1U);
    EXPECT_EQ(persistent.arc_receiver[32],persistent.arc_receiver[31]);
    EXPECT_DOUBLE_EQ(persistent.trench_m[32],0);
    for (auto& v:f.velocity) { v.y=v.x;v.x=0; }
    a=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,f.cells,f.velocity,p);
    for (size_t i=0;i<f.cells.size();++i) {
        EXPECT_DOUBLE_EQ(a.trench_m[i],0);
        EXPECT_DOUBLE_EQ(a.thermal_age_myr[i],-1);
        EXPECT_EQ(a.arc_receiver[i],std::numeric_limits<uint32_t>::max());
    }
}
TEST(TectonicMargins, ThermalHeatingPreservesGeologicAgeAndFlexureHasContinuousResponse) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;p.trench_response_myr=1.5;
    auto forcing=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,f.cells,f.velocity,p);
    tectonics::columns::Model a,b;
    a.initialize(f.w,f.h,50000,f.cells);b.initialize(f.w,f.h,50000,f.cells);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    a.advance(stopped,1,0,1,0,&forcing,p);
    for (int step=0;step<4;++step) b.advance(stopped,.25,0,1,0,&forcing,p);
    EXPECT_NEAR(a.cells()[31].flexure_m,b.cells()[31].flexure_m,1e-10);
    EXPECT_GT(a.cells()[31].flexure_m,0);
    EXPECT_LT(a.cells()[31].flexure_m,forcing.trench_m[31]);
    EXPECT_GT(a.cells()[0].thermal_rejuvenation_volume,0);
    EXPECT_NEAR(tectonics::columns::age_myr(a.cells()[0]),101,1e-10);
    EXPECT_DOUBLE_EQ(a.budget().generated_arc_kg,0);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual),1e-12);
}
TEST(TectonicMargins, ArcRockAccumulatesOnlyAfterSubductionAndCanEmerge) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    // Reference frame of a stationary overriding plate and fixed trench.
    // The older slab approaches at 4 cm/year; no island exists initially.
    for (size_t i=0;i<f.velocity.size();++i) f.velocity[i].x=f.owners[i]==0?.8:0;
    tectonics::columns::Model model;model.initialize(f.w,f.h,50000,f.cells);
    tectonics::margins::Polarity polarity;
    for (const auto& c:model.cells()) EXPECT_LT(tectonics::columns::elevation_m(c),0);
    for (int step=0;step<200;++step) {
        const auto forcing=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,model.cells(),f.velocity,p,&polarity);
        model.advance(f.velocity,1,0,1,0,&forcing,p);
    }
    EXPECT_GT(model.budget().generated_arc_kg,0);
    EXPECT_LE(model.budget().generated_arc_kg,model.budget().generated_ocean_kg);
    EXPECT_LE(model.budget().generated_arc_kg,model.budget().subducted_ocean_kg*p.arc_melt_ratio*(1+1e-12));
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-11);
    double peak=-1e10;
    for (const auto& c:model.cells()) peak=std::max(peak,tectonics::columns::elevation_m(c));
    EXPECT_GT(peak,0);
}
TEST(TectonicMargins, RiversDeliverSedimentOffshoreWithoutSmoothingQuietSeafloor) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;
    for (uint32_t i=0;i<f.w*f.h;++i) {
        const auto x=i%f.w;
        f.cells[i]=x<20?tectonics::columns::seed(2000-90*x,true,40):tectonics::columns::seed(-4000,false,40);
    }
    tectonics::columns::Model model;model.initialize(f.w,f.h,150000,f.cells);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    model.advance(stopped,1,0,1,0,nullptr,p);
    EXPECT_GT(model.budget().river_eroded_kg,0);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_GT(model.cells()[20].sediment,0);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    f.cells.assign(f.w*f.h,tectonics::columns::seed(-4000,false,40));
    model.initialize(f.w,f.h,150000,f.cells,2,1,{4000,600});
    const auto initial=model.cells();
    model.advance(stopped,1,0,0,0,nullptr,p);
    EXPECT_DOUBLE_EQ(model.budget().river_eroded_kg,0);
    for (size_t i=0;i<f.cells.size();++i) EXPECT_DOUBLE_EQ(initial[i].ocean,model.cells()[i].ocean);
}

namespace {
// The river fixture: a continental ramp 2000 - 90 x m on columns 0..19, ocean at -4000 m elsewhere.
std::vector<tectonics::columns::Cell> riverRamp(const MarginFixture& f) {
    std::vector<tectonics::columns::Cell> cells(f.w*f.h);
    for (uint32_t i=0;i<f.w*f.h;++i) {
        const auto x=i%f.w;
        cells[i]=x<20?tectonics::columns::seed(2000-90*x,true,40):tectonics::columns::seed(-4000,false,40);
    }
    return cells;
}
}

TEST(TectonicMargins, ImplicitRiversDeliverSedimentOffshoreAndCloseTheBudget) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_implicit=true;
    f.cells=riverRamp(f);
    tectonics::columns::Model model;model.initialize(f.w,f.h,150000,f.cells);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    model.advance(stopped,1,0,1,0,nullptr,p);
    EXPECT_GT(model.budget().river_eroded_kg,0);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_GT(model.cells()[20].sediment,0);
    EXPECT_EQ(model.budget().river_unsupplied_cells,0U);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    const double uniform_kg=model.budget().river_eroded_kg;
    EXPECT_NEAR(uniform_kg,3.8899922614932e16,1e6);
    // Per-cell erodibility from the process variation acts and keeps the budget closed.
    p.process_variation=.35;
    auto forcing=tectonics::margins::forcing(f.w,f.h,150000,1,f.owners,f.cells,f.velocity,p,nullptr,1,20);
    model.initialize(f.w,f.h,150000,f.cells);
    model.advance(stopped,1,0,1,0,&forcing,p);
    const double varied_kg=model.budget().river_eroded_kg;
    // Measured: 3.8895274412625e16 kg with the variation, 3.8899922614932e16
    // without (relative difference 1.19e-4: the multiplier averages near 1 on the ramp).
    EXPECT_NEAR(varied_kg,3.8895274412625e16,1e6);
    EXPECT_GT(std::abs(varied_kg-uniform_kg),1e-4*uniform_kg);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    p.process_variation=0;
    f.cells.assign(f.w*f.h,tectonics::columns::seed(-4000,false,40));
    model.initialize(f.w,f.h,150000,f.cells,2,1,{4000,600});
    const auto initial=model.cells();
    model.advance(stopped,1,0,0,0,nullptr,p);
    EXPECT_DOUBLE_EQ(model.budget().river_eroded_kg,0);
    EXPECT_DOUBLE_EQ(model.budget().deposited_kg,0);
    for (size_t i=0;i<f.cells.size();++i) EXPECT_DOUBLE_EQ(initial[i].ocean,model.cells()[i].ocean);
}

TEST(TectonicMargins, ImplicitAndExplicitRiversAgreeOnASmallStep) {
    // Same law, same receivers: with outlets at sea level the explicit step's
    // erosion is the implicit one's first-order limit as dt goes to zero.
    MarginFixture f;
    f.cells=riverRamp(f);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    const auto eroded=[&](bool implicit) {
        tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=15;
        p.river_outlet_sea_level=true;p.river_implicit=implicit;
        tectonics::columns::Model model;model.initialize(f.w,f.h,150000,f.cells);
        model.advance(stopped,.01,0,1,0,nullptr,p);
        return model.budget().river_eroded_kg;
    };
    const double explicit_kg=eroded(false),implicit_kg=eroded(true);
    EXPECT_GT(explicit_kg,0);
    // Measured: implicit / explicit = 0.99994325 at dt = 0.01 Myr (the implicit
    // step divides each drop by 1 + K sqrt(A) dt / L, at most 6.5e-4 here).
    EXPECT_LT(implicit_kg,explicit_kg);
    EXPECT_NEAR(implicit_kg/explicit_kg,1,1e-4);
}

TEST(TectonicMargins, ImplicitRiversKeepTheExplicitWavePass) {
    MarginFixture f;
    f.cells=riverRamp(f);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    tectonics::margins::Parameters p;p.enabled=true;p.river_implicit=true;
    tectonics::columns::Model rivers;rivers.initialize(f.w,f.h,150000,f.cells);
    rivers.advance(stopped,1,0,1,0,nullptr,p);
    p.wave_erosion_m_myr=50;
    tectonics::columns::Model waves;waves.initialize(f.w,f.h,150000,f.cells);
    waves.advance(stopped,1,0,1,0,nullptr,p);
    // river_eroded_kg also counts the explicit pass's wave attack.
    // Measured: 4.2952e15 kg of wave erosion beyond the rivers' 3.8900e16 kg.
    EXPECT_NEAR(waves.budget().river_eroded_kg-rivers.budget().river_eroded_kg,4.2952e15,1e11);
    EXPECT_NEAR(waves.budget().eroded_kg,waves.budget().deposited_kg,waves.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(waves.budget().relative_mass_residual),1e-12);
    p.alluvial_transport=1;
    EXPECT_THROW(waves.advance(stopped,1,0,1,0,nullptr,p),std::runtime_error);
}

TEST(TectonicMargins, ImplicitRiversPassAClosedDepressionsLoadToItsSpill) {
    // A pit at 500 m in row 8 of the ramp (neighbours 1010 to 1190 m) drains
    // through its rim to the sea; the explicit path's pit keeps what reaches it.
    MarginFixture f;
    f.cells=riverRamp(f);
    const size_t pit=8*f.w+10,shelf8=8*f.w+20,shelf3=3*f.w+20;
    f.cells[pit]=tectonics::columns::seed(500,true,40);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    const auto run=[&](bool implicit) {
        tectonics::margins::Parameters p;p.enabled=true;p.river_implicit=implicit;
        tectonics::columns::Model model;model.initialize(f.w,f.h,150000,f.cells);
        model.advance(stopped,1,0,1,0,nullptr,p);
        EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
        EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
        return model.cells();
    };
    const auto implicit=run(true),explicit_=run(false);
    // Measured: the pit's catchment (row 8 west of it and its neighbours in rows
    // 7 and 9) reaches the shelf of row 8, 82.939 m against 30.070 m on an
    // undisturbed row and 11.673 m on rows 7 and 9.
    EXPECT_DOUBLE_EQ(implicit[pit].sediment,f.cells[pit].sediment);
    EXPECT_NEAR(implicit[shelf8].sediment,82.939,1e-3);
    EXPECT_NEAR(implicit[shelf3].sediment,30.070,1e-3);
    EXPECT_GT(explicit_[pit].sediment,f.cells[pit].sediment);
}

TEST(TectonicMargins, ContinentalRiftsCreateOceanFloorWhenCrustBecomesThin) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    f.cells.assign(f.w*f.h,tectonics::columns::seed(100,true,40));
    for (size_t i=0;i<f.velocity.size();++i) f.velocity[i].x=f.owners[i]==0?-.15:.15;
    tectonics::columns::Model model;model.initialize(f.w,f.h,150000,f.cells);
    for (int step=0;step<100;++step) model.advance(f.velocity,1,0,1,0,nullptr,p);
    EXPECT_GT(model.budget().generated_ocean_kg,0);
    EXPECT_GT(model.cells()[31].ocean,0);
    EXPECT_LT(model.cells()[31].continental_area,.5);
    EXPECT_GT(tectonics::columns::elevation_m(model.cells()[31]),-6000);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-11);
}

TEST(TectonicInitialization, PreservesDetailedTerrainWithoutGloballyYoungOcean) {
    MarginFixture f;
    std::vector<double> z(f.w*f.h),age(z.size(),30);
    for (size_t i=0;i<z.size();++i) z[i]=-4200+600*std::sin(i*.7);
    const std::vector<tectonics::transport::Displacement> stopped(z.size());
    const auto cells=tectonics::initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1);
    double cooling=0;
    for (size_t i=0;i<z.size();++i) {
        EXPECT_NEAR(tectonics::columns::elevation_m(cells[i]),z[i],1e-8);
        auto older=cells[i];older.oceanic_age_volume+=older.ocean;
        cooling+=tectonics::columns::elevation_m(cells[i])-tectonics::columns::elevation_m(older);
    }
    auto young=tectonics::columns::seed(-4200,false,0);
    const double birth=tectonics::columns::elevation_m(young);
    young.oceanic_age_volume=young.ocean;
    EXPECT_LT(cooling/z.size(),.2*(birth-tectonics::columns::elevation_m(young)));
    const auto moving=tectonics::initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,f.velocity,1);
    // The periodic seam opens, the middle contact closes. Far-field relief is unchanged.
    EXPECT_GT(tectonics::columns::elevation_m(moving[0]),z[0]);
    EXPECT_NEAR(tectonics::columns::elevation_m(moving[31]),z[31],1e-8);
    EXPECT_EQ(cells[31].ocean,moving[31].ocean);
    z.assign(z.size(),-200);
    const auto shelf=tectonics::initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1);
    EXPECT_EQ(shelf[0].continental_area,1);
    EXPECT_NEAR(tectonics::columns::elevation_m(shelf[0]),-200,1e-8);
    z.assign(z.size(),1000);
    for (size_t i=0;i<z.size();++i) if (i%f.w<2) z[i]=-6000;
    const auto edge=tectonics::initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,f.velocity,1);
    for (size_t i=0;i<z.size();++i) {
        EXPECT_GE(tectonics::columns::elevation_m(edge[i]),-6000-1e-8);
        EXPECT_GT(edge[i].continent+edge[i].ocean,0);
    }
}

TEST(TectonicInitialization, GradualMarginsPreserveHeightAndAvoidCoolingDiscontinuity) {
    using namespace tectonics;
    MarginFixture f;
    std::vector<double> z(f.w*f.h),age(z.size(),40);
    for (size_t i=0;i<z.size();++i) z[i]=-1500-2000.0*(i%f.w)/(f.w-1);
    const std::vector<transport::Displacement> stopped(z.size());
    const auto cells=initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1,2000);
    for (size_t i=0;i<z.size();++i) {
        EXPECT_NEAR(columns::elevation_m(cells[i]),z[i],1e-8);
        EXPECT_NEAR(cells[i].continental_area+cells[i].oceanic_area,1,1e-14);
        EXPECT_GE(cells[i].continent,0); EXPECT_GE(cells[i].ocean,0);
        if (i%f.w>0) EXPECT_GE(cells[i].oceanic_area,cells[i-1].oceanic_area);
    }
    // A lithospheric root keeps the taper's heights: both end members carry it.
    const auto rooted=initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1,2000,1,246);
    for (size_t i=0;i<z.size();++i) EXPECT_NEAR(columns::elevation_m(rooted[i],1,246),z[i],1e-8);
    EXPECT_THROW(initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1,2000,1,500),std::runtime_error);
    const auto cooling=[&](double height) {
        z.assign(z.size(),height);
        auto c=initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1,2000)[0];
        c.oceanic_age_volume+=c.ocean*10;
        return height-columns::elevation_m(c);
    };
    EXPECT_NEAR(cooling(-2500.001),cooling(-2499.999),.01);
    EXPECT_DOUBLE_EQ(cells.front().oceanic_area,0);
    EXPECT_DOUBLE_EQ(cells[f.w-1].continental_area,0);
    EXPECT_THROW(initialization::terrain(f.w,f.h,150000,1,z,age,f.owners,stopped,1,4001),std::runtime_error);
    columns::Model model;model.initialize(f.w,f.h,150000,cells,4,1,{},true);
    model.advance(stopped,1,0,1);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
}

TEST(TectonicMargins, StrongerBendingIsLocalizedAndDoesNotConsumeStationaryRock) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    p.trench_load_n_m=3e13;p.trench_response_myr=1.5;
    auto load=tectonics::margins::forcing(f.w,f.h,156250,1,f.owners,f.cells,f.velocity,p);
    EXPECT_GT(load.trench_m[31],2000);
    EXPECT_DOUBLE_EQ(load.trench_m[20],0);
    EXPECT_DOUBLE_EQ(load.trench_m[32],0);
    tectonics::columns::Model model;model.initialize(f.w,f.h,156250,f.cells);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    for (int step=0;step<5;++step) model.advance(stopped,1,0,0,0,&load,p);
    EXPECT_LT(tectonics::columns::elevation_m(model.cells()[31],0),-4500);
    EXPECT_NEAR(tectonics::columns::elevation_m(model.cells()[20],0),-2500,1e-9);
    EXPECT_NEAR(tectonics::columns::elevation_m(model.cells()[32],0),-2500,1e-9);
    EXPECT_DOUBLE_EQ(model.budget().subducted_ocean_kg,0);
    EXPECT_NEAR(model.budget().relative_mass_residual,0,1e-12);
}

TEST(TectonicVariation, PeriodicDeterministicAndContinuousInAbsoluteTime) {
    const auto value=[](double x,double y,double t,uint32_t channel=0) {
        return tectonics::variation::multiplier(x,y,40000,20000,t,1,channel,.35,1800,80);
    };
    double low=2,high=0;
    for (int i=0;i<100;++i) {
        const double x=i*.013,y=i*.021,t=i*.9,a=value(x,y,t);
        EXPECT_DOUBLE_EQ(a,value(x,y,t));
        EXPECT_NEAR(a,value(x+1,y,t),1e-13);
        EXPECT_NEAR(a,value(x,y+1,t),1e-13);
        EXPECT_NEAR(a,value(x,y,t+80),1e-13);
        EXPECT_NEAR(a,value(x,y,t+.001),1e-4);
        EXPECT_GE(a,.65);EXPECT_LE(a,1.35);
        low=std::min(low,a);high=std::max(high,a);
    }
    EXPECT_GT(high-low,.2);
    EXPECT_NE(value(.2,.3,40,0),value(.2,.3,40,1));
    EXPECT_DOUBLE_EQ(tectonics::variation::multiplier(.2,.3,40000,20000,40,1,0,0,1800,80),1);
    EXPECT_THROW(tectonics::variation::validate(.8,1800,80),std::runtime_error);
}

TEST(TectonicVariation, ChangesProcessRatesWithoutBreakingSedimentMassBudget) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.process_variation=.35;
    for (uint32_t i=0;i<f.w*f.h;++i) {
        const auto x=i%f.w;
        f.cells[i]=x<20?tectonics::columns::seed(2000-90*x,true,40):tectonics::columns::seed(-4000,false,40);
    }
    auto forcing=tectonics::margins::forcing(f.w,f.h,150000,1,f.owners,f.cells,f.velocity,p,nullptr,1,20);
    const auto replay=tectonics::margins::forcing(f.w,f.h,150000,1,f.owners,f.cells,f.velocity,p,nullptr,1,20);
    EXPECT_EQ(forcing.trench_m,replay.trench_m);
    EXPECT_EQ(forcing.erosion_multiplier,replay.erosion_multiplier);
    EXPECT_NE(forcing.erosion_multiplier,forcing.sediment_multiplier);
    tectonics::columns::Model model;model.initialize(f.w,f.h,150000,f.cells);
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    model.advance(stopped,1,0,1,0,&forcing,p);
    EXPECT_GT(model.budget().river_eroded_kg,0);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
}

namespace {
// Earth-like continental growth keys (T4.12) on top of active margins without rivers.
tectonics::margins::Parameters growth_parameters() {
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    p.continental_min_thickness_m=20000;p.rift_breakup_thickness_m=30000;
    return p;
}
// An oceanic plate (0, left) subducting at 4 cm/yr under a stationary continental
// plate (1, right) with 500 m of sediment.
struct ContinentalMarginFixture : MarginFixture {
    ContinentalMarginFixture() {
        for (uint32_t i=0;i<w*h;++i) {
            if (owners[i]==1) { cells[i]=tectonics::columns::seed(500,true,100);cells[i].sediment=500; }
            velocity[i]={owners[i]==0?.8:0,0};
        }
    }
};
}

TEST(TectonicMargins, ContinentalMinimumThicknessCapsGhostAreaWithoutMovingRock) {
    MarginFixture f;
    for (auto& c:f.cells) {
        c=tectonics::columns::Cell{};
        c.continental_area=1;c.continent=1000;c.continental_age_volume=1000*50;
    }
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    tectonics::columns::Model legacy,capped;
    legacy.initialize(f.w,f.h,50000,f.cells);capped.initialize(f.w,f.h,50000,f.cells);
    tectonics::margins::Parameters off;off.enabled=true;off.river_erosion=0;
    const auto p=growth_parameters();
    legacy.advance(stopped,1,0,1,0,nullptr,off);
    capped.advance(stopped,1,0,1,0,nullptr,p);
    const double freed=(1-.05)*-std::expm1(-1.0/5);
    for (size_t i=0;i<f.cells.size();++i) {
        EXPECT_DOUBLE_EQ(legacy.cells()[i].continental_area,1);
        EXPECT_NEAR(capped.cells()[i].continental_area,1-freed,1e-12);
        EXPECT_NEAR(capped.cells()[i].oceanic_area,freed,1e-12);
        EXPECT_NEAR(capped.cells()[i].continent,1000,1e-9);
        EXPECT_DOUBLE_EQ(capped.cells()[i].ocean,0);
    }
    EXPECT_NEAR(capped.budget().oceanized_continental_area_m2,freed*50000.0*50000*f.cells.size(),1e-3*50000.0*50000);
    EXPECT_DOUBLE_EQ(legacy.budget().oceanized_continental_area_m2,0);
    for (int step=1;step<100;++step) capped.advance(stopped,1,0,1,0,nullptr,p);
    for (const auto& c:capped.cells()) {
        EXPECT_NEAR(c.continental_area,.05,1e-6);
        EXPECT_NEAR(c.continent/c.continental_area,20000,1);
    }
    EXPECT_LT(std::abs(capped.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(capped.budget().felsic_relative_residual),1e-12);
}

TEST(TectonicMargins, ContinentalMinimumThicknessMakesOceanOnceRiftedCrustReachesIt) {
    MarginFixture f;
    f.cells.assign(f.w*f.h,tectonics::columns::seed(100,true,40));
    for (size_t i=0;i<f.velocity.size();++i) f.velocity[i].x=f.owners[i]==0?-.15:.15;
    tectonics::margins::Parameters legacy_parameters;legacy_parameters.enabled=true;legacy_parameters.river_erosion=0;
    // The legacy share is the minimum thickness 0 written out explicitly.
    auto explicit_zero=legacy_parameters;explicit_zero.continental_min_thickness_m=0;
    const auto p=growth_parameters();
    tectonics::columns::Model legacy,zero,limited;
    for (auto* model:{&legacy,&zero,&limited}) model->initialize(f.w,f.h,150000,f.cells);
    const auto thinnest=[&](const tectonics::columns::Model& model) {
        double least=std::numeric_limits<double>::infinity();
        for (const auto& c:model.cells()) if (c.continental_area>.05) least=std::min(least,c.continent/c.continental_area);
        return least;
    };
    double limited_least=std::numeric_limits<double>::infinity();
    for (int step=0;step<100;++step) {
        legacy.advance(f.velocity,1,0,1,0,nullptr,legacy_parameters);
        zero.advance(f.velocity,1,0,1,0,nullptr,explicit_zero);
        limited.advance(f.velocity,1,0,1,0,nullptr,p);
        limited_least=std::min(limited_least,thinnest(limited));
    }
    for (size_t i=0;i<f.cells.size();++i) {
        EXPECT_EQ(legacy.cells()[i].continent,zero.cells()[i].continent);
        EXPECT_EQ(legacy.cells()[i].continental_area,zero.cells()[i].continental_area);
    }
    // Legacy rifts stretch the crust far below 20 km before it breaks up; with the
    // minimum, the opening is ocean once the crust reaches it; the cap relaxes the
    // overshoot of a fast opening (0.3 cells/Myr here) within its 5 Myr.
    EXPECT_LT(thinnest(legacy),15000);
    EXPECT_GT(limited_least,.7*p.continental_min_thickness_m);
    EXPECT_GT(thinnest(limited),.95*p.continental_min_thickness_m);
    EXPECT_GT(limited.budget().generated_ocean_kg,legacy.budget().generated_ocean_kg);
    EXPECT_LT(std::abs(limited.budget().relative_mass_residual),1e-11);
    EXPECT_LT(std::abs(limited.budget().felsic_relative_residual),1e-11);
}

TEST(TectonicMargins, ArcContinentShareKeepsMassAndBuildsIslandArcArea) {
    MarginFixture f;
    auto p=growth_parameters();p.arc_continent_fraction=.5;
    for (size_t i=0;i<f.velocity.size();++i) f.velocity[i].x=f.owners[i]==0?.8:0;
    tectonics::columns::Model model;model.initialize(f.w,f.h,50000,f.cells);
    tectonics::margins::Polarity polarity;
    for (int step=0;step<200;++step) {
        const auto forcing=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,model.cells(),f.velocity,p,&polarity);
        model.advance(f.velocity,1,0,1,0,&forcing,p);
    }
    const auto& b=model.budget();
    EXPECT_GT(b.generated_continent_kg,0);
    EXPECT_NEAR(b.generated_continent_kg,.5*b.generated_arc_kg,b.generated_arc_kg*1e-12);
    EXPECT_DOUBLE_EQ(b.arc_continent_fraction,.5);
    EXPECT_DOUBLE_EQ(b.initial_felsic_kg,0);
    EXPECT_LT(std::abs(b.relative_mass_residual),1e-11);
    EXPECT_NEAR(b.felsic_kg,b.generated_continent_kg,b.generated_continent_kg*1e-11);
    bool island=false;
    for (size_t i=0;i<f.cells.size();++i)
        island=island || (f.owners[i]==1 && model.cells()[i].continental_area>0 && model.cells()[i].continent>0);
    EXPECT_TRUE(island);
    auto invalid=p;invalid.continental_min_thickness_m=0;
    EXPECT_THROW(tectonics::margins::validate(invalid),std::runtime_error);
}

TEST(TectonicMargins, SubductionErosionAndSedimentSubductionAreReservoirSinks) {
    ContinentalMarginFixture f;
    auto p=growth_parameters();p.sediment_subduction_fraction=.75;p.subduction_erosion_ratio=.07;
    tectonics::columns::Model model;model.initialize(f.w,f.h,50000,f.cells);
    tectonics::margins::Polarity polarity;
    for (int step=0;step<100;++step) {
        const auto forcing=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,model.cells(),f.velocity,p,&polarity);
        model.advance(f.velocity,1,0,1,0,&forcing,p);
    }
    const auto& b=model.budget();
    EXPECT_GT(b.subducted_continent_kg,0);
    EXPECT_GT(b.subducted_sediment_kg,0);
    EXPECT_LE(b.subducted_continent_kg,.07*b.subducted_ocean_kg*2800/2900*(1+1e-12));
    EXPECT_DOUBLE_EQ(b.generated_continent_kg,0);
    EXPECT_LT(b.felsic_kg,b.initial_felsic_kg);
    EXPECT_LT(std::abs(b.relative_mass_residual),1e-11);
    EXPECT_LT(std::abs(b.felsic_relative_residual),1e-11);
}

TEST(TectonicMargins, GrowthControllerHoldsTheReservoirTargetWithAShareInZeroToOne) {
    const auto run=[](double growth,int steps) {
        ContinentalMarginFixture f;
        auto p=growth_parameters();p.sediment_subduction_fraction=.75;p.subduction_erosion_ratio=.07;
        p.continental_growth_control=true;p.continental_growth_per_gyr=growth;
        tectonics::columns::Model model;model.initialize(f.w,f.h,50000,f.cells);
        tectonics::margins::Polarity polarity;
        for (int step=0;step<steps;++step) {
            const auto forcing=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,model.cells(),f.velocity,p,&polarity);
            model.advance(f.velocity,1,0,1,0,&forcing,p);
            EXPECT_GE(model.budget().arc_continent_fraction,0);
            EXPECT_LE(model.budget().arc_continent_fraction,1);
        }
        EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-11);
        EXPECT_LT(std::abs(model.budget().felsic_relative_residual),1e-11);
        return model.budget();
    };
    // Steady reservoir: the arcs replace exactly what the trench takes.
    const auto steady=run(0,150);
    EXPECT_GT(steady.subducted_continent_kg+steady.subducted_sediment_kg,0);
    EXPECT_GT(steady.generated_continent_kg,0);
    EXPECT_DOUBLE_EQ(steady.continental_growth_shortfall_kg,0);
    EXPECT_NEAR(steady.felsic_kg,steady.initial_felsic_kg,steady.initial_felsic_kg*1e-9);
    // Growth of 0.2 per Gyr over 150 Myr: 3 % more.
    const auto growing=run(.2,150);
    EXPECT_DOUBLE_EQ(growing.continental_growth_shortfall_kg,0);
    EXPECT_NEAR(growing.felsic_kg,growing.initial_felsic_kg*1.03,growing.initial_felsic_kg*1e-9);
    // A rate the arcs cannot supply: the share saturates at 1 and the deficit is booked.
    const auto starved=run(50,150);
    EXPECT_DOUBLE_EQ(starved.arc_continent_fraction,1);
    EXPECT_GT(starved.continental_growth_shortfall_kg,0);
    EXPECT_LT(starved.felsic_kg,starved.initial_felsic_kg*(1+50*.15));
    // Only the physical demand is booked (the sinks and the target's growth), never a lasting deficit again.
    EXPECT_LE(starved.continental_growth_shortfall_kg,(starved.subducted_continent_kg+starved.subducted_sediment_kg
        +starved.initial_felsic_kg*50*.15)*(1+1e-12));
    // Feedback: 100 Myr of sinks without arcs leave the reservoir below its target; the
    // controller then closes the gap over its 50 Myr response, booking no shortfall.
    ContinentalMarginFixture f;
    auto p=growth_parameters();p.sediment_subduction_fraction=.75;p.subduction_erosion_ratio=.07;
    tectonics::columns::Model model;model.initialize(f.w,f.h,50000,f.cells);
    tectonics::margins::Polarity polarity;
    const auto step=[&] {
        const auto forcing=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,model.cells(),f.velocity,p,&polarity);
        model.advance(f.velocity,1,0,1,0,&forcing,p);
    };
    for (int i=0;i<100;++i) step();
    const double initial=model.budget().initial_felsic_kg;
    const double deficit=initial-model.budget().felsic_kg;
    ASSERT_GT(deficit,1e-3*initial);
    p.continental_growth_control=true;
    double error=deficit;
    for (int i=0;i<150;++i) {
        step();
        EXPECT_GE(model.budget().arc_continent_fraction,0);
        EXPECT_LE(model.budget().arc_continent_fraction,1);
        const double next=initial-model.budget().felsic_kg;
        EXPECT_LE(std::abs(next),error*(1+1e-9));
        error=std::abs(next);
    }
    EXPECT_LT(error,.25*deficit);
    EXPECT_DOUBLE_EQ(model.budget().continental_growth_shortfall_kg,0);
    EXPECT_LT(std::abs(model.budget().felsic_relative_residual),1e-11);
}

namespace {
// Uniform mixed columns: no motion, no erosion and no root flow between equal cells.
std::vector<tectonics::columns::Cell> mixed_cells(size_t count,double area,double continent,double ocean) {
    tectonics::columns::Cell c;
    c.continental_area=area;c.oceanic_area=1-area;c.continent=continent;c.ocean=ocean;
    c.continental_age_volume=continent*100;c.oceanic_age_volume=ocean*20;
    return std::vector<tectonics::columns::Cell>(count,c);
}
}

TEST(TectonicMargins, ContinentalMaximumThicknessGainsAreaWithoutMovingRock) {
    MarginFixture f;
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    const auto cells=mixed_cells(f.cells.size(),.5,30000,3500);
    tectonics::columns::Model reference,capped;
    reference.initialize(f.w,f.h,50000,cells);capped.initialize(f.w,f.h,50000,cells);
    auto p=growth_parameters();p.continental_max_thickness_m=40000;
    reference.advance(stopped,1,0,1,0,nullptr,growth_parameters());
    capped.advance(stopped,1,0,1,0,nullptr,p);
    const double gained=.25*-std::expm1(-1.0/10);
    double sum=0;
    for (size_t i=0;i<f.cells.size();++i) {
        const auto& c=capped.cells()[i];
        EXPECT_DOUBLE_EQ(reference.cells()[i].continental_area,.5);
        EXPECT_NEAR(c.continental_area,.5+gained,1e-12);
        EXPECT_NEAR(c.oceanic_area,.5-gained,1e-12);
        EXPECT_EQ(c.continent,reference.cells()[i].continent);
        EXPECT_EQ(c.ocean,reference.cells()[i].ocean);
        EXPECT_EQ(c.continental_age_volume,reference.cells()[i].continental_age_volume);
        EXPECT_EQ(c.oceanic_age_volume,reference.cells()[i].oceanic_age_volume);
        sum+=c.continental_area-.5;
    }
    const double cell_area=50000.0*50000;
    EXPECT_NEAR(capped.budget().continentalized_oceanic_area_m2,sum*cell_area,1e-12*sum*cell_area);
    EXPECT_DOUBLE_EQ(reference.budget().continentalized_oceanic_area_m2,0);
    for (int step=1;step<200;++step) capped.advance(stopped,1,0,1,0,nullptr,p);
    for (const auto& c:capped.cells()) {
        EXPECT_NEAR(c.continental_area,.75,1e-6);
        EXPECT_NEAR(c.continent/c.continental_area,40000,1);
    }
    EXPECT_LT(std::abs(capped.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(capped.budget().felsic_relative_residual),1e-12);
    // Room for more area than the oceanic area: the column saturates at area 1.
    tectonics::columns::Model saturated;
    saturated.initialize(f.w,f.h,50000,mixed_cells(f.cells.size(),.9,50000,350));
    for (int step=0;step<200;++step) saturated.advance(stopped,1,0,1,0,nullptr,p);
    // The cap moves no ocean; as the oceanic area vanishes, the remap's rounding-sized
    // overlap is a large share of it and subducts (budgeted) a little of the ocean volume.
    double lost=0;
    for (const auto& c:saturated.cells()) {
        EXPECT_NEAR(c.continental_area,1,1e-6);
        EXPECT_GE(c.oceanic_area,0);
        EXPECT_LE(c.ocean,350);
        EXPECT_GT(c.ocean,349);
        lost+=350-c.ocean;
    }
    EXPECT_NEAR(saturated.budget().subducted_ocean_kg-saturated.budget().generated_ocean_kg,lost*2900*cell_area,1e-6*350*2900*cell_area);
    EXPECT_LT(std::abs(saturated.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(saturated.budget().felsic_relative_residual),1e-12);
}

TEST(TectonicMargins, ContinentalThicknessBetweenTheCapsKeepsItsArea) {
    MarginFixture f;
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    tectonics::columns::Model model;
    model.initialize(f.w,f.h,50000,mixed_cells(f.cells.size(),.5,15000,3500));
    auto p=growth_parameters();p.continental_max_thickness_m=40000;p.orogenic_collapse_km2_myr=3000;
    for (int step=0;step<50;++step) model.advance(stopped,1,0,1,0,nullptr,p);
    for (const auto& c:model.cells()) EXPECT_DOUBLE_EQ(c.continental_area,.5);
    EXPECT_DOUBLE_EQ(model.budget().oceanized_continental_area_m2,0);
    EXPECT_DOUBLE_EQ(model.budget().continentalized_oceanic_area_m2,0);
    EXPECT_DOUBLE_EQ(model.budget().collapsed_continent_kg,0);
}

TEST(TectonicMargins, OrogenicCollapseSpreadsThickCrustIntoMixedAndOceanCells) {
    MarginFixture f;
    // A 60 km strip (columns 20 to 23, area 1), a mixed column (area 0.3, 10 km of
    // cell-mean continent, 33 km thick) east of it, ocean elsewhere.
    for (uint32_t i=0;i<f.w*f.h;++i) {
        const auto x=i%f.w;
        if (x>=20 && x<24) f.cells[i]=mixed_cells(1,1,60000,0)[0];
        else if (x==24) f.cells[i]=mixed_cells(1,.3,10000,4900)[0];
    }
    const std::vector<tectonics::transport::Displacement> stopped(f.cells.size());
    auto p=growth_parameters();p.continental_max_thickness_m=40000;p.orogenic_collapse_km2_myr=3000;
    p.root_mobility=0; // Collapse alone; root flow never moves volume into pure ocean.
    tectonics::columns::Model model;model.initialize(f.w,f.h,50000,f.cells);
    const auto sums=[&] {
        double area=0,continent=0;
        for (const auto& c:model.cells()) { area+=c.continental_area;continent+=c.continent; }
        return std::pair<double,double>{area,continent};
    };
    const auto [area0,continent0]=sums();
    for (int step=0;step<50;++step) {
        model.advance(stopped,1,0,1,0,nullptr,p);
        for (size_t i=0;i<f.cells.size();++i) if (i%f.w>=20 && i%f.w<24)
            EXPECT_GE(model.cells()[i].continent,p.continental_max_thickness_m*(1-1e-12));
    }
    const auto [area1,continent1]=sums();
    EXPECT_GT(area1,area0);
    EXPECT_NEAR(continent1,continent0,1e-12*continent0);
    EXPECT_LT(model.cells()[20].continent,60000);
    EXPECT_LT(model.cells()[23].continent,60000);
    // The ocean west of the strip gained continent of the strip's age.
    const auto& west=model.cells()[19],&source=model.cells()[21];
    ASSERT_GT(west.continent,0);
    EXPECT_GT(west.continental_area,0);
    EXPECT_NEAR(west.continental_age_volume/west.continent,source.continental_age_volume/source.continent,1e-9);
    for (const auto& c:model.cells()) EXPECT_FALSE(c.continental_area>=.5 && c.continent<1000);
    EXPECT_GT(model.budget().collapsed_continent_kg,0);
    EXPECT_GT(model.budget().continentalized_oceanic_area_m2,0);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(model.budget().felsic_relative_residual),1e-12);
}

TEST(TectonicMargins, ZeroMaximumThicknessAndCollapseAreTheLegacyEngine) {
    ContinentalMarginFixture f;
    auto p=growth_parameters();p.sediment_subduction_fraction=.75;p.subduction_erosion_ratio=.07;
    p.continental_growth_control=true;
    auto explicit_zero=p;explicit_zero.continental_max_thickness_m=0;explicit_zero.orogenic_collapse_km2_myr=0;
    tectonics::columns::Model legacy,zero;
    legacy.initialize(f.w,f.h,50000,f.cells);zero.initialize(f.w,f.h,50000,f.cells);
    tectonics::margins::Polarity legacy_polarity,zero_polarity;
    for (int step=0;step<100;++step) {
        const auto a=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,legacy.cells(),f.velocity,p,&legacy_polarity);
        const auto b=tectonics::margins::forcing(f.w,f.h,50000,1,f.owners,zero.cells(),f.velocity,explicit_zero,&zero_polarity);
        legacy.advance(f.velocity,1,0,1,0,&a,p);
        zero.advance(f.velocity,1,0,1,0,&b,explicit_zero);
    }
    for (size_t i=0;i<f.cells.size();++i) {
        const auto& a=legacy.cells()[i];const auto& b=zero.cells()[i];
        EXPECT_EQ(a.continent,b.continent);EXPECT_EQ(a.ocean,b.ocean);EXPECT_EQ(a.sediment,b.sediment);
        EXPECT_EQ(a.continental_area,b.continental_area);EXPECT_EQ(a.oceanic_area,b.oceanic_area);
        EXPECT_EQ(a.continental_age_volume,b.continental_age_volume);EXPECT_EQ(a.oceanic_age_volume,b.oceanic_age_volume);
        EXPECT_EQ(a.thermal_rejuvenation_volume,b.thermal_rejuvenation_volume);EXPECT_EQ(a.flexure_m,b.flexure_m);
    }
    EXPECT_EQ(legacy.budget().mass_kg,zero.budget().mass_kg);
    EXPECT_EQ(legacy.budget().felsic_kg,zero.budget().felsic_kg);
    EXPECT_DOUBLE_EQ(zero.budget().continentalized_oceanic_area_m2,0);
    EXPECT_DOUBLE_EQ(zero.budget().collapsed_continent_kg,0);
}

TEST(TectonicMargins, MaximumThicknessAndCollapseAreValidated) {
    auto p=growth_parameters();p.continental_max_thickness_m=40000;p.orogenic_collapse_km2_myr=3000;
    EXPECT_NO_THROW(tectonics::margins::validate(p));
    for (double h_max : {30000.0,25000.0,-1.0,std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid=p;invalid.continental_max_thickness_m=h_max;
        EXPECT_THROW(tectonics::margins::validate(invalid),std::runtime_error) << h_max;
    }
    for (double collapse : {-1.0,std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid=p;invalid.orogenic_collapse_km2_myr=collapse;
        EXPECT_THROW(tectonics::margins::validate(invalid),std::runtime_error) << collapse;
    }
    auto no_maximum=p;no_maximum.continental_max_thickness_m=0;
    EXPECT_THROW(tectonics::margins::validate(no_maximum),std::runtime_error);
    tectonics::margins::Parameters margins_off;margins_off.continental_max_thickness_m=40000;
    EXPECT_THROW(tectonics::margins::validate(margins_off),std::runtime_error);
    margins_off.enabled=true;
    EXPECT_NO_THROW(tectonics::margins::validate(margins_off));
}
