#include "gtest/gtest.h"
#include "tectonic_columns.hpp"
#include "tectonic_margins.hpp"
#include "tectonic_initialization.hpp"
#include "tectonic_variation.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

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
