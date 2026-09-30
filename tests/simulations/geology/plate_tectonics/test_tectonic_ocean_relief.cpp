#include "gtest/gtest.h"
#include "tectonic_columns.hpp"
#include "tectonic_ocean_relief.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

TEST(TectonicOceanRelief, SmoothPeriodicDeterministicNoiseAndPositiveCrust) {
    using namespace tectonics::ocean_relief;
    for (int i=-20;i<20;++i) {
        const double u=i*.137,v=i*.073;
        EXPECT_NEAR(noise(u,v,7,5,1),noise(u+1,v-1,7,5,1),1e-13);
        EXPECT_LE(std::abs(noise(u,v,7,5,1)),1);
    }
    EXPECT_NEAR(noise(1e-7,.23,7,5,1),noise(-1e-7,.23,7,5,1),1e-9);
    const auto a=crust(64,32,50000,1,{6000,400});
    EXPECT_EQ(a,crust(64,32,50000,1,{6000,400}));
    EXPECT_NE(a,crust(64,32,50000,2,{6000,400}));
    EXPECT_GT(*std::min_element(a.begin(),a.end()),0);
    for (double v : crust(64,32,50000,1,{})) EXPECT_DOUBLE_EQ(v,7000);
    EXPECT_THROW(validate({7000,400}),std::runtime_error);
    EXPECT_THROW(validate({100,-1}),std::runtime_error);
}

TEST(TectonicOceanRelief, ReliefIsMaterialAndDoesNotRegenerateOrDependOnUpdateCount) {
    using namespace tectonics;
    constexpr uint32_t w=32,h=16;
    std::vector<columns::Cell> cells(w*h,columns::seed(-4000,false,40));
    columns::Model a,b;
    a.initialize(w,h,150000,cells,2,1,{5000,600});
    b.initialize(w,h,150000,cells,2,1,{5000,600});
    const auto initial=a.cells();
    const std::vector<columns::Displacement> stopped(w*h);
    a.advance(stopped,1,0,0);
    for (int i=0;i<4;++i) b.advance(stopped,.25,0,0);
    for (size_t i=0;i<cells.size();++i) {
        EXPECT_DOUBLE_EQ(a.cells()[i].ocean,initial[i].ocean);
        EXPECT_DOUBLE_EQ(a.cells()[i].ocean,b.cells()[i].ocean);
        EXPECT_NEAR(columns::age_myr(a.cells()[i]),41,1e-12);
    }
    for (int i=0;i<20;++i) a.advance(std::vector<columns::Displacement>(w*h,{.25,0}),1,0,0);
    double error=0,signal=0;
    for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x) {
        const double expected=initial[y*w+(x+w-5)%w].ocean;
        error+=std::pow(a.cells()[y*w+x].ocean-expected,2);
        signal+=std::pow(expected-7000,2);
    }
    EXPECT_LT(error/signal,.05);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual),1e-12);
}

TEST(TectonicOceanRelief, SparseBasaltSeamountsFormIslandsAndRemainInTheMassBudget) {
    using namespace tectonics;
    constexpr uint32_t w=128,h=64;
    const ocean_relief::Parameters relief{0,600,45000,350,3500};
    const auto rocks=ocean_relief::crust(w,h,100000,1,relief);
    EXPECT_EQ(rocks,ocean_relief::crust(w,h,100000,1,relief));
    std::vector<columns::Cell> initial(w*h,columns::seed(-4000,false,40));
    columns::Model model; model.initialize(w,h,100000,initial,2,1,relief);
    size_t islands=0; double mass=0;
    for (const auto& c : model.cells()) {
        islands+=columns::elevation_m(c)>0;
        EXPECT_NEAR(columns::age_myr(c),40,1e-12);
        mass+=columns::mass_per_area(c)*1e10;
    }
    EXPECT_GT(islands,10U); EXPECT_LT(islands,initial.size()/10);
    EXPECT_NEAR(mass,model.budget().initial_mass_kg,mass*1e-13);
    for (int t=0;t<20;++t) model.advance(std::vector<columns::Displacement>(w*h,{.1,-.1}),1,1,1);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
}

TEST(TectonicOceanRelief, FineTransportSubstepsAgeRockOnlyOnceAndCloseMass) {
    using namespace tectonics;
    constexpr uint32_t w=32,h=16;
    std::vector<columns::Cell> cells(w*h,columns::seed(-4000,false,40));
    columns::Model model;
    model.initialize(w,h,150000,cells,4,1,{5000,600});
    model.advance(std::vector<columns::Displacement>(w*h,{.9,-.8}),2.5,0,0);
    for (const auto& c : model.cells()) {
        EXPECT_NEAR(columns::age_myr(c),42.5,1e-10);
        EXPECT_GT(c.ocean,0);
    }
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
}

TEST(TectonicOceanRelief, SubmarineTransportIsSeparateConservativeAndDoesNotChangeDryLand) {
    using namespace tectonics::columns;
    constexpr uint32_t w=32,h=16;
    std::vector<Cell> ocean(w*h,seed(-4000,false,40)),land(w*h);
    for (size_t i=0;i<land.size();++i) land[i]=seed(1000+100*std::sin(i%w*1.5),true,40);
    const std::vector<Displacement> stopped(w*h);
    Model a,b;
    a.initialize(w,h,150000,ocean,2,1,{5000,600});
    b.initialize(w,h,150000,ocean,2,1,{5000,600});
    a.advance(stopped,1,1,0,1); b.advance(stopped,1,1,0,.1);
    EXPECT_GT(a.budget().eroded_kg,0);
    EXPECT_NEAR(b.budget().eroded_kg,.1*a.budget().eroded_kg,a.budget().eroded_kg*1e-12);
    EXPECT_NEAR(b.budget().eroded_kg,b.budget().deposited_kg,b.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(b.budget().relative_mass_residual),1e-12);
    a.initialize(w,h,150000,land); b.initialize(w,h,150000,land);
    a.advance(stopped,1,1,1,1); b.advance(stopped,1,1,1,.1);
    for (size_t i=0;i<land.size();++i) {
        EXPECT_DOUBLE_EQ(a.cells()[i].continent,b.cells()[i].continent);
        EXPECT_DOUBLE_EQ(a.cells()[i].sediment,b.cells()[i].sediment);
    }
}

TEST(TectonicOceanRelief, IrregularSeamountProfilesArePositiveRepeatableMaterial) {
    using namespace tectonics;
    constexpr uint32_t w=64,h=32;
    ocean_relief::Parameters p{0,600,45000,350,3500,.65};
    const auto irregular=ocean_relief::crust(w,h,150000,1,p);
    EXPECT_EQ(irregular,ocean_relief::crust(w,h,150000,1,p));
    p.seamount_roughness=0;
    const auto smooth=ocean_relief::crust(w,h,150000,1,p);
    double difference=0;
    for (size_t i=0;i<irregular.size();++i) {
        EXPECT_GE(irregular[i],7000);
        if (smooth[i]==7000) EXPECT_DOUBLE_EQ(irregular[i],7000);
        difference+=std::abs(irregular[i]-smooth[i]);
    }
    EXPECT_GT(difference/irregular.size(),100);
    p.seamount_roughness=.65;
    columns::Model model;
    model.initialize(w,h,150000,std::vector<columns::Cell>(w*h,columns::seed(-4000,false,40)),4,1,p);
    const auto before=model.cells();
    model.advance(std::vector<columns::Displacement>(w*h),1,0,0);
    for (size_t i=0;i<before.size();++i) EXPECT_DOUBLE_EQ(model.cells()[i].ocean,before[i].ocean);
    p.seamount_roughness=1.01;
    EXPECT_THROW(ocean_relief::validate(p),std::runtime_error);
}
