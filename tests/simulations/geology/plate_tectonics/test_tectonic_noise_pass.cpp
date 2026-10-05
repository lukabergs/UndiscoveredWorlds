#include "gtest/gtest.h"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_noise_pass.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

TEST(TectonicNoisePass, ConservativeRepeatableContinuousAndNotAnAccumulatingRandomWalk) {
    using namespace tectonics;
    constexpr uint32_t w=64,h=32;
    const noise_pass::Parameters p{1000,600,180};
    std::vector<columns::Cell> initial(w*h,columns::seed(-4000,false,40));
    auto once=initial,steps=initial,replay=initial,other=initial;
    const double moved=noise_pass::apply(w,h,156250,once,1,0,600,p);
    EXPECT_GT(moved,0);
    for (int t=0;t<600;++t) noise_pass::apply(w,h,156250,steps,1,t,t+1,p);
    noise_pass::apply(w,h,156250,replay,1,0,600,p);
    noise_pass::apply(w,h,156250,other,2,0,600,p);
    double mass=0,original_mass=0,age_moment=0,original_moment=0,signal=0,difference=0;
    for (size_t i=0;i<once.size();++i) {
        EXPECT_DOUBLE_EQ(once[i].ocean,replay[i].ocean);
        EXPECT_NEAR(once[i].ocean,steps[i].ocean,1e-9);
        EXPECT_GE(once[i].ocean,0);
        EXPECT_DOUBLE_EQ(once[i].oceanic_area,1);
        mass+=columns::mass_per_area(once[i]);original_mass+=columns::mass_per_area(initial[i]);
        age_moment+=once[i].oceanic_age_volume;original_moment+=initial[i].oceanic_age_volume;
        signal+=std::pow(columns::elevation_m(once[i])-columns::elevation_m(initial[i]),2);
        difference+=std::abs(once[i].ocean-other[i].ocean);
    }
    EXPECT_NEAR(mass,original_mass,original_mass*1e-13);
    EXPECT_NEAR(age_moment,original_moment,original_moment*1e-13);
    EXPECT_GT(signal,0);EXPECT_GT(difference,0);
    auto continuous=once;
    noise_pass::apply(w,h,156250,continuous,1,600,600.00001,p);
    for (size_t i=0;i<once.size();++i) EXPECT_NEAR(once[i].ocean,continuous[i].ocean,.01);
    std::cout << "Noise potential 1000 m, stationary ocean height RMS " << std::sqrt(signal/once.size()) << " m\n";
}

TEST(TectonicNoisePass, DisabledIsExactAndModelIncludesTransfersInItsMassBalance) {
    using namespace tectonics;
    constexpr uint32_t w=16,h=8;
    std::vector<columns::Cell> cells(w*h,columns::seed(-4000,false,40));
    const auto initial=cells;
    EXPECT_DOUBLE_EQ(noise_pass::apply(w,h,150000,cells,1,0,1,{}),0);
    for (size_t i=0;i<cells.size();++i) EXPECT_DOUBLE_EQ(cells[i].ocean,initial[i].ocean);
    columns::Model a,b;
    a.initialize(w,h,150000,cells,4,1);b.initialize(w,h,150000,cells,4,1);
    const std::vector<columns::Displacement> velocity(w*h,{.1,-.1});
    for (int step=0;step<12;++step) {
        a.advance(velocity,1,0,0,0,nullptr,{}, {4000,600,180});
        b.advance(velocity,1,0,0,0,nullptr,{}, {4000,600,180});
    }
    EXPECT_GT(a.budget().noise_redistributed_kg,0);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual),1e-12);
    for (size_t i=0;i<cells.size();++i) EXPECT_DOUBLE_EQ(a.cells()[i].ocean,b.cells()[i].ocean);
    EXPECT_THROW(noise_pass::validate({-1,600,180}),std::runtime_error);
    EXPECT_THROW(noise_pass::validate({1,0,180}),std::runtime_error);
    EXPECT_THROW(noise_pass::validate({1,600,std::numeric_limits<double>::infinity()}),std::runtime_error);
    EXPECT_THROW(noise_pass::validate({1,600,180,-1}),std::runtime_error);
}

TEST(TectonicNoisePass, CoastalTaperPreservesNearSeaLevelReliefAndMass) {
    using namespace tectonics;
    constexpr uint32_t w=64,h=32;
    std::vector<columns::Cell> cells;
    for (uint32_t y=0;y<h;++y) for (uint32_t x=0;x<w;++x)
        cells.push_back(columns::seed((x%2 ? 1 : -1)*.01,true,40));
    const auto initial=cells;
    const double moved=noise_pass::apply(w,h,156250,cells,1,0,1,{4000,600,180,500});
    double before=0,after=0;
    EXPECT_GT(moved,0);
    for (size_t i=0;i<cells.size();++i) {
        const double z=columns::elevation_m(cells[i]);
        const double original=columns::elevation_m(initial[i]);
        EXPECT_EQ(z>=0,original>=0);
        EXPECT_NEAR(z,original,1e-6);
        before+=columns::mass_per_area(initial[i]);after+=columns::mass_per_area(cells[i]);
    }
    EXPECT_NEAR(before,after,before*1e-13);
}
