#include "gtest/gtest.h"
#include "tectonic_deformation.hpp"
#include "tectonic_columns.hpp"
#include <algorithm>
#include <cmath>

TEST(TectonicDeformation, DistanceCurveHasAnInlandStrainPeakAndCompactSupport) {
    using tectonics::deformation::taper;
    EXPECT_DOUBLE_EQ(taper(0,300000),1);
    EXPECT_DOUBLE_EQ(taper(600000,300000),0);
    EXPECT_DOUBLE_EQ(taper(900000,300000),0);
    const double shore=taper(0,300000)-taper(10000,300000);
    const double inland=taper(295000,300000)-taper(305000,300000);
    EXPECT_GT(inland,shore*10);
    EXPECT_NEAR(taper(150000,300000),taper(300000,600000),1e-15);
}

TEST(TectonicDeformation, TransformDivergenceAndDisabledProfileLeaveMaterialMotionUnchanged) {
    using namespace tectonics;
    constexpr uint32_t w=32,h=8;
    std::vector<uint32_t> owners(w*h); std::vector<double> continent(w*h,1);
    std::vector<transport::Displacement> motion(w*h);
    for (uint32_t i=0;i<w*h;++i) { owners[i]=i%w<16?0:1; motion[i]={0,owners[i]?.1:-.1}; }
    const auto shear=deformation::distribute(w,h,50000,owners,continent,motion,1,{});
    for (size_t i=0;i<motion.size();++i) { EXPECT_DOUBLE_EQ(shear[i].x,motion[i].x); EXPECT_DOUBLE_EQ(shear[i].y,motion[i].y); }
    for (uint32_t i=0;i<w*h;++i) motion[i]={owners[i]?-.1:.1,0};
    const auto disabled=deformation::distribute(w,h,50000,owners,continent,motion,1,{0,.2,.2});
    for (size_t i=0;i<motion.size();++i) EXPECT_DOUBLE_EQ(disabled[i].x,motion[i].x);
    // Seam is divergent, centre is convergent; only the latter spreads shortening.
    const auto active=deformation::distribute(w,h,50000,owners,continent,motion,1,{200,0,.2});
    EXPECT_DOUBLE_EQ(active[0].x,motion[0].x);
    EXPECT_LT(active[15].x,motion[15].x); EXPECT_GT(active[16].x,motion[16].x);
}

TEST(TectonicDeformation, SeededWidthIsRepeatableAndDoesNotChangeWithACommonTranslationVelocity) {
    using namespace tectonics;
    constexpr uint32_t w=64,h=32;
    std::vector<uint32_t> owners(w*h); std::vector<double> continent(w*h,1);
    std::vector<transport::Displacement> v(w*h),offset(w*h);
    for (uint32_t i=0;i<w*h;++i) { owners[i]=i%w<32?0:1; v[i]={owners[i]?-.1:.1,0}; offset[i]={v[i].x+.2,.3}; }
    const auto a=deformation::distribute(w,h,50000,owners,continent,v,1,{});
    const auto b=deformation::distribute(w,h,50000,owners,continent,v,1,{});
    const auto moved=deformation::distribute(w,h,50000,owners,continent,offset,1,{});
    double along_strike_range=0;
    for (size_t i=0;i<v.size();++i) {
        EXPECT_DOUBLE_EQ(a[i].x,b[i].x);
        EXPECT_NEAR(moved[i].x,a[i].x+.2,1e-14); EXPECT_NEAR(moved[i].y,.3,1e-14);
        along_strike_range=std::max(along_strike_range,std::abs(a[i].x-a[i%w].x));
    }
    EXPECT_GT(along_strike_range,1e-4);
}

TEST(TectonicDeformation, OceanConvergenceBuildsAnInlandBeltWithoutCreatingContinentalMass) {
    using namespace tectonics;
    constexpr uint32_t w=128,h=8;
    std::vector<uint32_t> owners(w*h); std::vector<double> continent(w*h);
    std::vector<transport::Displacement> rigid(w*h);
    std::vector<columns::Cell> initial(w*h);
    double mass=0;
    for (uint32_t i=0;i<w*h;++i) {
        const bool land=i%w>=32 && i%w<96;
        owners[i]=land?1:0; continent[i]=land?1:0; rigid[i]={land?-.02:.02,0};
        initial[i]=columns::seed(1000,land,100); mass+=initial[i].continent;
    }
    const auto v=deformation::distribute(w,h,50000,owners,continent,rigid,1,{300,0,.2});
    columns::Model a,b; a.initialize(w,h,50000,initial); b.initialize(w,h,50000,initial);
    for (int t=0;t<100;++t) { a.advance(v,.25,0,1); b.advance(rigid,.25,0,1); }
    uint32_t peak=32; double final_mass=0,control_peak=0;
    for (uint32_t i=0;i<w*h;++i) final_mass+=a.cells()[i].continent;
    for (uint32_t x=32;x<64;++x) {
        if (columns::elevation_m(a.cells()[x])>columns::elevation_m(a.cells()[peak])) peak=x;
        control_peak=std::max(control_peak,columns::elevation_m(b.cells()[x]));
    }
    EXPECT_GT(peak,34U); EXPECT_LT(peak,46U);
    EXPECT_GT(columns::elevation_m(a.cells()[peak]),1100);
    EXPECT_LT(columns::elevation_m(a.cells()[32]),control_peak);
    EXPECT_NEAR(final_mass,mass,mass*1e-12);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual),1e-12);
    RecordProperty("peak_cell_from_initial_contact",peak-32);
}

TEST(TectonicColumns, PersistentSubcellsRetainFourCellTerrainWithoutNoiseOrExtraMass) {
    using namespace tectonics::columns;
    constexpr uint32_t w=64,h=8; constexpr double pi=3.141592653589793;
    std::vector<Cell> initial(w*h);
    for (uint32_t i=0;i<w*h;++i) initial[i]=seed(1000+100*std::sin(2*pi*(i%w)/4),true,100);
    double retained[2]{};
    for (int mode=0;mode<2;++mode) {
        Model model; model.initialize(w,h,1e7,initial,mode?2:1); // Isolate transport from root relaxation.
        for (int t=0;t<200;++t) model.advance(std::vector<Displacement>(w*h,{.16,-.12}),1,0,0);
        double signal=0;
        for (uint32_t i=0;i<w*h;++i) {
            signal+=(elevation_m(model.cells()[i])-1000)*std::sin(2*pi*(i%w)/4);
            EXPECT_NEAR(age_myr(model.cells()[i]),300,1e-8);
            EXPECT_GE(elevation_m(model.cells()[i]),899.9);
            EXPECT_LE(elevation_m(model.cells()[i]),1100.1);
        }
        retained[mode]=signal/(w*h*50);
        EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-11);
    }
    RecordProperty("native_signal_retained",retained[0]); RecordProperty("subcell_signal_retained",retained[1]);
    EXPECT_LT(std::abs(retained[0]),.01); EXPECT_GT(retained[1],.40);
}
