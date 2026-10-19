// Ocean relief (W2-E, E3): the periodic noise and the crust pattern (a torus raster, read on the
// cube through the 4N x 2N lat-lon proxy until E4), material relief that moves with the rock.
#include "gtest/gtest.h"
#include "cube_fixture.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_ocean_relief.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
using cubefixture::CellId;
using cubefixture::Vec3;
}

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

// The model's newborn crust is the proxy pattern at each cell: it varies over the sphere, stays
// put while the rock moves, and a rigid rotation carries the relief with the rock.
TEST(TectonicOceanRelief, ReliefIsMaterialAndDoesNotRegenerateOrDependOnUpdateCount) {
    using namespace tectonics;
    const auto grid = cubefixture::cube(16, cubefixture::radiusfor_km(16, 150000));
    std::vector<columns::Cell> cells(grid->cellcount(),columns::seed(-4000,false,40));
    columns::Model a,b;
    a.initialize(grid,cells,2,1,{5000,600});
    b.initialize(grid,cells,2,1,{5000,600});
    const auto initial=a.cells();
    double lowest=1e300,highest=0;
    for (const auto& c:initial) { lowest=std::min(lowest,c.ocean); highest=std::max(highest,c.ocean); }
    EXPECT_GT(highest-lowest,1000);
    const auto stopped=cubefixture::still(*grid);
    a.advance(stopped,1,0,0);
    for (int i=0;i<4;++i) b.advance(stopped,.25,0,0);
    for (size_t i=0;i<cells.size();++i) {
        EXPECT_DOUBLE_EQ(a.cells()[i].ocean,initial[i].ocean);
        // One update or four: the same up to the transport's total-restoring scale per update.
        EXPECT_NEAR(a.cells()[i].ocean,b.cells()[i].ocean,1e-12*a.cells()[i].ocean);
        EXPECT_NEAR(columns::age_myr(a.cells()[i]),41,1e-12);
    }
    // Twenty quarter columns east (five columns at the equator) about the z axis.
    const double column=cubefixture::pi/32;
    const Vec3 rotation{0,0,.25*column};
    for (int i=0;i<20;++i) a.advance(cubefixture::rotating(*grid,rotation),1,0,0);
    double error=0,signal=0;
    const auto& mesh=*grid->mesh();
    for (CellId c=0;c<cells.size();++c) {
        const Vec3 p=grid->position(c);
        // The rock now at p came from p turned back by the total rotation.
        const double angle=-5*column;
        const Vec3 from{std::cos(angle)*p.x-std::sin(angle)*p.y,std::sin(angle)*p.x+std::cos(angle)*p.y,p.z};
        const double expected=initial[mesh.locate(from)].ocean;
        error+=std::pow(a.cells()[c].ocean-expected,2);
        signal+=std::pow(expected-7000,2);
    }
    RecordProperty("relative_error",std::to_string(error/signal));
    // Measured 0.074 (N 16): the remap's smoothing of the proxy pattern over five columns.
    EXPECT_LT(error/signal,.1);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual),1e-12);
}

TEST(TectonicOceanRelief, SparseBasaltSeamountsFormIslandsAndRemainInTheMassBudget) {
    using namespace tectonics;
    const ocean_relief::Parameters relief{0,600,45000,350,3500};
    const auto rocks=ocean_relief::crust(128,64,100000,1,relief);
    EXPECT_EQ(rocks,ocean_relief::crust(128,64,100000,1,relief));
    const auto grid = cubefixture::cube(32, cubefixture::radiusfor_km(32, 100000));
    std::vector<columns::Cell> initial(grid->cellcount(),columns::seed(-4000,false,40));
    columns::Model model; model.initialize(grid,initial,2,1,relief);
    size_t islands=0; double mass=0;
    for (CellId c=0;c<model.cells().size();++c) {
        const auto& cell=model.cells()[c];
        islands+=columns::elevation_m(cell)>0;
        EXPECT_NEAR(columns::age_myr(cell),40,1e-12);
        mass+=columns::mass_per_area(cell)*grid->area_m2(c);
    }
    EXPECT_GT(islands,10U); EXPECT_LT(islands,initial.size()/10);
    EXPECT_NEAR(mass,model.budget().initial_mass_kg,mass*1e-13);
    const auto motion=cubefixture::rotating(*grid,{.1*cubefixture::pi/64,-.1*cubefixture::pi/64,0});
    for (int t=0;t<20;++t) model.advance(motion,1,1,1);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
}

TEST(TectonicOceanRelief, FineTransportSubstepsAgeRockOnlyOnceAndCloseMass) {
    using namespace tectonics;
    const auto grid = cubefixture::cube(8, cubefixture::radiusfor_km(8, 150000));
    std::vector<columns::Cell> cells(grid->cellcount(),columns::seed(-4000,false,40));
    columns::Model model;
    model.initialize(grid,cells,4,1,{5000,600});
    // Half a column per update over 2.5 Myr: four subcells per cell edge, several substeps.
    const double column=cubefixture::pi/16;
    model.advance(cubefixture::rotating(*grid,{.35*column,-.28*column,.21*column}),2.5,0,0);
    EXPECT_GT(model.transportDiagnostics().substeps,1u);
    for (const auto& c : model.cells()) {
        EXPECT_NEAR(columns::age_myr(c),42.5,1e-10);
        EXPECT_GT(c.ocean,0);
    }
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
}

TEST(TectonicOceanRelief, SubmarineTransportIsSeparateConservativeAndDoesNotChangeDryLand) {
    using namespace tectonics::columns;
    const auto grid = cubefixture::cube(8, cubefixture::radiusfor_km(8, 150000));
    std::vector<Cell> ocean(grid->cellcount(),seed(-4000,false,40)),land(grid->cellcount());
    for (size_t i=0;i<land.size();++i) land[i]=seed(1000+100*std::sin(i%32*1.5),true,40);
    const auto stopped=cubefixture::still(*grid);
    Model a,b;
    a.initialize(grid,ocean,2,1,{5000,600});
    b.initialize(grid,ocean,2,1,{5000,600});
    a.advance(stopped,1,1,0,1); b.advance(stopped,1,1,0,.1);
    EXPECT_GT(a.budget().eroded_kg,0);
    EXPECT_NEAR(b.budget().eroded_kg,.1*a.budget().eroded_kg,a.budget().eroded_kg*1e-12);
    EXPECT_NEAR(b.budget().eroded_kg,b.budget().deposited_kg,b.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(b.budget().relative_mass_residual),1e-12);
    a.initialize(grid,land); b.initialize(grid,land);
    a.advance(stopped,1,1,1,1); b.advance(stopped,1,1,1,.1);
    for (size_t i=0;i<land.size();++i) {
        EXPECT_DOUBLE_EQ(a.cells()[i].continent,b.cells()[i].continent);
        EXPECT_DOUBLE_EQ(a.cells()[i].sediment,b.cells()[i].sediment);
    }
}

TEST(TectonicOceanRelief, IrregularSeamountProfilesArePositiveRepeatableMaterial) {
    using namespace tectonics;
    ocean_relief::Parameters p{0,600,45000,350,3500,.65};
    const auto irregular=ocean_relief::crust(64,32,150000,1,p);
    EXPECT_EQ(irregular,ocean_relief::crust(64,32,150000,1,p));
    p.seamount_roughness=0;
    const auto smooth=ocean_relief::crust(64,32,150000,1,p);
    double difference=0;
    for (size_t i=0;i<irregular.size();++i) {
        EXPECT_GE(irregular[i],7000);
        if (smooth[i]==7000) EXPECT_DOUBLE_EQ(irregular[i],7000);
        difference+=std::abs(irregular[i]-smooth[i]);
    }
    EXPECT_GT(difference/irregular.size(),100);
    p.seamount_roughness=.65;
    const auto grid = cubefixture::cube(16, cubefixture::radiusfor_km(16, 150000));
    columns::Model model;
    model.initialize(grid,std::vector<columns::Cell>(grid->cellcount(),columns::seed(-4000,false,40)),4,1,p);
    const auto before=model.cells();
    model.advance(cubefixture::still(*grid),1,0,0);
    for (size_t i=0;i<before.size();++i) EXPECT_DOUBLE_EQ(model.cells()[i].ocean,before[i].ocean);
    p.seamount_roughness=1.01;
    EXPECT_THROW(ocean_relief::validate(p),std::runtime_error);
}
