// Active margins on the cube (W2-E, E3): trenches, arcs and ridge heat from the margin forcing,
// rivers and waves, rifts, the continental thickness caps and the orogenic collapse, the initial
// terrain, on the equatorial band of the cube (cube_fixture.hpp).
#include "gtest/gtest.h"
#include "cube_fixture.hpp"
#include "simulations/geology/plate_tectonics/tectonic_columns.hpp"
#include "simulations/geology/plate_tectonics/tectonic_margins.hpp"
#include "simulations/geology/plate_tectonics/tectonic_initialization.hpp"
#include "simulations/geology/plate_tectonics/tectonic_variation.hpp"
#include "core/mesh/cube_geometry.hpp"
#include "core/mesh/mesh_geometry.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using cubefixture::Band;
using cubefixture::CellId;
using cubefixture::Vec3;
namespace st = tectonics::sphere_transport;
using tectonics::columns::Model;

// Two oceanic plates on a band of N 16 with 50 km columns at the equator (64 columns), split by
// the plane y = 0: plate 0 (columns 56..63 and 0..23, old floor of 100 Myr) moves east and
// plate 1 (columns 24..55, 20 Myr) west by 0.2 column per update, so they close at longitude 0
// (between columns 23 and 24) and part at 180 degrees (between 55 and 56). Rows 7 and 8
// straddle the equator.
struct MarginFixture {
    Band band{16, cubefixture::radiusfor_km(16, 50000)};
    tectonics::columns::Metric metric{*band.grid};
    std::vector<uint32_t> owners = cubefixture::hemispheres(*band.grid);
    std::vector<tectonics::columns::Cell> cells;
    st::Motion motion;
    std::vector<Vec3> velocity;
    static constexpr int west_of_trench = 23, east_of_trench = 24, east_of_ridge = 56, row = 8;
    MarginFixture() {
        cells.resize(band.cellcount());
        for (CellId c = 0; c < cells.size(); ++c) cells[c] = tectonics::columns::seed(-4000, false, owners[c] == 0 ? 100 : 20);
        move(.2, -.2);
    }
    // Plate 0 east by u0, plate 1 east by u1 (columns per update).
    void move(double u0, double u1) {
        motion = cubefixture::plates(owners, {band.eastward(u0), band.eastward(u1)});
        velocity = Model::displacement(*band.grid, motion);
    }
    CellId at(int x, int y = row) const { return band.cell(x, y); }
    double area() const { return band.spacing_m() * band.spacing_m(); }
    tectonics::margins::Forcing forcing(const std::vector<tectonics::columns::Cell>& state,
                                        const tectonics::margins::Parameters& p,
                                        tectonics::margins::Polarity* polarity = nullptr, long seed = 0,
                                        double time = 0) const {
        return tectonics::margins::forcing(metric, 1, owners, state, velocity, p, polarity, seed, time);
    }
    double total(const Model& model, double tectonics::columns::Cell::* member) const {
        double sum = 0;
        for (CellId c = 0; c < model.cells().size(); ++c) sum += model.cells()[c].*member * band.grid->area_m2(c);
        return sum;
    }
};
}

TEST(TectonicMargins, SubductionPolarityPlacesTrenchAndArcOnOppositePlates) {
    MarginFixture f;
    tectonics::margins::Parameters p; p.enabled=true;
    auto a=f.forcing(f.cells,p);
    const CellId west=f.at(f.west_of_trench), east=f.at(f.east_of_trench);
    EXPECT_GT(a.trench_m[west],0);
    EXPECT_DOUBLE_EQ(a.trench_m[east],0);
    ASSERT_NE(a.arc_receiver[west],std::numeric_limits<uint32_t>::max());
    EXPECT_EQ(f.owners[a.arc_receiver[west]],1U);
    EXPECT_GT(f.band.column(a.arc_receiver[west]),f.east_of_trench);
    EXPECT_LT(a.thermal_age_myr[f.at(f.east_of_ridge)],100);
    EXPECT_DOUBLE_EQ(a.thermal_age_myr[west],-1);
    const auto b=f.forcing(f.cells,p);
    EXPECT_EQ(a.trench_m,b.trench_m);EXPECT_EQ(a.arc_receiver,b.arc_receiver);
    tectonics::margins::Polarity polarity;
    f.forcing(f.cells,p,&polarity);
    for (size_t i=0;i<f.cells.size();++i) f.cells[i].oceanic_age_volume=f.cells[i].ocean*(f.owners[i]==0?10:200);
    const auto persistent=f.forcing(f.cells,p,&polarity);
    EXPECT_EQ(f.owners[persistent.arc_receiver[west]],1U);
    EXPECT_EQ(persistent.arc_receiver[east],persistent.arc_receiver[west]);
    EXPECT_DOUBLE_EQ(persistent.trench_m[east],0);
    // Pure slip along the boundary (the two plates turning about the boundary plane's normal,
    // y): no trench, no arc and no ridge heat.
    f.motion=cubefixture::plates(f.owners,{{0,.2*cubefixture::pi/32,0},{0,-.2*cubefixture::pi/32,0}});
    f.velocity=Model::displacement(*f.band.grid,f.motion);
    a=f.forcing(f.cells,p);
    for (size_t i=0;i<f.cells.size();++i) {
        EXPECT_NEAR(a.trench_m[i],0,1e-6);
        EXPECT_TRUE(a.thermal_age_myr[i]<0 || a.thermal_age_myr[i]>1e6) << i;
    }
}

TEST(TectonicMargins, ThermalHeatingPreservesGeologicAgeAndFlexureHasContinuousResponse) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;p.trench_response_myr=1.5;
    auto forcing=f.forcing(f.cells,p);
    Model a,b;
    a.initialize(f.band.grid,f.cells);b.initialize(f.band.grid,f.cells);
    const auto stopped=cubefixture::still(*f.band.grid);
    a.advance(stopped,1,0,1,0,&forcing,p);
    for (int step=0;step<4;++step) b.advance(stopped,.25,0,1,0,&forcing,p);
    const CellId west=f.at(f.west_of_trench), ridge=f.at(f.east_of_ridge);
    EXPECT_NEAR(a.cells()[west].flexure_m,b.cells()[west].flexure_m,1e-10);
    EXPECT_GT(a.cells()[west].flexure_m,0);
    EXPECT_LT(a.cells()[west].flexure_m,forcing.trench_m[west]);
    EXPECT_GT(a.cells()[ridge].thermal_rejuvenation_volume,0);
    EXPECT_NEAR(tectonics::columns::age_myr(a.cells()[ridge]),101,1e-10);
    EXPECT_DOUBLE_EQ(a.budget().generated_arc_kg,0);
    EXPECT_LT(std::abs(a.budget().relative_mass_residual),1e-12);
}

TEST(TectonicMargins, ArcRockAccumulatesOnlyAfterSubductionAndCanEmerge) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    // Reference frame of a stationary overriding plate and fixed trench.
    // The older slab approaches at 3 cm/year (0.6 column per update: the cube allows 0.95 of a
    // cell's smallest edge, about 0.71 column at a face edge); no island exists initially.
    f.move(.6,0);
    Model model;model.initialize(f.band.grid,f.cells);
    tectonics::margins::Polarity polarity;
    for (const auto& c:model.cells()) EXPECT_LT(tectonics::columns::elevation_m(c),0);
    for (int step=0;step<267;++step) {
        const auto forcing=f.forcing(model.cells(),p,&polarity);
        model.advance(f.motion,1,0,1,0,&forcing,p);
    }
    EXPECT_GT(model.budget().generated_arc_kg,0);
    EXPECT_LE(model.budget().generated_arc_kg,model.budget().generated_ocean_kg);
    EXPECT_LE(model.budget().generated_arc_kg,model.budget().subducted_ocean_kg*p.arc_melt_ratio*(1+1e-12));
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-11);
    double peak=-1e10;
    for (const auto& c:model.cells()) peak=std::max(peak,tectonics::columns::elevation_m(c));
    EXPECT_GT(peak,0);
}

namespace {
// The river fixture: a continental ramp 2000 - 90 x m on columns 0..19 of a band with 150 km
// columns, ocean at -4000 m elsewhere (the polar faces too).
struct RiverFixture {
    Band band{16, cubefixture::radiusfor_km(16, 150000)};
    std::vector<tectonics::columns::Cell> cells = cubefixture::layout(band, tectonics::columns::seed(-4000,false,40),
        [](int x, int) { return x<20?tectonics::columns::seed(2000-90*x,true,40):tectonics::columns::seed(-4000,false,40); });
    st::Motion stopped = cubefixture::still(*band.grid);
    CellId at(int x, int y = 8) const { return band.cell(x, y); }
};
}

TEST(TectonicMargins, RiversDeliverSedimentOffshoreWithoutSmoothingQuietSeafloor) {
    RiverFixture f;
    tectonics::margins::Parameters p;p.enabled=true;
    Model model;model.initialize(f.band.grid,f.cells);
    model.advance(f.stopped,1,0,1,0,nullptr,p);
    EXPECT_GT(model.budget().river_eroded_kg,0);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_GT(model.cells()[f.at(20)].sediment,0);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    f.cells.assign(f.band.cellcount(),tectonics::columns::seed(-4000,false,40));
    model.initialize(f.band.grid,f.cells,2,1,{4000,600});
    const auto initial=model.cells();
    model.advance(f.stopped,1,0,0,0,nullptr,p);
    EXPECT_DOUBLE_EQ(model.budget().river_eroded_kg,0);
    for (size_t i=0;i<f.cells.size();++i) EXPECT_DOUBLE_EQ(initial[i].ocean,model.cells()[i].ocean);
}

TEST(TectonicMargins, ImplicitRiversDeliverSedimentOffshoreAndCloseTheBudget) {
    RiverFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_implicit=true;
    Model model;model.initialize(f.band.grid,f.cells);
    model.advance(f.stopped,1,0,1,0,nullptr,p);
    EXPECT_GT(model.budget().river_eroded_kg,0);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_GT(model.cells()[f.at(20)].sediment,0);
    EXPECT_EQ(model.budget().river_unsupplied_cells,0U);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    const double uniform_kg=model.budget().river_eroded_kg;
    RecordProperty("uniform_river_kg",std::to_string(uniform_kg));
    // Per-cell erodibility from the process variation acts and keeps the budget closed.
    p.process_variation=.35;
    const tectonics::columns::Metric metric(*f.band.grid);
    const auto owners=cubefixture::hemispheres(*f.band.grid);
    const std::vector<Vec3> rigid(f.cells.size());
    auto forcing=tectonics::margins::forcing(metric,1,owners,f.cells,rigid,p,nullptr,1,20);
    model.initialize(f.band.grid,f.cells);
    model.advance(f.stopped,1,0,1,0,&forcing,p);
    const double varied_kg=model.budget().river_eroded_kg;
    RecordProperty("varied_river_kg",std::to_string(varied_kg));
    EXPECT_GT(std::abs(varied_kg-uniform_kg),1e-5*uniform_kg);
    EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    p.process_variation=0;
    f.cells.assign(f.band.cellcount(),tectonics::columns::seed(-4000,false,40));
    model.initialize(f.band.grid,f.cells,2,1,{4000,600});
    const auto initial=model.cells();
    model.advance(f.stopped,1,0,0,0,nullptr,p);
    EXPECT_DOUBLE_EQ(model.budget().river_eroded_kg,0);
    EXPECT_DOUBLE_EQ(model.budget().deposited_kg,0);
    for (size_t i=0;i<f.cells.size();++i) EXPECT_DOUBLE_EQ(initial[i].ocean,model.cells()[i].ocean);
}

TEST(TectonicMargins, ImplicitAndExplicitRiversAgreeOnASmallStep) {
    // Same law, same receivers: with outlets at sea level the explicit step's
    // erosion is the implicit one's first-order limit as dt goes to zero.
    RiverFixture f;
    const auto eroded=[&](bool implicit) {
        tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=15;
        p.river_outlet_sea_level=true;p.river_implicit=implicit;
        Model model;model.initialize(f.band.grid,f.cells);
        model.advance(f.stopped,.01,0,1,0,nullptr,p);
        return model.budget().river_eroded_kg;
    };
    const double explicit_kg=eroded(false),implicit_kg=eroded(true);
    EXPECT_GT(explicit_kg,0);
    RecordProperty("implicit_over_explicit",std::to_string(implicit_kg/explicit_kg));
    // The implicit step divides each drop by 1 + K sqrt(A) dt / L; the routes agree on the
    // band's meridian columns. Measured on the cube (N 16, 150 km): 0.999923.
    EXPECT_LT(implicit_kg,explicit_kg);
    EXPECT_NEAR(implicit_kg/explicit_kg,1,1e-4);
}

TEST(TectonicMargins, ImplicitRiversKeepTheExplicitWavePass) {
    RiverFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_implicit=true;
    Model rivers;rivers.initialize(f.band.grid,f.cells);
    rivers.advance(f.stopped,1,0,1,0,nullptr,p);
    p.wave_erosion_m_myr=50;
    Model waves;waves.initialize(f.band.grid,f.cells);
    waves.advance(f.stopped,1,0,1,0,nullptr,p);
    // river_eroded_kg also counts the explicit pass's wave attack.
    const double wave_kg=waves.budget().river_eroded_kg-rivers.budget().river_eroded_kg;
    RecordProperty("wave_kg",std::to_string(wave_kg));
    EXPECT_GT(wave_kg,1e-3*rivers.budget().river_eroded_kg);
    EXPECT_NEAR(waves.budget().eroded_kg,waves.budget().deposited_kg,waves.budget().eroded_kg*1e-12);
    EXPECT_LT(std::abs(waves.budget().relative_mass_residual),1e-12);
    p.alluvial_transport=1;
    EXPECT_THROW(waves.advance(f.stopped,1,0,1,0,nullptr,p),std::runtime_error);
}

TEST(TectonicMargins, ImplicitRiversPassAClosedDepressionsLoadToItsSpill) {
    // A pit at 500 m in row 8 of the ramp (neighbours 1010 to 1190 m) drains through its rim to
    // the sea; the explicit path's pit keeps what reaches it.
    RiverFixture f;
    const CellId pit=f.at(10),shelf8=f.at(20),shelf3=f.at(20,3);
    f.cells[pit]=tectonics::columns::seed(500,true,40);
    const auto run=[&](bool implicit) {
        tectonics::margins::Parameters p;p.enabled=true;p.river_implicit=implicit;
        Model model;model.initialize(f.band.grid,f.cells);
        model.advance(f.stopped,1,0,1,0,nullptr,p);
        EXPECT_NEAR(model.budget().eroded_kg,model.budget().deposited_kg,model.budget().eroded_kg*1e-12);
        EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
        return model.cells();
    };
    const auto implicit=run(true),explicit_=run(false);
    // The pit's catchment (row 8 west of it and its neighbours) reaches the shelf of row 8.
    EXPECT_DOUBLE_EQ(implicit[pit].sediment,f.cells[pit].sediment);
    RecordProperty("shelf_row8_m",std::to_string(implicit[shelf8].sediment));
    RecordProperty("shelf_row3_m",std::to_string(implicit[shelf3].sediment));
    EXPECT_GT(implicit[shelf8].sediment,implicit[shelf3].sediment);
    EXPECT_GT(explicit_[pit].sediment,f.cells[pit].sediment);
}

TEST(TectonicMargins, ContinentalRiftsCreateOceanFloorWhenCrustBecomesThin) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    f.cells.assign(f.band.cellcount(),tectonics::columns::seed(100,true,40));
    f.move(-.15,.15);  // apart at longitude 0, together at 180 degrees
    Model model;model.initialize(f.band.grid,f.cells);
    const auto area_terms=[](const tectonics::columns::Budget& b) {
        return b.stretched_continental_area_m2+b.arc_continental_area_m2+b.continentalized_oceanic_area_m2
            -b.collided_continental_area_m2-b.oceanized_continental_area_m2;
    };
    for (int step=0;step<100;++step) {
        // The continental area changes by the budget's area terms only (the transport moves it).
        const double area=f.total(model,&tectonics::columns::Cell::continental_area),terms=area_terms(model.budget());
        model.advance(f.motion,1,0,1,0,nullptr,p);
        EXPECT_NEAR(f.total(model,&tectonics::columns::Cell::continental_area)-area,area_terms(model.budget())-terms,
                    1e-9*f.area());
    }
    // The rift stretches the continent, the far contact collides; no arcs, no caps.
    EXPECT_GT(model.budget().stretched_continental_area_m2,0);
    EXPECT_GT(model.budget().collided_continental_area_m2,0);
    EXPECT_DOUBLE_EQ(model.budget().arc_continental_area_m2,0);
    EXPECT_DOUBLE_EQ(model.budget().oceanized_continental_area_m2,0);
    EXPECT_DOUBLE_EQ(model.budget().continentalized_oceanic_area_m2,0);
    EXPECT_GT(model.budget().generated_ocean_kg,0);
    const CellId rift=f.at(f.west_of_trench);
    EXPECT_GT(model.cells()[rift].ocean,0);
    EXPECT_LT(model.cells()[rift].continental_area,.5);
    EXPECT_GT(tectonics::columns::elevation_m(model.cells()[rift]),-6000);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-11);
}

TEST(TectonicInitialization, PreservesDetailedTerrainWithoutGloballyYoungOcean) {
    MarginFixture f;
    std::vector<double> z(f.cells.size()),age(z.size(),30);
    for (size_t i=0;i<z.size();++i) z[i]=-4200+600*std::sin(i*.7);
    const std::vector<Vec3> stopped(z.size());
    const auto cells=tectonics::initialization::terrain(f.metric,1,z,age,f.owners,stopped,1);
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
    const auto moving=tectonics::initialization::terrain(f.metric,1,z,age,f.owners,f.velocity,1);
    // The ridge at 180 degrees opens, the contact at 0 closes. Far-field relief is unchanged.
    const CellId ridge=f.at(f.east_of_ridge),trench=f.at(f.west_of_trench);
    EXPECT_GT(tectonics::columns::elevation_m(moving[ridge]),z[ridge]);
    EXPECT_NEAR(tectonics::columns::elevation_m(moving[trench]),z[trench],1e-8);
    EXPECT_EQ(cells[trench].ocean,moving[trench].ocean);
    z.assign(z.size(),-200);
    const auto shelf=tectonics::initialization::terrain(f.metric,1,z,age,f.owners,stopped,1);
    EXPECT_EQ(shelf[0].continental_area,1);
    EXPECT_NEAR(tectonics::columns::elevation_m(shelf[0]),-200,1e-8);
    z.assign(z.size(),1000);
    for (int y=0;y<f.band.rows();++y) for (int x : {f.east_of_ridge, f.east_of_ridge+1}) z[f.band.cell(x,y)]=-6000;
    const auto edge=tectonics::initialization::terrain(f.metric,1,z,age,f.owners,f.velocity,1);
    for (size_t i=0;i<z.size();++i) {
        EXPECT_GE(tectonics::columns::elevation_m(edge[i]),-6000-1e-8);
        EXPECT_GT(edge[i].continent+edge[i].ocean,0);
    }
}

TEST(TectonicInitialization, GradualMarginsPreserveHeightAndAvoidCoolingDiscontinuity) {
    using namespace tectonics;
    MarginFixture f;
    std::vector<double> z(f.cells.size(),-1500),age(z.size(),40);
    for (int y=0;y<f.band.rows();++y) for (int x=0;x<f.band.columns();++x)
        z[f.band.cell(x,y)]=-1500-2000.0*x/(f.band.columns()-1);
    const std::vector<Vec3> stopped(z.size());
    const auto cells=initialization::terrain(f.metric,1,z,age,f.owners,stopped,1,2000);
    for (size_t i=0;i<z.size();++i) {
        EXPECT_NEAR(columns::elevation_m(cells[i]),z[i],1e-8);
        EXPECT_NEAR(cells[i].continental_area+cells[i].oceanic_area,1,1e-14);
        EXPECT_GE(cells[i].continent,0); EXPECT_GE(cells[i].ocean,0);
    }
    for (int x=1;x<f.band.columns();++x)
        EXPECT_GE(cells[f.at(x)].oceanic_area,cells[f.at(x-1)].oceanic_area);
    // A lithospheric root keeps the taper's heights: both end members carry it.
    const auto rooted=initialization::terrain(f.metric,1,z,age,f.owners,stopped,1,2000,1,246);
    for (size_t i=0;i<z.size();++i) EXPECT_NEAR(columns::elevation_m(rooted[i],1,246),z[i],1e-8);
    EXPECT_THROW(initialization::terrain(f.metric,1,z,age,f.owners,stopped,1,2000,1,500),std::runtime_error);
    const auto cooling=[&](double height) {
        z.assign(z.size(),height);
        auto c=initialization::terrain(f.metric,1,z,age,f.owners,stopped,1,2000)[0];
        c.oceanic_age_volume+=c.ocean*10;
        return height-columns::elevation_m(c);
    };
    EXPECT_NEAR(cooling(-2500.001),cooling(-2499.999),.01);
    EXPECT_DOUBLE_EQ(cells[f.at(0)].oceanic_area,0);
    EXPECT_DOUBLE_EQ(cells[f.at(f.band.columns()-1)].continental_area,0);
    EXPECT_THROW(initialization::terrain(f.metric,1,z,age,f.owners,stopped,1,4001),std::runtime_error);
    columns::Model model;model.initialize(f.band.grid,cells,4,1,{},true);
    model.advance(cubefixture::still(*f.band.grid),1,0,1);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
}

TEST(TectonicMargins, StrongerBendingIsLocalizedAndDoesNotConsumeStationaryRock) {
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.river_erosion=0;
    p.trench_load_n_m=3e13;p.trench_response_myr=1.5;
    auto load=f.forcing(f.cells,p);
    const CellId trench=f.at(f.west_of_trench),interior=f.at(12),east=f.at(f.east_of_trench);
    EXPECT_GT(load.trench_m[trench],2000);
    EXPECT_DOUBLE_EQ(load.trench_m[interior],0);
    EXPECT_DOUBLE_EQ(load.trench_m[east],0);
    Model model;model.initialize(f.band.grid,f.cells);
    const auto stopped=cubefixture::still(*f.band.grid);
    for (int step=0;step<5;++step) model.advance(stopped,1,0,0,0,&load,p);
    EXPECT_LT(tectonics::columns::elevation_m(model.cells()[trench],0),-4500);
    EXPECT_NEAR(tectonics::columns::elevation_m(model.cells()[interior],0),-2500,1e-9);
    EXPECT_NEAR(tectonics::columns::elevation_m(model.cells()[east],0),-2500,1e-9);
    EXPECT_DOUBLE_EQ(model.budget().subducted_ocean_kg,0);
    EXPECT_NEAR(model.budget().relative_mass_residual,0,1e-12);
}

// On the sphere (E4): deterministic, periodic in absolute time, bounded, and continuous across
// the cube's face edges and corners: a pure function of the unit vector whose 3-D gradient is
// at most amplitude x 2 pi / wavelength x (0.5 sqrt 2 + 0.5 sqrt 5) per km, so the two cells of
// every edge (those across face edges and around the 24 corner cells included) differ by at
// most that bound times their distance.
TEST(TectonicVariation, PeriodicDeterministicAndContinuousInAbsoluteTime) {
    const int n=16;
    const auto grid=cubefixture::cube(n);
    const double radius_km=grid->radius_m()/1000;
    const auto value=[&](const Vec3& p,double t,uint32_t channel=0) {
        return tectonics::variation::multiplier(p,radius_km,t,1,channel,.35,1800,80);
    };
    double low=2,high=0;
    for (CellId c=0;c<grid->cellcount();++c) {
        const Vec3 p=grid->position(c);
        const double t=(c%100)*.9,a=value(p,t);
        EXPECT_DOUBLE_EQ(a,value(p,t));
        EXPECT_NEAR(a,value(p,t+80),1e-12);
        EXPECT_NEAR(a,value(p,t+.001),1e-4);
        EXPECT_GE(a,.65);EXPECT_LE(a,1.35);
        low=std::min(low,a);high=std::max(high,a);
    }
    EXPECT_GT(high-low,.2);
    const double bound=.35*2*cubefixture::pi/1800*(.5*std::sqrt(2.)+.5*std::sqrt(5.));
    size_t faceedges=0;
    double largest=0;
    for (const auto& e : grid->edges()) {
        const double step=std::abs(value(grid->position(e.a),40)-value(grid->position(e.b),40));
        const double limit=bound*e.distance_m/1000;
        EXPECT_LE(step,limit*(1+1e-9));
        largest=std::max(largest,step/limit);
        faceedges+=e.a/(n*n)!=e.b/(n*n);
    }
    EXPECT_EQ(faceedges,12u*n);
    RecordProperty("largest_step_over_bound",std::to_string(largest));
    // The same pattern from either side of a cube corner (the unit vector alone decides).
    const Vec3 corner=mesh::normalise(Vec3{1,1,1});
    EXPECT_NEAR(value(corner,40),value(mesh::normalise(Vec3{1+1e-12,1,1}),40),1e-9);
    const Vec3 p=grid->position(5);
    EXPECT_NE(value(p,40,0),value(p,40,1));
    EXPECT_DOUBLE_EQ(tectonics::variation::multiplier(p,radius_km,40,1,0,0,1800,80),1);
    EXPECT_THROW(tectonics::variation::validate(.8,1800,80),std::runtime_error);
}

TEST(TectonicVariation, ChangesProcessRatesWithoutBreakingSedimentMassBudget) {
    RiverFixture rivers;
    MarginFixture f;
    tectonics::margins::Parameters p;p.enabled=true;p.process_variation=.35;
    const tectonics::columns::Metric metric(*rivers.band.grid);
    const auto owners=cubefixture::hemispheres(*rivers.band.grid);
    const auto velocity=Model::displacement(*rivers.band.grid,cubefixture::plates(owners,{rivers.band.eastward(.2),rivers.band.eastward(-.2)}));
    auto forcing=tectonics::margins::forcing(metric,1,owners,rivers.cells,velocity,p,nullptr,1,20);
    const auto replay=tectonics::margins::forcing(metric,1,owners,rivers.cells,velocity,p,nullptr,1,20);
    EXPECT_EQ(forcing.trench_m,replay.trench_m);
    EXPECT_EQ(forcing.erosion_multiplier,replay.erosion_multiplier);
    EXPECT_NE(forcing.erosion_multiplier,forcing.sediment_multiplier);
    Model model;model.initialize(rivers.band.grid,rivers.cells);
    model.advance(rivers.stopped,1,0,1,0,&forcing,p);
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
// An oceanic plate (0) subducting at 3 cm/yr under a stationary continental plate (1) with
// 500 m of sediment.
struct ContinentalMarginFixture : MarginFixture {
    ContinentalMarginFixture() {
        for (CellId i=0;i<cells.size();++i)
            if (owners[i]==1) { cells[i]=tectonics::columns::seed(500,true,100);cells[i].sediment=500; }
        move(.6,0);
    }
};
}

TEST(TectonicMargins, ContinentalMinimumThicknessCapsGhostAreaWithoutMovingRock) {
    MarginFixture f;
    for (auto& c:f.cells) {
        c=tectonics::columns::Cell{};
        c.continental_area=1;c.continent=1000;c.continental_age_volume=1000*50;
    }
    const auto stopped=cubefixture::still(*f.band.grid);
    Model legacy,capped;
    legacy.initialize(f.band.grid,f.cells);capped.initialize(f.band.grid,f.cells);
    tectonics::margins::Parameters off;off.enabled=true;off.river_erosion=0;
    const auto p=growth_parameters();
    legacy.advance(stopped,1,0,1,0,nullptr,off);
    capped.advance(stopped,1,0,1,0,nullptr,p);
    const double freed=(1-.05)*-std::expm1(-1.0/5);
    double total_area=0;
    for (size_t i=0;i<f.cells.size();++i) {
        EXPECT_DOUBLE_EQ(legacy.cells()[i].continental_area,1);
        EXPECT_NEAR(capped.cells()[i].continental_area,1-freed,1e-12);
        EXPECT_NEAR(capped.cells()[i].oceanic_area,freed,1e-12);
        EXPECT_NEAR(capped.cells()[i].continent,1000,1e-9);
        EXPECT_DOUBLE_EQ(capped.cells()[i].ocean,0);
        total_area+=f.band.grid->area_m2(static_cast<CellId>(i));
    }
    EXPECT_NEAR(capped.budget().oceanized_continental_area_m2,freed*total_area,1e-9*total_area);
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
    f.cells.assign(f.band.cellcount(),tectonics::columns::seed(100,true,40));
    f.move(-.15,.15);
    tectonics::margins::Parameters legacy_parameters;legacy_parameters.enabled=true;legacy_parameters.river_erosion=0;
    // The legacy share is the minimum thickness 0 written out explicitly.
    auto explicit_zero=legacy_parameters;explicit_zero.continental_min_thickness_m=0;
    const auto p=growth_parameters();
    Model legacy,zero,limited;
    for (auto* model:{&legacy,&zero,&limited}) model->initialize(f.band.grid,f.cells);
    const auto thinnest=[&](const Model& model) {
        double least=std::numeric_limits<double>::infinity();
        for (const auto& c:model.cells()) if (c.continental_area>.05) least=std::min(least,c.continent/c.continental_area);
        return least;
    };
    double limited_least=std::numeric_limits<double>::infinity();
    for (int step=0;step<100;++step) {
        legacy.advance(f.motion,1,0,1,0,nullptr,legacy_parameters);
        zero.advance(f.motion,1,0,1,0,nullptr,explicit_zero);
        limited.advance(f.motion,1,0,1,0,nullptr,p);
        limited_least=std::min(limited_least,thinnest(limited));
    }
    for (size_t i=0;i<f.cells.size();++i) {
        EXPECT_EQ(legacy.cells()[i].continent,zero.cells()[i].continent);
        EXPECT_EQ(legacy.cells()[i].continental_area,zero.cells()[i].continental_area);
    }
    // Legacy rifts stretch the crust far below 20 km before it breaks up; with the
    // minimum, the opening is ocean once the crust reaches it; the cap relaxes the
    // overshoot of a fast opening (0.3 columns/Myr here) within its 5 Myr.
    EXPECT_LT(thinnest(legacy),15000);
    EXPECT_GT(limited_least,.7*p.continental_min_thickness_m);
    RecordProperty("thinnest_limited_m",std::to_string(thinnest(limited)));
    EXPECT_GT(thinnest(limited),.9*p.continental_min_thickness_m);
    EXPECT_GT(limited.budget().generated_ocean_kg,legacy.budget().generated_ocean_kg);
    EXPECT_LT(std::abs(limited.budget().relative_mass_residual),1e-11);
    EXPECT_LT(std::abs(limited.budget().felsic_relative_residual),1e-11);
}

TEST(TectonicMargins, ArcContinentShareKeepsMassAndBuildsIslandArcArea) {
    MarginFixture f;
    auto p=growth_parameters();p.arc_continent_fraction=.5;
    f.move(.6,0);
    Model model;model.initialize(f.band.grid,f.cells);
    tectonics::margins::Polarity polarity;
    for (int step=0;step<200;++step) {
        const auto forcing=f.forcing(model.cells(),p,&polarity);
        model.advance(f.motion,1,0,1,0,&forcing,p);
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
    Model model;model.initialize(f.band.grid,f.cells);
    tectonics::margins::Polarity polarity;
    for (int step=0;step<100;++step) {
        const auto forcing=f.forcing(model.cells(),p,&polarity);
        model.advance(f.motion,1,0,1,0,&forcing,p);
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
        Model model;model.initialize(f.band.grid,f.cells);
        tectonics::margins::Polarity polarity;
        for (int step=0;step<steps;++step) {
            const auto forcing=f.forcing(model.cells(),p,&polarity);
            model.advance(f.motion,1,0,1,0,&forcing,p);
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
    Model model;model.initialize(f.band.grid,f.cells);
    tectonics::margins::Polarity polarity;
    const auto step=[&] {
        const auto forcing=f.forcing(model.cells(),p,&polarity);
        model.advance(f.motion,1,0,1,0,&forcing,p);
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
    // Erosion turns continent into sediment: the felsic reservoir keeps its mass (above a
    // slowly falling target, no sinks) while the continent falls short. Only the
    // continent-only reservoir sees the deficit and makes arc melt continent.
    const auto eroding=[](bool continent_only) {
        ContinentalMarginFixture f;
        auto p=growth_parameters();
        p.continental_growth_control=true;p.continental_growth_per_gyr=-.01;
        p.continental_growth_continent_only=continent_only;
        Model model;model.initialize(f.band.grid,f.cells);
        tectonics::margins::Polarity polarity;
        std::vector<double> shares;
        for (int step=0;step<150;++step) {
            const auto forcing=f.forcing(model.cells(),p,&polarity);
            // Erosion 0.25: the explicit diffusion limit on the cube is 0.2 over the
            // non-orthogonal correction factor (50 km columns).
            model.advance(f.motion,1,.25,1,0,&forcing,p);
            shares.push_back(model.budget().arc_continent_fraction);
        }
        EXPECT_LT(std::abs(model.budget().felsic_relative_residual),1e-11);
        return std::pair{model.budget(),shares};
    };
    const auto [felsic,felsic_shares]=eroding(false);
    const auto [continent,continent_shares]=eroding(true);
    for (double share:felsic_shares) EXPECT_DOUBLE_EQ(share,0);
    EXPECT_DOUBLE_EQ(felsic.generated_continent_kg,0);
    EXPECT_GE(felsic.felsic_kg,felsic.initial_felsic_kg*(1-1e-12));
    EXPECT_LT(felsic.continent_kg,felsic.initial_continent_kg);
    EXPECT_GT(continent.generated_continent_kg,0);
    // The share rises as the continent's deficit grows.
    RecordProperty("continent_only_share_150",std::to_string(continent_shares.back()));
    EXPECT_GT(continent_shares.back(),continent_shares[10]);
    EXPECT_GT(continent.continent_kg,felsic.continent_kg);
}

namespace {
// Uniform mixed columns: no motion, no erosion and no root flow between equal cells.
std::vector<tectonics::columns::Cell> mixed_cells(size_t count,double area,double continent,double ocean) {
    tectonics::columns::Cell c;
    c.continental_area=area;c.oceanic_area=1-area;c.continent=continent;c.ocean=ocean;
    c.continental_age_volume=continent*100;c.oceanic_age_volume=ocean*20;
    return std::vector<tectonics::columns::Cell>(count,c);
}
// A strip of columns 20..23 at `strip` (area 1) on the margin fixture's ocean, optionally an
// interior cell (area 1, 30 km) at column 19 and a margin cell (area 0.5, 12 km) at column 24.
std::vector<tectonics::columns::Cell> strip_cells(const MarginFixture& f,double strip,bool neighbours) {
    auto c=f.cells;
    for (int y=0;y<f.band.rows();++y) for (int x=0;x<f.band.columns();++x) {
        const CellId i=f.band.cell(x,y);
        if (x>=20 && x<24) c[i]=mixed_cells(1,1,strip,0)[0];
        else if (neighbours && x==19) c[i]=mixed_cells(1,1,30000,0)[0];
        else if (neighbours && x==24) c[i]=mixed_cells(1,.5,12000,3500)[0];
    }
    return c;
}
// The strip's rows away from its ends: band rows 2..13. There the strip is uniform along its
// length and the corrected differences across its sides are the plain ones; its end rows (0
// and 15) also feed the polar faces through the non-orthogonal correction.
constexpr int inner_first=2, inner_last=13;
// The continent (m x m^2) columns first..last lost in the inner rows.
double inner_loss(const MarginFixture& f,const Model& model,double initial,int first=20,int last=23) {
    double lost=0;
    for (int y=inner_first;y<=inner_last;++y) for (int x=first;x<=last;++x) {
        const CellId i=f.band.cell(x,y);
        lost+=(initial-model.cells()[i].continent)*f.band.grid->area_m2(i);
    }
    return lost;
}
}

TEST(TectonicMargins, ContinentalMaximumThicknessGainsAreaWithoutMovingRock) {
    MarginFixture f;
    const auto stopped=cubefixture::still(*f.band.grid);
    const auto cells=mixed_cells(f.cells.size(),.5,30000,3500);
    Model reference,capped;
    reference.initialize(f.band.grid,cells);capped.initialize(f.band.grid,cells);
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
        sum+=(c.continental_area-.5)*f.band.grid->area_m2(static_cast<CellId>(i));
    }
    EXPECT_NEAR(capped.budget().continentalized_oceanic_area_m2,sum,1e-12*sum);
    EXPECT_DOUBLE_EQ(reference.budget().continentalized_oceanic_area_m2,0);
    for (int step=1;step<200;++step) capped.advance(stopped,1,0,1,0,nullptr,p);
    for (const auto& c:capped.cells()) {
        EXPECT_NEAR(c.continental_area,.75,1e-6);
        EXPECT_NEAR(c.continent/c.continental_area,40000,1);
    }
    EXPECT_LT(std::abs(capped.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(capped.budget().felsic_relative_residual),1e-12);
    // Room for more area than the oceanic area: the column saturates at area 1.
    Model saturated;
    saturated.initialize(f.band.grid,mixed_cells(f.cells.size(),.9,50000,350));
    for (int step=0;step<200;++step) saturated.advance(stopped,1,0,1,0,nullptr,p);
    // The cap moves no ocean; as the oceanic area vanishes, a rounding-sized overlap is a large
    // share of it and subducts (budgeted) a little of the ocean volume.
    double lost=0;
    for (CellId c=0;c<saturated.cells().size();++c) {
        const auto& cell=saturated.cells()[c];
        EXPECT_NEAR(cell.continental_area,1,1e-6);
        EXPECT_GE(cell.oceanic_area,0);
        EXPECT_LE(cell.ocean,350);
        EXPECT_GT(cell.ocean,349);
        lost+=(350-cell.ocean)*f.band.grid->area_m2(c);
    }
    EXPECT_NEAR(saturated.budget().subducted_ocean_kg-saturated.budget().generated_ocean_kg,lost*2900,
                1e-6*350*2900*f.band.grid->totalarea_m2());
    EXPECT_LT(std::abs(saturated.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(saturated.budget().felsic_relative_residual),1e-12);
}

TEST(TectonicMargins, ContinentalThicknessBetweenTheCapsKeepsItsArea) {
    MarginFixture f;
    const auto stopped=cubefixture::still(*f.band.grid);
    Model model;
    model.initialize(f.band.grid,mixed_cells(f.cells.size(),.5,15000,3500));
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
    for (int y=0;y<f.band.rows();++y) for (int x=0;x<f.band.columns();++x) {
        const CellId i=f.band.cell(x,y);
        if (x>=20 && x<24) f.cells[i]=mixed_cells(1,1,60000,0)[0];
        else if (x==24) f.cells[i]=mixed_cells(1,.3,10000,4900)[0];
    }
    const auto stopped=cubefixture::still(*f.band.grid);
    auto p=growth_parameters();p.continental_max_thickness_m=40000;p.orogenic_collapse_km2_myr=3000;
    p.root_mobility=0; // Collapse alone; root flow never moves volume into pure ocean.
    Model model;model.initialize(f.band.grid,f.cells);
    const auto [area0,continent0]=std::pair{f.total(model,&tectonics::columns::Cell::continental_area),
                                            f.total(model,&tectonics::columns::Cell::continent)};
    for (int step=0;step<50;++step) {
        model.advance(stopped,1,0,1,0,nullptr,p);
        for (int y=0;y<f.band.rows();++y) for (int x=20;x<24;++x)
            EXPECT_GE(model.cells()[f.band.cell(x,y)].continent,p.continental_max_thickness_m*(1-1e-12));
    }
    const double area1=f.total(model,&tectonics::columns::Cell::continental_area);
    const double continent1=f.total(model,&tectonics::columns::Cell::continent);
    EXPECT_GT(area1,area0);
    EXPECT_NEAR(continent1,continent0,1e-12*continent0);
    EXPECT_LT(model.cells()[f.at(20)].continent,60000);
    EXPECT_LT(model.cells()[f.at(23)].continent,60000);
    // The ocean west of the strip gained continent of the strip's age.
    const auto& west=model.cells()[f.at(19)],&source=model.cells()[f.at(21)];
    ASSERT_GT(west.continent,0);
    EXPECT_GT(west.continental_area,0);
    EXPECT_NEAR(west.continental_age_volume/west.continent,source.continental_age_volume/source.continent,1e-9);
    for (const auto& c:model.cells()) EXPECT_FALSE(c.continental_area>=.5 && c.continent<1000);
    EXPECT_GT(model.budget().collapsed_continent_kg,0);
    EXPECT_GT(model.budget().continentalized_oceanic_area_m2,0);
    EXPECT_LT(std::abs(model.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(model.budget().felsic_relative_residual),1e-12);
}

// The collapse onset (orogenic_collapse_onset_thickness_m): crust between the maximum
// thickness and the onset is held by its strength; above the onset only the excess over
// it spreads; an onset of 0 or equal to the maximum is the legacy collapse bit for bit.
TEST(TectonicMargins, OrogenicCollapseStartsAtItsOnsetThickness) {
    MarginFixture f;
    const auto strip=[&](double continent) { return strip_cells(f,continent,false); };
    const auto stopped=cubefixture::still(*f.band.grid);
    auto p=growth_parameters();p.continental_max_thickness_m=40000;p.orogenic_collapse_km2_myr=3000;
    p.root_mobility=0; // Collapse alone.
    const auto run=[&](Model& model,const std::vector<tectonics::columns::Cell>& cells,
                       const tectonics::margins::Parameters& q,int steps) {
        model.initialize(f.band.grid,cells);
        for (int step=0;step<steps;++step) model.advance(stopped,1,0,1,0,nullptr,q);
    };
    const auto same_cells=[&](const Model& a,const Model& b) {
        for (size_t i=0;i<f.cells.size();++i) {
            const auto& x=a.cells()[i];const auto& y=b.cells()[i];
            EXPECT_EQ(x.continent,y.continent) << i;EXPECT_EQ(x.continental_area,y.continental_area) << i;
            EXPECT_EQ(x.continental_age_volume,y.continental_age_volume) << i;EXPECT_EQ(x.ocean,y.ocean) << i;
        }
    };
    auto onset=p;onset.orogenic_collapse_onset_thickness_m=60000;
    auto no_collapse=p;no_collapse.orogenic_collapse_km2_myr=0;
    // 55 km: above the maximum, below the onset: held, the run without collapse bit for bit.
    {
        Model held,still,legacy;
        run(held,strip(55000),onset,50);run(still,strip(55000),no_collapse,50);run(legacy,strip(55000),p,50);
        same_cells(held,still);
        EXPECT_DOUBLE_EQ(held.budget().collapsed_continent_kg,0);
        EXPECT_GT(legacy.budget().collapsed_continent_kg,0);
        EXPECT_LT(legacy.cells()[f.at(20)].continent,55000);
    }
    // 70 km: only the excess over the onset spreads. One sub-step: each of the strip's two
    // edges per row passes kappa (70 - 60) km against kappa (70 - 40) km (the corrected
    // differences of a strip along meridian columns are nearly the plain ones).
    {
        auto slow=p;slow.orogenic_collapse_km2_myr=25;
        auto slow_onset=onset;slow_onset.orogenic_collapse_km2_myr=25;
        Model a,b;
        run(a,strip(70000),slow_onset,1);run(b,strip(70000),slow,1);
        ASSERT_GT(b.budget().collapsed_continent_kg,0);
        // Measured 0.3437: the corrected differences on the cube's face read the strip's
        // neighbours across rows too (3 per cent).
        RecordProperty("onset_over_legacy",std::to_string(inner_loss(f,a,70000)/inner_loss(f,b,70000)));
        EXPECT_NEAR(inner_loss(f,a,70000)/inner_loss(f,b,70000),1./3,2e-2);
        Model spread;
        for (int step=0;step<50;++step) {
            if (step==0) run(spread,strip(70000),onset,1);
            else spread.advance(stopped,1,0,1,0,nullptr,onset);
            for (int y=0;y<f.band.rows();++y) for (int x=20;x<24;++x)
                EXPECT_GE(spread.cells()[f.band.cell(x,y)].continent,60000*(1-1e-12));
        }
        EXPECT_LT(spread.cells()[f.at(20)].continent,70000);
        EXPECT_GT(spread.cells()[f.at(19)].continent,0);
        EXPECT_LT(std::abs(spread.budget().relative_mass_residual),1e-12);
        EXPECT_LT(std::abs(spread.budget().felsic_relative_residual),1e-12);
    }
    // Onset 0 and onset at the maximum are the legacy collapse.
    {
        auto at_maximum=p;at_maximum.orogenic_collapse_onset_thickness_m=40000;
        auto zero=p;zero.orogenic_collapse_onset_thickness_m=0;
        Model legacy,a,b;
        run(legacy,strip(70000),p,20);run(a,strip(70000),zero,20);run(b,strip(70000),at_maximum,20);
        same_cells(legacy,a);same_cells(legacy,b);
        EXPECT_EQ(legacy.budget().collapsed_continent_kg,b.budget().collapsed_continent_kg);
    }
    // Validation: 0, or at least the maximum thickness (which must be set).
    EXPECT_NO_THROW(tectonics::margins::validate(onset));
    for (double value : {30000.0,39999.0,-1.0,std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::infinity()}) {
        auto invalid=p;invalid.orogenic_collapse_onset_thickness_m=value;
        EXPECT_THROW(tectonics::margins::validate(invalid),std::runtime_error) << value;
    }
    auto no_maximum=growth_parameters();no_maximum.orogenic_collapse_onset_thickness_m=60000;
    EXPECT_THROW(tectonics::margins::validate(no_maximum),std::runtime_error);
}

// With an onset: the margin cap (orogenic_collapse_margin_cap) lets thick crust spread into
// a receiver with oceanic area above the maximum thickness and holds it against interior
// neighbours until the onset; the creep below the onset (orogenic_collapse_below_onset_km2_myr)
// adds D_low (min(c_a, onset) - max(c_b, cap))+ to the full D (c_a - max(c_b, onset))+.
// One sub-step (D 25 km²/Myr), so the fluxes scale with the coefficients and thresholds.
TEST(TectonicMargins, OrogenicCollapseMarginCapAndCreepBelowTheOnset) {
    MarginFixture f;
    const auto cells=[&](double strip,bool neighbours) { return strip_cells(f,strip,neighbours); };
    auto legacy=growth_parameters();legacy.continental_max_thickness_m=40000;legacy.orogenic_collapse_km2_myr=25;
    legacy.root_mobility=0;
    auto onset=legacy;onset.orogenic_collapse_onset_thickness_m=65000;
    const auto stopped=cubefixture::still(*f.band.grid);
    const auto collapsed=[&](const std::vector<tectonics::columns::Cell>& start,const tectonics::margins::Parameters& q,
                             Model& model) {
        model.initialize(f.band.grid,start);
        model.advance(stopped,1,0,1,0,nullptr,q);
        return model.budget().collapsed_continent_kg;
    };
    Model m0,m1,m2,m3,m4;
    // A: a 50 km strip loses toward the margin cell (50 - 40 km), not toward the interior
    // cell (held until 65 km): about half of the legacy collapse, which also feeds the interior.
    auto margin=onset;margin.orogenic_collapse_margin_cap=true;
    const double legacy_kg=collapsed(cells(50000,true),legacy,m0),margin_kg=collapsed(cells(50000,true),margin,m1);
    ASSERT_GT(legacy_kg,0);
    ASSERT_GT(margin_kg,0);
    EXPECT_NEAR(inner_loss(f,m1,50000)/inner_loss(f,m0,50000),.5,2e-2);
    for (int y=inner_first;y<=inner_last;++y) {
        EXPECT_EQ(m1.cells()[f.band.cell(19,y)].continent,30000);
        EXPECT_EQ(m1.cells()[f.band.cell(20,y)].continent,50000);
        EXPECT_LT(m1.cells()[f.band.cell(23,y)].continent,50000);
        EXPECT_GT(m1.cells()[f.band.cell(24,y)].continent,12000);
        EXPECT_GT(m0.cells()[f.band.cell(19,y)].continent,30000);
    }
    // Without any margin receiver above the cap the margin cap is the onset alone, bit for bit.
    Model held_a,held_b;
    const auto interior=[&] { auto c=cells(70000,false);
        for (auto& cell:c) if (cell.continental_area<1) cell=mixed_cells(1,1,25000,0)[0];
        return c; }();
    collapsed(interior,margin,held_a);collapsed(interior,onset,held_b);
    for (size_t i=0;i<f.cells.size();++i) EXPECT_EQ(held_a.cells()[i].continent,held_b.cells()[i].continent) << i;
    EXPECT_EQ(held_a.budget().collapsed_continent_kg,held_b.budget().collapsed_continent_kg);
    // B: a 70 km strip beside ocean: 0.4 x 25 km + 5 km against 5 km (D_low 10 of D 25).
    auto creep=onset;creep.orogenic_collapse_below_onset_km2_myr=10;
    const double onset_kg=collapsed(cells(70000,false),onset,m2),creep_kg=collapsed(cells(70000,false),creep,m3);
    ASSERT_GT(onset_kg,0);
    ASSERT_GT(creep_kg,0);
    EXPECT_NEAR(inner_loss(f,m3,70000)/inner_loss(f,m2,70000),3,.1);
    // At 50 km (below the onset) only the creep acts: 0.4 x 10 km against legacy 10 km.
    Model c_a,c_b;
    EXPECT_NEAR(collapsed(cells(50000,false),creep,c_a)/collapsed(cells(50000,false),legacy,c_b),.4,1e-2);
    // A and B: the margin receiver takes the full collapse above the cap (10 km), the
    // interior receiver the creep (0.4 x (50 - 40) km): 14 against legacy 20 per row.
    auto both=margin;both.orogenic_collapse_below_onset_km2_myr=10;
    EXPECT_NEAR(collapsed(cells(50000,true),both,m4)/legacy_kg,.7,2e-2);
    EXPECT_GT(m4.cells()[f.at(19)].continent,30000);
    EXPECT_LT(std::abs(m4.budget().relative_mass_residual),1e-12);
    // Validation: both options require an onset; the creep must be finite and not negative.
    EXPECT_NO_THROW(tectonics::margins::validate(both));
    auto no_onset=legacy;no_onset.orogenic_collapse_margin_cap=true;
    EXPECT_THROW(tectonics::margins::validate(no_onset),std::runtime_error);
    no_onset=legacy;no_onset.orogenic_collapse_below_onset_km2_myr=10;
    EXPECT_THROW(tectonics::margins::validate(no_onset),std::runtime_error);
    // The onset and both options require the collapse itself.
    for (auto invalid : {onset,margin,creep,both}) {
        invalid.orogenic_collapse_km2_myr=0;
        EXPECT_THROW(tectonics::margins::validate(invalid),std::runtime_error);
    }
    for (double value : {-1.0,std::numeric_limits<double>::quiet_NaN()}) {
        auto invalid=creep;invalid.orogenic_collapse_below_onset_km2_myr=value;
        EXPECT_THROW(tectonics::margins::validate(invalid),std::runtime_error) << value;
    }
    // The creep is weak: at most the collapse diffusivity.
    auto equal=creep;equal.orogenic_collapse_below_onset_km2_myr=equal.orogenic_collapse_km2_myr;
    EXPECT_NO_THROW(tectonics::margins::validate(equal));
    auto stronger=creep;stronger.orogenic_collapse_below_onset_km2_myr=stronger.orogenic_collapse_km2_myr*1.01;
    EXPECT_THROW(tectonics::margins::validate(stronger),std::runtime_error);
}

// The margin cap's receiver alone decides: a margin cell (oceanic area 0.1) above the cap
// feeds its ocean neighbour from the cap but its interior neighbour only above the onset;
// two sub-steps with both coefficients close the mass budget.
TEST(TectonicMargins, OrogenicCollapseReceiverDecidesAndTwoRegimesSubStep) {
    MarginFixture f;
    const auto stopped=cubefixture::still(*f.band.grid);
    auto legacy=growth_parameters();legacy.continental_max_thickness_m=40000;legacy.orogenic_collapse_km2_myr=25;
    legacy.root_mobility=0;
    auto margin=legacy;margin.orogenic_collapse_onset_thickness_m=65000;margin.orogenic_collapse_margin_cap=true;
    // Column 20: a margin cell (area 0.9, 50 km of cell-mean continent); 21: an interior cell
    // (area 1, 30 km); ocean elsewhere (19 and 22 included).
    auto cells=f.cells;
    for (int y=0;y<f.band.rows();++y) {
        cells[f.band.cell(20,y)]=mixed_cells(1,.9,50000,700)[0];
        cells[f.band.cell(21,y)]=mixed_cells(1,1,30000,0)[0];
    }
    Model a,b;
    a.initialize(f.band.grid,cells);b.initialize(f.band.grid,cells);
    a.advance(stopped,1,0,1,0,nullptr,margin);b.advance(stopped,1,0,1,0,nullptr,legacy);
    ASSERT_GT(b.budget().collapsed_continent_kg,0);
    EXPECT_NEAR(inner_loss(f,a,50000,20,20)/inner_loss(f,b,50000,20,20),.5,2e-2);
    for (int y=inner_first;y<=inner_last;++y) {
        // Held (the transport's total-restoring scale moves it by rounding only).
        EXPECT_NEAR(a.cells()[f.band.cell(21,y)].continent,30000,1e-9);
        EXPECT_GT(a.cells()[f.band.cell(19,y)].continent,0);
        EXPECT_GT(b.cells()[f.band.cell(21,y)].continent,30000);
    }
    // Two regimes over two sub-steps, a 70 km strip (columns 20 to 23) beside ocean: the inner
    // columns lose only after the edge columns have dropped below them.
    auto creep=legacy;creep.orogenic_collapse_onset_thickness_m=65000;
    creep.orogenic_collapse_km2_myr=400;creep.orogenic_collapse_below_onset_km2_myr=160;
    Model two;two.initialize(f.band.grid,strip_cells(f,70000,false));
    two.advance(stopped,1,0,1,0,nullptr,creep);
    EXPECT_GT(two.budget().collapsed_continent_kg,0);
    for (int y=0;y<f.band.rows();++y) {
        EXPECT_LT(two.cells()[f.band.cell(20,y)].continent,70000);
        EXPECT_LT(two.cells()[f.band.cell(20,y)].continent,two.cells()[f.band.cell(21,y)].continent);
        EXPECT_GT(two.cells()[f.band.cell(19,y)].continent,0);
    }
    EXPECT_LT(std::abs(two.budget().relative_mass_residual),1e-12);
    EXPECT_LT(std::abs(two.budget().felsic_relative_residual),1e-12);
}

TEST(TectonicMargins, ZeroMaximumThicknessAndCollapseAreTheLegacyEngine) {
    ContinentalMarginFixture f;
    auto p=growth_parameters();p.sediment_subduction_fraction=.75;p.subduction_erosion_ratio=.07;
    p.continental_growth_control=true;
    auto explicit_zero=p;explicit_zero.continental_max_thickness_m=0;explicit_zero.orogenic_collapse_km2_myr=0;
    Model legacy,zero;
    legacy.initialize(f.band.grid,f.cells);zero.initialize(f.band.grid,f.cells);
    tectonics::margins::Polarity legacy_polarity,zero_polarity;
    for (int step=0;step<100;++step) {
        const auto a=f.forcing(legacy.cells(),p,&legacy_polarity);
        const auto b=f.forcing(zero.cells(),explicit_zero,&zero_polarity);
        legacy.advance(f.motion,1,0,1,0,&a,p);
        zero.advance(f.motion,1,0,1,0,&b,explicit_zero);
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

// The arc walk on the cube (W2-E, E2B): an ocean plate on face 0 closes on a continent that
// holds the face's last two columns and every other face; the walk from the overriding
// boundary cell goes straight on across the face edge, by the cell table's ahead.
TEST(TectonicMarginsCube, ArcWalkCrossesAFaceEdgeStraightOn) {
    constexpr int n = 32;
    constexpr double earth_km = 6371.0;
    const auto grid = tectonics::grid::TectonicGrid::cube(n, earth_km);
    const tectonics::columns::Metric metric(*grid);
    const auto& mesh = *grid->mesh();
    const size_t count = grid->cellcount();
    std::vector<uint32_t> owners(count, 1);
    std::vector<tectonics::columns::Cell> cells(count, tectonics::columns::seed(500, true, 100));
    std::vector<mesh::Vec3> rigid(count);
    const mesh::Vec3 right = mesh::cubeface(0).right;
    for (CellId c = 0; c < count; ++c) {
        int face = 0, column = 0, row = 0;
        mesh::cubeindices(n, c, face, column, row);
        if (face != 0 || column >= n - 2) continue;
        owners[c] = 0;
        cells[c] = tectonics::columns::seed(-4000, false, 60);
        // 30 km per update toward the face's east edge.
        const mesh::Vec3 r = mesh.centre(c);
        rigid[c] = 30000.0 * mesh::normalise(right - mesh::dot(right, r) * r);
    }
    tectonics::margins::Parameters p;
    p.enabled = true;
    p.arc_distance_km = 700;
    const auto forcing = tectonics::margins::forcing(metric, 1, owners, cells, rigid, p, nullptr, 1, 0);
    const CellId lower = mesh::cubeid(n, 0, n - 3, n / 2), upper = mesh::cubeid(n, 0, n - 2, n / 2);
    const CellId receiver = forcing.arc_receiver[lower];
    ASSERT_NE(receiver, std::numeric_limits<uint32_t>::max());
    EXPECT_EQ(forcing.arc_receiver[upper], receiver);
    EXPECT_GT(forcing.trench_m[lower], 0);
    int face = 0, column = 0, row = 0;
    mesh::cubeindices(n, receiver, face, column, row);
    EXPECT_NE(face, 0);
    EXPECT_EQ(owners[receiver], 1u);
    // Three steps of about 313 km (the walk stops at the first beyond 700 km), straight on:
    // the receiver lies on the great circle through the trench and its overriding cell.
    const double spacing_km = mesh.distance_m(upper, 3) / 1000;
    const double walked_km = mesh::greatcircle_rad(mesh.centre(upper), mesh.centre(receiver)) * earth_km;
    EXPECT_NEAR(walked_km, 3 * spacing_km, .05 * spacing_km);
    const mesh::Vec3 pole = mesh::normalise(mesh::cross(mesh.centre(lower), mesh.centre(upper)));
    EXPECT_LT(std::abs(mesh::dot(pole, mesh.centre(receiver))) * earth_km, .15 * spacing_km);
}
