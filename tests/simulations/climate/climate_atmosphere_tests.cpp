#include "climate_atmosphere.hpp"
#include "climate_weather.hpp"
#include "detail/surface_wind_remap.hpp"
#include "detail/surface_drag_mobility.hpp"
#include "detail/thermal_trough.hpp"
#include "detail/surface_momentum_advection.hpp"
#include "detail/surface_850_exchange.hpp"
#include "detail/moisture_carrier.hpp"
#include "detail/katabatic_drainage.hpp"
#include "detail/upper_thermal_projection.hpp"
#include "detail/mountain_wave_height.hpp"
#include "detail/stationary_replay.hpp"
#include "detail/seasonal_overturning.hpp"
#include "detail/frozen_surface_drag.hpp"
#include "detail/vegetation_drag.hpp"
#include "../../../src/io/climate_surface_input.hpp"
#include <sstream>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <chrono>
#include <iomanip>
#include <string>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (condition)
        return;

    std::cerr << "FAIL: " << message << '\n';
    failures++;
}

void testSurface850Exchange()
{
    using namespace climateatmosphere::detail;
    expect(lowLevelMoistureWeight(1.0f, 1.0f, 0.0f, 290.0f, 9.81f) == 1.0f &&
        lowLevelMoistureWeight(1.0f, 1.0f, 2500.0f, 290.0f, 9.81f) == 0.0f &&
        lowLevelMoistureWeight(1.0f, 0.0f, 0.0f, 290.0f, 9.81f) == 0.0f,
        "850-hPa moisture carriers must exclude buried and unavailable pressure levels");
    expect(moistureCarrierMode(-8.0f, 12.0f, 0.0f) == -8.0f &&
        moistureCarrierMode(-8.0f, 12.0f, 1.0f) == 12.0f &&
        moistureCarrierMode(-8.0f, 12.0f, .25f) == -3.0f,
        "carrier interpolation must preserve signed vectors and both layer endpoints");
    const SurfaceMomentumForce fs{.00004f,.00008f}, fa{.0001f,.0003f};
    constexpr float cd=.0013f, depth=300, density=1.225f, mass=1162, rayleigh=1.0f/432000;
    for (float f : {-1e-4f, 0.0f, 1e-4f})
        for (float exchange : {0.0f,5e-6f,2e-5f,1e-3f})
        {
            const auto wind=coupleSurface850(fs,fa,f,cd,depth,density,mass,rayleigh,exchange);
            expect(wind.valid && wind.maximumResidualMps2 < 1e-7,
                "coupled surface/850 winds must satisfy both momentum equations in both hemispheres");
            const double us=wind.surface.eastMetresPerSecond, ns=-wind.surface.southMetresPerSecond;
            const double ua=wind.lowLevel.eastMetresPerSecond, na=-wind.lowLevel.southMetresPerSecond;
            const double ms=depth*density;
            const double forcing=ms*(fs.eastMps2*us+fs.northMps2*ns)+mass*(fa.eastMps2*ua+fa.northMps2*na);
            const double drag=ms*cd/depth*std::pow(std::hypot(us,ns),3)+mass*rayleigh*(ua*ua+na*na);
            expect(wind.mixingLossWm2 >= 0 && std::abs(forcing-drag-wind.mixingLossWm2) < 1e-5*std::max(1.0,std::abs(forcing)),
                "equal/opposite layer stress must dissipate shear and close the local steady work budget");
            const double totalMixingWork=wind.stressEastNm2*(us-ua)+wind.stressNorthNm2*(ns-na);
            expect(std::abs(totalMixingWork+wind.mixingLossWm2) < 1e-8,
                "surface stress and its low-level reaction must use the same interface force");
            if (exchange == 0)
            {
                const auto isolated=climateatmosphere::steadyMixedDragCoriolisWind(fs.eastMps2,fs.northMps2,f,0,cd/depth);
                expect(wind.surface.eastMetresPerSecond == isolated.eastMetresPerSecond &&
                    wind.surface.southMetresPerSecond == isolated.southMetresPerSecond && wind.mixingLossWm2 == 0,
                    "disabled exchange must recover the isolated surface balance");
            }
        }
    expect(!coupleSurface850(fs,fa,0,cd,depth,density,-1,rayleigh,1e-5f).valid,
        "invalid reservoir mass must reject momentum exchange");
}

void testSurfaceMomentumAdvection()
{
    using climateatmosphere::detail::advectSurfaceMomentum;
    constexpr int columns = 64, rows = 32, cells = columns * rows;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 6371000.0);
    climateatmosphere::ModeSeparatedCirculationConfig config;
    config.rotationRatePerSecond = 0;
    config.surfaceBoundaryLayerDepthMetres = 300;
    std::vector<float> zero(cells, 0), cd(cells, .01f), north(cells, .01f * 9 / 300), south(cells, -3);
    const auto calm = advectSurfaceMomentum(grid, zero, zero, cd, zero, zero, config, 64);
    expect(calm.east == zero && calm.south == zero, "advection must preserve a calm unforced state");
    const auto uniform = advectSurfaceMomentum(grid, zero, north, cd, zero, south, config, 64);
    expect(uniform.east == zero && uniform.south == south && uniform.relativeMomentumResidual < 1e-6,
        "uniform meridional drag balance must survive horizontal momentum transport");
    std::vector<float> forcing(cells, .00015f), initial(cells);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const int i = y * columns + x;
            cd[i] = x < columns / 2 ? .0013f : .05f;
            initial[i] = std::sqrt(forcing[i] * 300 / cd[i]);
        }
    const auto flow = advectSurfaceMomentum(grid, forcing, zero, cd, initial, zero, config, 128);
    expect(flow.east.size() == cells && flow.relativeMomentumResidual < .005,
        "the stationary upwind roughness-step solve must converge to its momentum equation");
    expect(flow.dragWorkWm2 >= 0 && flow.upwindMixingWorkWm2 >= 0 &&
        std::abs(flow.forcingWorkWm2 - flow.dragWorkWm2 - flow.advectionWorkWm2) < .005 * flow.forcingWorkWm2,
        "force work must balance positive drag and diagnosed advective work within the solve tolerance");
    const int ocean = rows / 2 * columns;
    expect(flow.east[ocean] < flow.east[ocean + 5] && flow.east[ocean + columns / 2 + 5] < flow.east[ocean],
        "ocean wind must recover over downwind fetch after a rough continent");
    std::vector<float> shiftedCd(cells), shiftedInitial(cells), du(cells), dv(cells);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto i = grid.index(x, y), j = grid.index(x + 7, y);
            shiftedCd[j] = cd[i]; shiftedInitial[j] = initial[i];
            du[i] = flow.east[i] - initial[i]; dv[i] = flow.south[i];
            if (std::abs(grid.latitudeCentresRadians[y]) > 35 * 3.141592653589793 / 180)
                expect(flow.east[i] == initial[i] && flow.south[i] == 0,
                    "tropical advection must retain the outside-domain wind exactly");
        }
    const auto shifted = advectSurfaceMomentum(grid, forcing, zero, shiftedCd, shiftedInitial, zero, config, 128);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto i = grid.index(x, y), j = grid.index(x + 7, y);
            expect(shifted.east[j] == flow.east[i] && shifted.south[j] == flow.south[i],
                "momentum transport must be invariant under a longitude-seam rotation");
        }
    const auto ascent = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(grid, du, dv, 150);
    expect(std::abs(climategrid::areaWeightedIntegral(columns, rows, ascent)) /
        (4 * 3.141592653589793) < 1e-6,
        "advection ascent correction must preserve the global pressure-mass integral");
    cd[0] = -1;
    expect(advectSurfaceMomentum(grid, forcing, zero, cd, initial, zero, config, 128).east.empty(),
        "invalid momentum coefficients must reject the diagnostic solve");
}

void testFrozenSurfaceDrag()
{
    using climateatmosphere::detail::vegetationSurfaceDragCoefficient;
    for (float land : {0.0f, .3f, 1.0f})
        for (float height : {0.0f, 2000.0f, 5000.0f})
            for (float ratio : {.25f, .5f, .9f, 1.0f, 2.0f, 4.0f})
            {
                const auto drag = [&](float ice, float strength) {
                    return vegetationSurfaceDragCoefficient(land,height,ice,.0013f,.018f,.036f,.003f,ratio,strength);
                };
                const float old = climateatmosphere::detail::frozenSurfaceDragCoefficient(land,height,0,.0013f,.018f,.036f,.003f);
                expect(drag(0,0) == old, "disabled vegetation contrast must preserve the old drag exactly");
                expect(drag(0,1) > 0 && std::isfinite(drag(0,1)), "vegetation stress must remain finite and dissipative");
                expect(drag(land,1) == drag(land,0), "fully frozen land and pure ocean must be unaffected");
                expect(drag(0,.5f) == old + land*.018f*.5f*(ratio-1),
                    "the previous half-strength calibration must remain exactly reproducible");
                for (float strong : {6.0f, 8.0f, 12.0f, 16.0f, 100.0f})
                {
                    expect(drag(0,strong)>0 && std::isfinite(drag(0,strong)),
                        "amplified contrast must remain dissipative, including smooth surfaces");
                    expect(drag(land,strong)==drag(land,0),
                        "strong vegetation stress must not alter pure ocean or frozen land");
                    if (ratio>1 && land>0)
                        expect(drag(0,strong)>drag(0,1), "large gains must actually amplify vegetation drag");
                    if (ratio<1)
                        expect(drag(0,strong)==drag(0,1),
                            "protection gains must not collapse all exposed surfaces to the drag floor");
                }
                expect(drag(0,100)==drag(0,16), "experimental vegetation gain has a finite upper bound");
                if (ratio>1 && land>0)
                {
                    expect(drag(0,12)>drag(0,6), "the new stress test must be stronger than run 294");
                    expect(drag(0,6)==old+land*.018f*6*(ratio-1),
                        "run 294 remains reproducible with gain six and its increase-only input");
                }
                expect(std::abs(drag(0,1)-old-land*.018f*(ratio-1)) < 1e-7f,
                    "vegetation correction must be area weighted without replacing relief drag");
            }
    std::istringstream good("UWROUGH1 4 2 1 2 .3 1 1 1 1 1");
    expect(climateio::readSurfaceDragContrast(good,4,2).size()==8, "explicit surface input must read the registered grid");
    for (const char* bad : {"BAD 4 2", "UWROUGH1 8 4", "UWROUGH1 4 2 1", "UWROUGH1 4 2 -1 1 1 1 1 1 1 1",
        "UWROUGH1 4 2 1 1 1 1 1 1 1 1 extra", "UWROUGH1 4 2 nan 1 1 1 1 1 1 1"})
    {
        bool rejected=false;
        try { std::istringstream stream(bad); climateio::readSurfaceDragContrast(stream,4,2); }
        catch (const std::runtime_error&) { rejected=true; }
        expect(rejected,"malformed, wrong-grid and nonfinite surface inputs must fail explicitly");
    }
    using climateatmosphere::detail::frozenSurfaceDragCoefficient;
    using climateatmosphere::detail::persistentSnowCover;
    const std::vector<float> winter{1, 1, 0, .5f}, spring{.8f, .7f, 1, .6f},
        summer{.4f, 0, 1, .7f}, autumn{.7f, .1f, 1, .8f};
    std::array<const std::vector<float>*, 4> quarters{&winter, &spring, &summer, &autumn};
    expect(persistentSnowCover(quarters, 4) == std::vector<float>({.4f, 0, 0, .5f}),
        "only cover retained in the same cell through all seasons may reduce persistent-snow drag");
    quarters[2] = nullptr;
    expect(persistentSnowCover(quarters, 4) == std::vector<float>(4, 0),
        "missing seasonal evidence must not manufacture persistent snow");
    quarters[2] = &summer;
    expect(persistentSnowCover(quarters, 3) == std::vector<float>(3, 0),
        "mismatched seasonal grids must not manufacture persistent snow");
    for (float land : {0.0f, 0.25f, 0.7f, 1.0f})
        for (float height : {0.0f, 2000.0f, 7000.0f})
        {
            const float base = .0013f + land * (.018f - .0013f);
            const float legacy = base + std::clamp(height / 4500.0f, 0.0f, 1.0f) * (.036f - base);
            expect(frozenSurfaceDragCoefficient(land, height, 0, .0013f, .018f, .036f, .003f) == legacy,
                "ice-free drag must preserve the previous closure exactly");
            if (land == 0) continue;
            const float frozen = .0013f + land * (.003f - .0013f);
            expect(std::abs(frozenSurfaceDragCoefficient(land, height, land, .0013f, .018f, .036f, .003f) - frozen) < 1e-8f,
                "fully frozen land must use its own drag without changing ocean area");
            const float partial = frozenSurfaceDragCoefficient(land, height, land / 2, .0013f, .018f, .036f, .003f);
            expect(partial >= std::min(legacy, frozen) && partial <= std::max(legacy, frozen) && partial > 0,
                "partial cover must blend positive dissipative stress coefficients");
        }
    expect(frozenSurfaceDragCoefficient(0, 0, 1, .0013f, .018f, .036f, .003f) == .0013f,
        "land ice drag must not alter pure ocean cells");
}

void testKatabaticDrainage()
{
    using namespace climateatmosphere::detail;
    for (double depth : {50.0, 75.0, 100.0, 150.0})
    {
        double integral = 0;
        for (int i = 0; i < 10000; ++i) integral += coldJetProfileWeight((i + .5) * depth / 10000, depth) / 10000;
        expect(std::abs(integral - coldJetMeanToTenMetres(depth)) < 1e-7,
            "10 m jet profile must integrate to the mean used by momentum");
        expect(std::abs(integral * depth / 300 - coldJetTransportFraction(depth, 300)) < 1e-7,
            "column transport must equal the same jet's integrated volume flux");
        expect(coldJetProfileWeight(10, depth) == 1 && coldJetProfileWeight(0, depth) == 0 &&
            coldJetProfileWeight(depth, depth) == 0 && coldJetProfileWeight(depth/2, depth) >= 1,
            "jet increment must be normalized at 10 m and vanish at ground/top");
        for (float latitude : {-85.0f, -45.0f, 0.0f, 45.0f, 85.0f})
            for (float cd : {.001f, .003f, .012f})
                for (float rotation : {-1.0f, 1.0f})
                {
                    const climateatmosphere::HorizontalWind carrier{3, -2};
                    const SurfaceMomentumForce cold{.001f, -.0007f};
                    const auto unchanged = coldJetTenMetreWind(carrier, {}, latitude, cd, depth, 7.2921159e-5f, rotation);
                    expect(unchanged.eastMetresPerSecond == 3 && unchanged.southMetresPerSecond == -2,
                        "zero cold forcing must preserve the background exactly");
                    const auto wind = coldJetTenMetreWind(carrier, cold, latitude, cd, depth, 7.2921159e-5f, rotation);
                    const double u=wind.eastMetresPerSecond, v=-wind.southMetresPerSecond, u0=3, v0=2;
                    const double f=climateatmosphere::coriolisParameterPerSecond(latitude,7.2921159e-5f,rotation);
                    const double r=cd/depth, mean=coldJetMeanToTenMetres(depth);
                    const double east=r*(std::hypot(u,v)*u-std::hypot(u0,v0)*u0)-f*mean*(v-v0)-cold.eastMps2;
                    const double north=r*(std::hypot(u,v)*v-std::hypot(u0,v0)*v0)+f*mean*(u-u0)-cold.northMps2;
                    expect(std::hypot(east,north)<1e-8,
                        "shallow incremental force balance must close in either hemisphere/rotation");
                    expect(cold.eastMps2*(u-u0)+cold.northMps2*(v-v0)>0,
                        "incremental quadratic stress must dissipate cold forcing work");
                }
    }
    expect(coldJetMeanToTenMetres(10)==0 && coldJetTransportFraction(100,75)==0,
        "invalid or uncontained cold layers must be rejected");
    for (int columns : {64,128,256})
    {
        const auto grid=climategrid::makeSphericalGrid(columns,columns/2,6371000);
        std::vector<float> dome(columns*columns/2);
        for (int y=0;y<grid.rows;++y) for (int x=0;x<columns;++x)
            dome[grid.index(x,y)]=3000*std::sin(grid.latitudeCentresRadians[y]);
        const auto old=terrainSlopes(grid,dome), fixed=terrainSlopes(grid,dome,true);
        for (int y : {0,grid.rows-1})
        {
            const double exact=3000*std::cos(grid.latitudeCentresRadians[y])/6371000;
            expect(std::abs(fixed[grid.index(0,y)].north/exact-1)<.003,
                "pole-crossing centred scalar gradient must recover a smooth dome's analytic slope");
            expect(std::abs(old[grid.index(0,y)].north/exact-1)>.9,
                "diagnostic fixture must expose the historical doubled polar slope");
        }
    }
    for (float control : {-6.0f, 0.0f, 6.0f})
        for (float latitude = -90; latitude <= 90; latitude += .25f)
        {
            const float original = climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(latitude, control, 30, 10);
            expect(confinedOverturningPressureHpa(latitude, control, 30, 0, control, 10) == original,
                "disabled confined migration must preserve every pressure value");
            const float changed = confinedOverturningPressureHpa(latitude, control, 30, 0, control * 2.5f, 10);
            expect(std::isfinite(changed) && changed >= -8.501f && changed <= 10.001f,
                "confined migration must preserve bounded finite pressure");
            if (std::abs(latitude - control) >= 30)
                expect(changed == original, "confined migration must leave extratropical forcing unchanged");
            expect(std::abs(changed - confinedOverturningPressureHpa(-latitude, -control, 30, 0, -control * 2.5f, 10)) < 1e-5f,
                "confined migration must reflect exactly between hemispheres");
        }
    for (float target : {-15.0f, 0.0f, 15.0f})
    {
        expect(std::abs(confinedOverturningPressureHpa(target, 6, 30, 0, target, 10) + 5.5f) < 1e-5f,
            "confined trough must reach the intended minimum");
        for (float edge : {-24.0f, 36.0f})
            expect(std::abs(confinedOverturningPressureHpa(edge, 6, 30, 0, target, 10) - 10) < 1e-5f,
                "subtropical maxima must remain fixed");
    }
    for (float latitude = -90; latitude <= 90; latitude += .25f)
    {
        const float p = climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(latitude,0,30,10);
        expect(polarBranchPressureHpa(p,latitude,0,30,0,10,.5f,0)==p,
            "disabled polar thermal contrast must preserve every forcing value");
        const float q=polarBranchPressureHpa(p,latitude,0,30,0,10,.5f,1);
        if (std::abs(latitude)<=69) expect(q==p,"polar thermal contrast must preserve all values through the subpolar low");
        expect(std::abs(q-polarBranchPressureHpa(p,-latitude,0,30,0,10,.5f,1))<1e-6f,
            "polar thermal closure must reflect with its thermal inputs");
        expect(std::isfinite(q),"polar pressure must remain finite");
    }
    expect(std::abs(polarBranchPressureHpa(4.5f,90,0,30,0,10,.5f,1)+2.0f)<1e-6f,
        "half polar contrast must halve the high-minus-low pressure difference");
    for (float latitude : {-80.0f, -60.0f, 0.0f, 60.0f, 80.0f})
        for (float land : {0.0f, 0.3f, 1.0f})
        {
            const float delta = continentalSubpolarCorrectionHpa(-12, latitude, 0, 30, 0, land, 1);
            expect(delta >= 0 && delta <= 12 && std::isfinite(delta),
                "marine weighting must not reverse the negative eddy lobe or deepen land lows");
            expect(delta == continentalSubpolarCorrectionHpa(-12, -latitude, 0, 30, 0, land, 1),
                "marine weighting must reflect between hemispheres");
            if (std::abs(latitude) <= 30 || land == 0)
                expect(delta == 0, "marine weighting must preserve tropical troughs and all-ocean forcing");
            else expect(std::abs(delta - 12 * land) < 1e-6f,
                "coastal fractional land must partition the eddy forcing continuously");
            expect(continentalSubpolarCorrectionHpa(12, latitude, 0, 30, 0, land, 1) == 0 &&
                continentalSubpolarCorrectionHpa(-12, latitude, 0, 30, 0, land, 0) == 0,
                "positive lobes and disabled marine weighting must retain the original template");
        }
    for (float latitude : {-90.0f, -65.0f, -20.0f, 0.0f, 20.0f, 65.0f, 90.0f})
        for (float ratio : {0.25f, 0.5f, 1.0f, 1.5f, 2.0f})
        {
            const float p = hemisphericSubpolarPressureHpa(-12, latitude, 5, 30, .5f, ratio, 2);
            expect(std::isfinite(p) && p <= 0 && p >= -48,
                "hemispheric scaling must retain finite bounded negative pressure lobes");
            expect(p == hemisphericSubpolarPressureHpa(-12, -latitude, -5, 30, .5f, ratio, 2),
                "hemispheric scaling must reflect with the seasonal trough");
            expect(hemisphericSubpolarPressureHpa(-12, latitude, 5, 30, .5f, ratio, 0) == -12 &&
                hemisphericSubpolarPressureHpa(-12, latitude, 5, 30, .5f, 1, 2) == -12 &&
                hemisphericSubpolarPressureHpa(12, latitude, 5, 30, .5f, ratio, 2) == 12,
                "disabled scaling, equal contrasts and positive pressure must preserve the template");
            if (std::abs(latitude) <= 20)
                expect(p == -12, "hemispheric eddy scaling must leave tropical forcing unchanged");
        }
    expect(hemisphericSubpolarPressureHpa(-12, 65, 0, 30, 0, .5f, 2) == -3 &&
        hemisphericSubpolarPressureHpa(-12, 65, 0, 30, 0, 1.5f, 2) == -27,
        "hemispheric eddy lobes must weaken or deepen with their own thermal contrast");
    expect(snowPatchKatabaticDeficitK(5, 4, 1, .5f, 15) == 2.5f &&
        snowPatchKatabaticDeficitK(5, 4, .2f, .5f, 15) == .5f,
        "remaining melting snow must have a bounded cold-patch deficit weighted by cover");
    expect(snowPatchKatabaticDeficitK(5, 4, 0, .5f, 15) == 0 &&
        snowPatchKatabaticDeficitK(-5, 0, 1, .5f, 15) == 0,
        "bare ground and air colder than the snow patch must not gain drainage");
    expect(snowPatchKatabaticDeficitK(-10, -20, 1, .5f, 15) ==
        katabaticDeficitK(-10, -20, 1, .5f, 15),
        "fully frozen cold-skin conditions must preserve the original deficit");
    expect(radiativeKatabaticDeficitK(-10, -10, 1, 1, 1000, 100, 10000, 1, 15) == 1,
        "radiative inversion must convert joules per square metre to kelvin using layer heat capacity");
    expect(radiativeKatabaticDeficitK(-10, -10, .5f, 1, 1000, 100, 10000, 1, 15) == .5f,
        "radiative inversion must scale with frozen area");
    expect(radiativeKatabaticDeficitK(10, -10, 1, 1, 1000, 100, 10000, 1, 15) == 0 &&
        radiativeKatabaticDeficitK(-10, -10, 0, 1, 1000, 100, 10000, 1, 15) == 0 &&
        radiativeKatabaticDeficitK(-10, -10, 1, 1, 1000, 100, 10000, 0, 15) == 0,
        "warming, absent snow and disabled cooling must not create an inversion");
    expect(radiativeKatabaticDeficitK(-500, -10, 1, 1, 1000, 100, 10000, 1, 3) == 3 &&
        radiativeKatabaticDeficitK(-10, -10, 1, 1, 1000, 0, 10000, 1, 15) == 0,
        "inversion estimate must be bounded and reject zero layer capacity");
    for (float trough : {-15.0f, 0.0f, 15.0f})
        for (float latitude : {-70.0f, -10.0f, 10.0f, 70.0f})
        {
            expect(winterAnchoredHadleyWidth(latitude, trough, 30, 0) == 30,
                "disabled winter anchoring must preserve the original width");
            expect(winterAnchoredHadleyWidth(latitude, trough, 30, 1) ==
                winterAnchoredHadleyWidth(-latitude, -trough, 30, 1),
                "seasonal geometry must reflect between hemispheres");
        }
    for (float trough : {-15.0f, 15.0f})
    {
        const float edge = trough > 0 ? -30.0f : 30.0f;
        const float width = winterAnchoredHadleyWidth(edge, trough, 30, 1);
        expect(std::abs(climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(edge, trough, width, 10) - 10) < 1e-6,
            "winter pressure maximum must remain at the diagnosed subtropical edge");
        expect(winterAnchoredHadleyWidth(-edge, trough, 30, 1) == 30,
            "winter anchoring must leave the summer branch width unchanged");
    }
    for (double decay : {0.01, 25.0, 75.0, 300.0, 1e9})
    {
        constexpr double depth = 300.0;
        const double mean = katabaticLayerMeanFraction(decay, depth);
        double integrated = 0.0;
        constexpr int samples = 100000;
        for (int i = 0; i < samples; ++i)
            integrated += katabaticProfileWeight((i + .5) * depth / samples, decay, depth) / samples;
        expect(std::abs(integrated - mean) < 1e-6 && mean >= 0 && mean <= .5,
            "bulk drainage velocity must equal the integrated, bounded vertical profile");
        expect(katabaticProfileWeight(0, decay, depth) == 1.0 &&
            katabaticProfileWeight(depth, decay, depth) == 0.0 &&
            katabaticProfileWeight(1000, decay, depth) == 0.0,
            "surface drainage must be retained locally and vanish above the boundary layer");
    }
    expect(katabaticLayerMeanFraction(0, 300) == 0 && katabaticLayerMeanFraction(75, 0) == 0,
        "zero thickness profiles must not introduce a transport increment");
    expect(std::abs(katabaticLayerMeanFraction(75, 300) - .23134264f) < 1e-7,
        "the trial profile must retain its documented mass-weighted velocity fraction");
    expect(std::abs(frozenLandAlbedo(0, 0, .2f, .55f, .8f) - .2f) < 1e-6f &&
        std::abs(frozenLandAlbedo(1, 0, .2f, .55f, .8f) - .55f) < 1e-6f &&
        std::abs(frozenLandAlbedo(1, 1, .2f, .55f, .8f) - .8f) < 1e-6f,
        "frozen surface mixtures must preserve bare, snow and permanent ice endpoints");
    expect(std::abs(frozenLandAlbedo(.75f, .5f, .2f, .55f, .8f) - .5875f) < 1e-6f,
        "snow and diagnosed ice overlap must not count the frozen area twice");
    for (float snow : {-1.0f, 0.0f, .5f, 1.0f, 2.0f})
        for (float ice : {-1.0f, 0.0f, .5f, 1.0f, 2.0f})
        {
            const float albedo = frozenLandAlbedo(snow, ice, .2f, .55f, .8f);
            expect(albedo >= .2f && albedo <= .800001f,
                "fractional frozen surface albedo must remain between material endpoints");
        }
    expect(terrainSlopes(climategrid::makeSphericalGrid(2, 1, 6371000), {1, 1}).empty(),
        "a one-row grid cannot support meridional terrain gradients");
    expect(katabaticDeficitK(-10, -20, 1, 0, 15) == 0, "disabled drainage must exert no force");
    expect(katabaticDeficitK(-10, -20, 0, 0.5f, 15) == 0, "snow-free terrain must not receive ice drainage");
    expect(katabaticDeficitK(-20, -10, 1, 0.5f, 15) == 0, "warm skin must not drive cold drainage");
    expect(katabaticDeficitK(5, 1, 1, 0.5f, 15) == 0, "thawed surfaces must not drive ice drainage");
    expect(katabaticDeficitK(-10, -20, 1, 0.5f, 15) == 5, "half skin deficit must yield half layer deficit");
    expect(katabaticDeficitK(0, -100, 1, 1, 15) == 15, "cold-layer deficit must stay bounded");
    const auto flat = katabaticForce(0, 0, 250, 10, 9.81f, 0.05f);
    expect(flat.eastMps2 == 0 && flat.northMps2 == 0, "flat ice must not generate drainage");
    const auto force = katabaticForce(0.01f, 0.02f, 250, 10, 9.81f, 0.05f);
    expect(force.eastMps2 < 0 && force.northMps2 < 0, "buoyancy force must point downhill");
    const auto capped = katabaticForce(1, 1, 250, 10, 9.81f, 0.05f);
    expect(std::hypot(capped.eastMps2, capped.northMps2) <= 9.81f*10/250*0.050001f,
        "unresolved steep slopes must not produce unbounded acceleration");
    for (float latitude : {-70.0f, 0.0f, 70.0f})
    {
        const auto wind = climateatmosphere::steadyQuadraticDragCoriolisWind(
            force.eastMps2, force.northMps2, latitude, 0.018f, 300, 7.2921159e-5f, 1);
        expect(force.eastMps2 * wind.eastMetresPerSecond - force.northMps2 * wind.southMetresPerSecond > 0,
            "drag/Coriolis balance must retain a downhill component in either hemisphere");
    }
    for (int columns : {32, 64})
    {
        const auto grid = climategrid::makeSphericalGrid(columns, columns/2, 6371000);
        std::vector<float> height(columns*columns/2);
        for (int y = 0; y < grid.rows; ++y)
            for (int x = 0; x < columns; ++x)
                height[grid.index(x,y)] = static_cast<float>(10000 + 0.001 * grid.radiusMetres * grid.latitudeCentresRadians[y]);
        const auto slopes = terrainSlopes(grid, height);
        for (const auto& slope : slopes)
            expect(std::abs(slope.north - 0.001f) < 1e-8f && slope.east == 0,
                "metric planar slope must be independent of resolution and finite at polar rows");
        std::fill(height.begin(), height.end(), 3000.0f);
        const auto level = terrainSlopes(grid, height);
        for (const auto& slope : level)
            expect(slope.east == 0 && slope.north == 0, "constant elevation must have no slope including at the seam");
    }
}

void testSubpolarDeepening()
{
    const auto deepen = climateatmosphere::detail::deepenSubpolarLowPressureHpa;
    for (int latitude = -90; latitude <= 90; ++latitude)
    {
        const float p = deepen(7.0f, static_cast<float>(latitude), 5.0f, 30.0f, 0.0f, 10.0f, 0.5f, false);
        const float reflected = deepen(7.0f, static_cast<float>(-latitude), -5.0f, 30.0f, 0.0f, 10.0f, 0.5f, false);
        expect(std::isfinite(p) && p >= 2.0f && p <= 7.0f && std::abs(p - reflected) < 1e-6f,
            "subpolar deepening must be bounded and symmetric when hemisphere and seasonal trough reflect");
        if ((latitude >= -25 && latitude <= 35) || std::abs(latitude) == 90)
            expect(p == 7.0f, "subpolar deepening must preserve tropical forcing and both polar endpoints");
        expect(deepen(7.0f, static_cast<float>(latitude), 5.0f, 30.0f, 0.0f, 10.0f, 0.0f, false) == 7.0f,
            "disabled subpolar deepening must exactly preserve baseline pressure");
    }
    expect(std::abs(deepen(7.0f, 69.0f, 0.0f, 30.0f, 0.0f, 10.0f, 0.5f, false) - 2.0f) < 1e-6f,
        "additional pressure depression must have its specified amplitude at the existing subpolar low");
    for (int latitude = 69; latitude <= 90; ++latitude)
    {
        const float original = 0.2f * latitude;
        const float changed = deepen(original, static_cast<float>(latitude), 0.0f, 30.0f, 0.0f, 10.0f, 0.5f, true);
        expect(std::abs(changed - original + 5.0f) < 1e-6f,
            "polar-gradient preservation must add the same pressure offset everywhere poleward of the low");
        expect(deepen(original, static_cast<float>(-latitude), 0.0f, 30.0f, 0.0f, 10.0f, 0.5f, true) == changed,
            "the preserved polar branch must reflect symmetrically between hemispheres");
    }
}

void testUpperThermalProjection()
{
    const auto grid = climategrid::makeSphericalGrid(32, 16, 6371000.0);
    std::vector<float> source(16);
    double mean = 0.0, area = 0.0;
    for (int y = 0; y < 16; ++y)
    {
        source[y] = 250.0f + 40.0f * std::cos(grid.latitudeCentresRadians[y]);
        mean += source[y] * grid.cellAreasSquareMetres[y];
        area += grid.cellAreasSquareMetres[y];
    }
    mean /= area;
    for (float scale : {0.0f, 0.35f, 0.5f, 1.0f})
    {
        const auto projected = climateatmosphere::detail::projectUpperThermalGradient(
            source, grid.cellAreasSquareMetres, scale, scale, 20.0f, 50.0f);
        for (int y = 0; y < 16; ++y)
            expect(std::abs(projected[y] - scale * (source[y] - mean)) < 3e-6,
                "uniform gradient projection must recover uniformly scaled thermal anomalies");
    }
    const auto projected = climateatmosphere::detail::projectUpperThermalGradient(
        source, grid.cellAreasSquareMetres, 0.35f, 1.0f, 20.0f, 50.0f);
    double projectedMean = 0.0;
    for (int y = 0; y < 16; ++y)
    {
        projectedMean += projected[y] * grid.cellAreasSquareMetres[y] / area;
        if (y > 0 && source[y] != source[y - 1])
        {
            const double ratio = (projected[y] - projected[y - 1]) / (source[y] - source[y - 1]);
            expect(ratio >= 0.34999 && ratio <= 1.00001,
                "latitude blending must neither reverse nor amplify a thermal gradient");
        }
    }
    expect(std::abs(projectedMean) < 1e-6, "projected height forcing must have zero area mean");
    const auto flat = climateatmosphere::detail::projectUpperThermalGradient(
        std::vector<float>(16, 280.0f), grid.cellAreasSquareMetres, 0.35f, 1.0f, 20.0f, 50.0f);
    expect(std::all_of(flat.begin(), flat.end(), [](float x) { return x == 0.0f; }),
        "latitude blending must not create thermal winds from uniform temperature");

    const auto polar = climateatmosphere::detail::projectUpperThermalGradient(
        source, grid.cellAreasSquareMetres, 0.70f, 1.0f, 20.0f, 50.0f, 0.25f, 60.0f, 80.0f);
    double polarMean = 0.0;
    for (int y = 0; y < 16; ++y)
    {
        polarMean += polar[y] * grid.cellAreasSquareMetres[y] / area;
        expect(std::abs(polar[y] - polar[15 - y]) < 3e-6,
            "polar thermal projection must preserve hemispheric reflection");
        if (y > 0 && source[y] != source[y - 1])
        {
            const double ratio = (polar[y] - polar[y - 1]) / (source[y] - source[y - 1]);
            expect(ratio >= 0.24999 && ratio <= 1.00001,
                "polar attenuation must neither reverse nor amplify the source gradient");
            const double latitude = std::abs(90.0 - 180.0 * y / 16.0);
            if (latitude >= 60.0)
                expect(ratio < 1.0, "polar projection must attenuate high-latitude thermal gradients");
            if (latitude >= 50.0 && latitude <= 60.0)
                expect(std::abs(ratio - 1.0) < 1e-5,
                    "polar attenuation must preserve the unmodified midlatitude gradient");
        }
    }
    expect(std::abs(polarMean) < 1e-6, "polar thermal projection must retain zero area mean");
    const auto polarFlat = climateatmosphere::detail::projectUpperThermalGradient(
        std::vector<float>(16, 280.0f), grid.cellAreasSquareMetres, 0.7f, 1.0f,
        20.0f, 50.0f, 0.25f, 60.0f, 80.0f);
    expect(std::all_of(polarFlat.begin(), polarFlat.end(), [](float x) { return x == 0.0f; }),
        "polar attenuation must not create winds from uniform temperature");
}

void testReceiverSurfaceDrag()
{
    climateatmosphere::ModeSeparatedCirculationConfig config;
    const climateatmosphere::HorizontalWind original{8.0f, -3.0f};
    constexpr float oceanCd = 0.0013f;
    for (float latitude : {-60.0f, -15.0f, 0.0f, 15.0f, 60.0f})
    {
        const auto same = climateatmosphere::detail::rediagnoseSurfaceWindForDrag(
            original, oceanCd, oceanCd, latitude, config);
        expect(same.eastMetresPerSecond == original.eastMetresPerSecond &&
            same.southMetresPerSecond == original.southMetresPerSecond,
            "receiver drag must preserve the exact wind when source and recipient drag match");
        const auto calm = climateatmosphere::detail::rediagnoseSurfaceWindForDrag(
            {}, oceanCd, 0.003f, latitude, config);
        expect(calm.eastMetresPerSecond == 0.0f && calm.southMetresPerSecond == 0.0f,
            "receiver roughness must not accelerate a calm flow, including at the equator");
        const double f = climateatmosphere::coriolisParameterPerSecond(latitude,
            config.rotationRatePerSecond, config.rotationDirection);
        const double originalDrag = oceanCd * std::hypot(original.eastMetresPerSecond,
            original.southMetresPerSecond) / config.surfaceBoundaryLayerDepthMetres;
        const double eastForce = originalDrag * original.eastMetresPerSecond + f * original.southMetresPerSecond;
        const double northForce = f * original.eastMetresPerSecond - originalDrag * original.southMetresPerSecond;
        double previousSpeed = std::hypot(original.eastMetresPerSecond, original.southMetresPerSecond);
        for (float cd : {0.003f, 0.006f})
        {
            const auto rough = climateatmosphere::detail::rediagnoseSurfaceWindForDrag(
                original, oceanCd, cd, latitude, config);
            const double speed = std::hypot(rough.eastMetresPerSecond, rough.southMetresPerSecond);
            const double drag = cd * speed / config.surfaceBoundaryLayerDepthMetres;
            expect(speed < previousSpeed,
                "stronger receiver drag must monotonically slow fixed-force flow in either hemisphere");
            expect(std::hypot(drag * rough.eastMetresPerSecond + f * rough.southMetresPerSecond - eastForce,
                f * rough.eastMetresPerSecond - drag * rough.southMetresPerSecond - northForce) < 2.0e-9,
                "receiver roughness must retain the effective momentum forcing in either hemisphere");
            previousSpeed = speed;
        }
    }
    const auto equatorial = climateatmosphere::detail::rediagnoseSurfaceWindForDrag(
        original, oceanCd, 0.003f, 0.0f, config);
    const float scale = std::sqrt(oceanCd / 0.003f);
    expect(std::abs(equatorial.eastMetresPerSecond - scale * original.eastMetresPerSecond) < 2.0e-6f &&
        std::abs(equatorial.southMetresPerSecond - scale * original.southMetresPerSecond) < 2.0e-6f,
        "equatorial fixed-force flow must obey the quadratic-drag speed ratio without changing direction");
}

void testConstantForceAcrossCoast()
{
    climateatmosphere::ModeSeparatedCirculationConfig config;
    constexpr float east = 0.0008f, north = -0.0003f;
    for (float latitude : {-60.0f, -15.0f, 0.0f, 15.0f, 60.0f})
    {
        climateatmosphere::detail::SurfaceMomentumForce source[2];
        for (int side = 0; side < 2; ++side)
        {
            const float cd = side ? 0.018f : 0.0013f;
            const auto wind = climateatmosphere::steadyQuadraticDragCoriolisWind(
                east, north, latitude, cd, config.surfaceBoundaryLayerDepthMetres,
                config.rotationRatePerSecond, config.rotationDirection);
            source[side] = climateatmosphere::detail::surfaceMomentumForce(wind, cd, latitude, config);
        }
        for (float fraction : {0.0f, 0.25f, 0.5f, 0.75f, 1.0f})
            for (float cd : {0.0013f, 0.006f, 0.018f})
            {
                const float blendedEast = source[0].eastMps2 * (1.0f - fraction) + source[1].eastMps2 * fraction;
                const float blendedNorth = source[0].northMps2 * (1.0f - fraction) + source[1].northMps2 * fraction;
                const auto actual = climateatmosphere::steadyQuadraticDragCoriolisWind(
                    blendedEast, blendedNorth, latitude, cd, config.surfaceBoundaryLayerDepthMetres,
                    config.rotationRatePerSecond, config.rotationDirection);
                const auto expected = climateatmosphere::steadyQuadraticDragCoriolisWind(
                    east, north, latitude, cd, config.surfaceBoundaryLayerDepthMetres,
                    config.rotationRatePerSecond, config.rotationDirection);
                expect(std::hypot(actual.eastMetresPerSecond - expected.eastMetresPerSecond,
                    actual.southMetresPerSecond - expected.southMetresPerSecond) < 5e-6,
                    "constant native forcing must not acquire an artificial coastal acceleration on refinement");
            }
    }
}

void testParallelSurfaceDrag()
{
    constexpr int columns = 32, rows = 16;
    std::vector<float> drag(columns * rows);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            drag[y * columns + x] = x % 2 ? 0.1f : 0.0013f;
    const auto effective = climateatmosphere::detail::parallelSurfaceDragCoefficients(
        columns, rows, drag, columns / 2, rows / 2);
    for (float east : {-0.0008f, 0.0f, 0.0008f})
        for (float north : {-0.0003f, 0.0f, 0.0003f})
        {
            const auto wind = [&](float cd) {
                return climateatmosphere::steadyQuadraticDragCoriolisWind(
                    east, north, 0.0f, cd, 300.0f, 7.2921159e-5f, 1.0f);
            };
            const auto ocean = wind(0.0013f), land = wind(0.1f);
            for (float cd : effective)
            {
                const auto bulk = wind(cd);
                expect(std::hypot(bulk.eastMetresPerSecond - 0.5f *
                    (ocean.eastMetresPerSecond + land.eastMetresPerSecond),
                    bulk.southMetresPerSecond - 0.5f *
                    (ocean.southMetresPerSecond + land.southMetresPerSecond)) < 4e-6f,
                    "mixed-cell mobility must reproduce area-mean patch winds under signed common equatorial forcing");
            }
        }
    std::fill(drag.begin(), drag.end(), 0.018f);
    const auto uniform = climateatmosphere::detail::parallelSurfaceDragCoefficients(
        columns, rows, drag, columns / 2, rows / 2);
    expect(std::all_of(uniform.begin(), uniform.end(), [](float cd) { return std::abs(cd - 0.018f) < 1e-8f; }),
        "uniform roughness must be invariant under mobility aggregation");
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            drag[y * columns + x] *= 1.0f + 0.01f * x;
    auto rotated = drag;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            rotated[y * columns + (x + 2) % columns] = drag[y * columns + x];
    const auto original = climateatmosphere::detail::parallelSurfaceDragCoefficients(columns, rows, drag, columns / 2, rows / 2);
    const auto shifted = climateatmosphere::detail::parallelSurfaceDragCoefficients(columns, rows, rotated, columns / 2, rows / 2);
    for (int y = 0; y < rows / 2; ++y)
        for (int x = 0; x < columns / 2; ++x)
            expect(std::abs(original[y * (columns / 2) + x] - shifted[y * (columns / 2) + (x + 1) % (columns / 2)]) < 1e-8f,
                "mixed roughness must rotate through the longitude seam with the surface");
    drag.front() = 0.0f;
    expect(climateatmosphere::detail::parallelSurfaceDragCoefficients(columns, rows, drag, columns / 2, rows / 2).empty(),
        "nonpositive drag must be rejected before computing mobility");
}

void testRegionalThermalTrough()
{
    constexpr int columns = 64, rows = 32;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 6371000.0);
    std::vector<float> temperature(columns * rows, 300.0f);
    const auto diagnose = [&](const auto& field, float origin) {
        return climateatmosphere::detail::regionalThermalEquatorDegrees(
            columns, rows, field, grid.latitudeCentresRadians, origin, 30.0f);
    };
    const auto flat = diagnose(temperature, -4.0f);
    expect(std::all_of(flat.begin(), flat.end(), [](float x) { return x == -4.0f; }),
        "uniform temperature must not create a regional trough shift");
    for (int sign : {-1, 1})
    {
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const double latitude = grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793;
                temperature[y * columns + x] = static_cast<float>(300.0 - 0.05 * std::pow(latitude - sign * 24.0, 2));
            }
        const auto trough = diagnose(temperature, sign * 15.0f);
        expect(std::all_of(trough.begin(), trough.end(), [sign](float x) { return std::abs(x - sign * 24.0f) < 0.001f; }),
            "continental thermal peak must be recovered beyond the old 0.65-width search in either hemisphere");
    }
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const double latitude = grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793;
            const double peak = 5.0 + 2.0 * std::sin(2.0 * 3.141592653589793 * x / columns);
            temperature[y * columns + x] = static_cast<float>(300.0 - 0.05 * std::pow(latitude - peak, 2));
        }
    auto rotated = temperature;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
            rotated[y * columns + (x + 17) % columns] = temperature[y * columns + x];
    const auto trough = diagnose(temperature, -4.0f), shifted = diagnose(rotated, -4.0f);
    for (int x = 0; x < columns; ++x)
        expect(trough[x] >= -14.5f && trough[x] <= 6.5f && trough[x] == shifted[(x + 17) % columns],
            "regional thermal trough must remain bounded and rotate with the geography");
}

void testReceiverAscentCorrection()
{
    constexpr int columns = 32, rows = 16, cells = columns * rows;
    constexpr double pi = 3.14159265358979323846;
    constexpr float depth = 7.5f;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, 6371000.0);
    std::vector<float> east(cells), south(cells), shiftedEast(cells), shiftedSouth(cells);
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            east[grid.index(x, y)] = static_cast<float>(std::sin(2.0 * pi * x / columns) * (1.0 + y));
            south[grid.index(x, y)] = static_cast<float>(std::cos(4.0 * pi * x / columns) + 0.1 * y);
        }
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            shiftedEast[grid.index(x + 1, y)] = east[grid.index(x, y)];
            shiftedSouth[grid.index(x + 1, y)] = south[grid.index(x, y)];
        }
    const auto correction = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(grid, east, south, depth);
    const auto shifted = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(grid, shiftedEast, shiftedSouth, depth);
    double integral = 0.0, magnitude = 0.0;
    for (int y = 0; y < rows; ++y)
        for (int x = 0; x < columns; ++x)
        {
            const auto cell = grid.index(x, y);
            integral += grid.cellAreasSquareMetres[y] * correction[cell];
            magnitude += grid.cellAreasSquareMetres[y] * std::abs(correction[cell]);
            expect(shifted[grid.index(x + 1, y)] == correction[cell],
                "ascent corrections must remain periodic and unchanged by a longitude seam shift");
        }
    expect(std::abs(integral) < 1.0e-7 * magnitude,
        "receiver wind corrections must not introduce globally integrated ascent");
    const std::vector<float> zero(cells, 0.0f), uniformSouth(cells, 3.0f);
    const auto unchanged = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(grid, zero, zero, depth);
    expect(std::all_of(unchanged.begin(), unchanged.end(), [](float value) { return value == 0.0f; }),
        "unchanged winds must leave the legacy ascent field exactly unchanged");
    const auto polar = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(grid, zero, uniformSouth, depth);
    const double northExpected = -86400.0 * depth * 3.0 * grid.southFaceLengthsMetres.front() / grid.cellAreasSquareMetres.front();
    const double southExpected = 86400.0 * depth * 3.0 * grid.northFaceLengthsMetres.back() / grid.cellAreasSquareMetres.back();
    expect(std::abs(polar.front() - northExpected) < 1.0e-5 &&
        std::abs(polar.back() - southExpected) < 1.0e-5,
        "nonzero polar-row winds must still have closed outer polar faces");
    climateatmosphere::ModeSeparatedCirculationConfig config;
    config.enabled = {true, false, true, false};
    climateatmosphere::ModeSeparatedCirculation flow;
    flow.surfacePressureAnomalyHpa = east;
    flow.upperHeightAnomalyMetres = zero;
    climateatmosphere::diagnoseModeWinds(columns, rows, config, flow);
    const auto lowerAscent = climateatmosphere::detail::surfaceAscentCorrectionHpaPerDay(grid,
        flow.surfaceEastWindMps, flow.surfaceSouthWindMps, 0.5f * config.surfaceEquivalentPressureDepthHpa);
    for (int cell = 0; cell < cells; ++cell)
        expect(std::abs(lowerAscent[cell] - flow.ascentHpaPerDay[cell]) < 1.0e-5f,
            "receiver corrections must use the legacy atmosphere's spherical mass-flux divergence operator");
}
}

void testStationaryPreconditioning()
{
    using Solver = climateatmosphere::StationarySolver;
    constexpr double pi = 3.14159265358979323846;
    for (int columns : {12, 32})
        for (float direction : {-1.0f, 1.0f})
            for (bool rowProjection : {false, true})
            {
                const int rows = columns / 2;
                std::vector<float> forcing(columns * rows), drag(columns * rows, 43200.0f);
                for (int y = 0; y < rows; ++y)
                    for (int x = 0; x < columns; ++x)
                        forcing[y * columns + x] = static_cast<float>(std::cos(2 * pi * (x + 0.5) / columns) *
                            (1.0 + 0.6 * std::sin(pi * (y + 0.5) / rows)));
                const auto solve = [&](Solver solver)
                {
                    return climateatmosphere::solveSteadyStationaryWavePressure(columns, rows, forcing, drag,
                        48.0f, 190080.0f, 1.225f, 6371000.0f, 7.2921159e-5f, direction,
                        rowProjection, 5000, 1.0e-8f, 60, solver);
                };
                const auto uniform = solve(Solver::Zonal);
                expect(uniform.converged && uniform.iterations <= 2,
                    "longitude-uniform mobility must be inverted through both polar boundaries, rotations and projection modes");
                for (int y = 0; y < rows; ++y)
                    for (int x = columns / 4; x < columns / 2; ++x) drag[y * columns + x] = 8640.0f;
                const auto reference = solve(Solver::ProjectedJacobi), spatial = solve(Solver::Zonal);
                double difference = 0.0;
                for (std::size_t i = 0; i < forcing.size(); ++i)
                    difference = std::max(difference, std::abs(double(reference.pressureAnomalyHpa[i]) - spatial.pressureAnomalyHpa[i]));
                expect(reference.converged && spatial.converged && difference < 2.0e-6,
                    "spatial preconditioning must solve the same constrained equation as strict projected Jacobi with variable drag");
                expect(spatial.returnedRelativeResidual < 3.0e-6,
                    "float pressure returned after projection must still satisfy the original operator within storage precision");
                if (rowProjection)
                {
                    const auto legacy = solve(Solver::LegacyJacobi);
                    expect(legacy.returnedRelativeResidual > 100 * spatial.returnedRelativeResidual,
                        "variable-drag regression must expose the historical post-projection residual corruption");
                }
            }
}

int main(int argc, char** argv)
{
    testSurfaceMomentumAdvection();
    testSurface850Exchange();
    {
        constexpr int columns = 32, rows = 16;
        std::vector<float> forcing(columns * rows), zero(columns * rows, 0.0f), zonal(rows, 0.0f);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
                forcing[y * columns + x] = static_cast<float>(std::cos(2.0 * 3.141592653589793 * 2 * x / columns) *
                    std::cos(3.141592653589793 * 3 * (y + 0.5) / rows));
        climateatmosphere::ModeSeparatedCirculationConfig config;
        config.enabled = {false, true, true, false};
        config.maximumZonalWavenumber = 7;
        config.maximumMeridionalWavenumber = 7;
        const auto result = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonal, zero, forcing, config);
        double error = 0.0;
        for (int cell = 0; cell < columns * rows; ++cell)
            error = std::max(error, std::abs(static_cast<double>(result.surfaceStationarySolver.equilibriumPressureAnomalyHpa[cell]) - 0.46875 * forcing[cell]));
        expect(error < 1.0e-4, "stationary forcing must preserve modal phase and apply its documented spectral gain");
        if (true)
        {
            for (int y = 0; y < rows; ++y)
                for (int x = 0; x < columns; ++x) forcing[y * columns + x] = x < columns / 4 ? 1.0f : 0.0f;
            const auto step = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonal, zero, forcing, config);
            const auto& values = step.surfaceStationarySolver.equilibriumPressureAnomalyHpa;
            expect(*std::min_element(values.begin(), values.end()) >= -0.25001f &&
                *std::max_element(values.begin(), values.end()) <= 0.75001f,
                "Fejer filtering must not overshoot a sharp forcing step after removal of its zonal mean");
        }
    }

    if ((argc == 3 || argc == 4) && std::string(argv[1]) == "--replay")
    {
        auto f = climateatmosphere::detail::StationaryReplay::read(argv[2]);
        if (argc == 4) f.tolerance = std::stof(argv[3]);
        std::vector<climateatmosphere::StationaryWaveResponse> results;
        std::cout << std::setprecision(12) << "{\"fixture\":" << std::quoted(argv[2])
            << ",\"tolerance\":" << f.tolerance << ",\"cases\":[";
        for (int method = 0; method < 3; ++method)
        {
            const auto started = std::chrono::steady_clock::now();
            results.push_back(climateatmosphere::solveSteadyStationaryWavePressure(f.columns, f.rows, f.forcing, f.drag,
                f.depth, f.damping, f.density, f.radius, f.rotation, f.direction, f.rowProjection, f.limit,
                f.tolerance, f.restart, static_cast<climateatmosphere::StationarySolver>(method)));
            const auto& r = results.back();
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            if (method) std::cout << ',';
            std::cout << "{\"solver\":" << method << ",\"converged\":" << r.converged << ",\"seconds\":" << seconds
                << ",\"iterations\":" << r.iterations << ",\"residual\":" << r.relativeResidual
                << ",\"returned_residual\":" << r.returnedRelativeResidual << '}';
        }
        double difference = 0.0;
        for (std::size_t i = 0; i < f.forcing.size(); ++i)
            difference = std::max(difference, std::abs(double(results[1].pressureAnomalyHpa[i]) - results[2].pressureAnomalyHpa[i]));
        std::cout << "],\"projected_pressure_difference_hpa\":" << difference << "}\n";
        return results[2].converged ? 0 : 1;
    }
    if (argc == 3 && std::string(argv[1]) == "--benchmark")
    {
        const int columns = std::stoi(argv[2]), rows = columns / 2;
        if (!climategrid::validGlobalGridDimensions(columns, rows) || columns > 512) return 2;
        std::vector<float> forcing(columns * rows), drag(columns * rows);
        constexpr double pi = 3.14159265358979323846;
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const double lon = 2.0 * pi * (x + 0.5) / columns, lat = pi * (0.5 - (y + 0.5) / rows);
                forcing[y * columns + x] = static_cast<float>(12.0 * std::cos(lat) * std::cos(lat) *
                    (std::cos(3 * lon) + 0.3 * std::sin(7 * lon)) * std::cos(4 * lat));
                drag[y * columns + x] = (x > columns / 5 && x < columns / 2 && y > rows / 5 && y < 4 * rows / 5) ? 8640.0f : 43200.0f;
            }
        std::cout << std::setprecision(12) << "{\"columns\":" << columns << ",\"cases\":[";
        for (int method = 0; method < 3; ++method)
        {
            const auto start = std::chrono::steady_clock::now();
            const auto r = climateatmosphere::solveSteadyStationaryWavePressure(columns, rows, forcing, drag,
                48.0f, 190080.0f, 1.225f, 6371000.0f, 7.2921159e-5f, 1.0f, true, 2000, 1.0e-4f, 60,
                static_cast<climateatmosphere::StationarySolver>(method));
            const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            if (method) std::cout << ',';
            std::cout << "{\"solver\":" << method << ",\"converged\":" << r.converged
                << ",\"seconds\":" << seconds << ",\"iterations\":" << r.iterations << ",\"residual\":" << r.relativeResidual
                << ",\"returned_residual\":" << r.returnedRelativeResidual << ",\"operator_seconds\":" << r.operatorSeconds
                << ",\"preconditioner_seconds\":" << r.preconditionerSeconds << ",\"orthogonalization_seconds\":" << r.orthogonalizationSeconds << '}';
            if (method == 2 && !r.converged) return 1;
        }
        std::cout << "]}\n";
        return 0;
    }
    if (argc != 1) return 2;
    testStationaryPreconditioning();
    testSubpolarDeepening();
    testUpperThermalProjection();
    testReceiverSurfaceDrag();
    testKatabaticDrainage();
    testFrozenSurfaceDrag();
    testConstantForceAcrossCoast();
    testParallelSurfaceDrag();
    testRegionalThermalTrough();
    testReceiverAscentCorrection();
    constexpr float pi = 3.14159265358979323846f;
    constexpr float earthRotation = 7.2921159e-5f;
    expect(
        std::abs(climateatmosphere::coriolisParameterPerSecond(0.0f, earthRotation)) < 1.0e-9f,
        "Coriolis acceleration must vanish at the equator");
    expect(
        climateatmosphere::coriolisParameterPerSecond(45.0f, earthRotation) > 0.0f &&
            climateatmosphere::coriolisParameterPerSecond(-45.0f, earthRotation) < 0.0f,
        "Coriolis sign must reverse across the equator");

    const float heightResponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        100000.0f, 50000.0f);
    expect(
        std::abs(heightResponse - 20.27f) < 0.05f,
        "1000-to-500 hPa thickness must change by about 20.27 metres per kelvin");

    const float earthHadleyEdge = climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
        60.0f, 10000.0f, 288.0f, 9.80665f, earthRotation, 6371000.0f);
    expect(
        earthHadleyEdge > 20.0f && earthHadleyEdge < 26.0f,
        "Held-Hou scaling must place the Earth-like Hadley edge in the subtropics");
    expect(
        climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
            60.0f, 10000.0f, 288.0f, 9.80665f, 0.0f, 6371000.0f) == 90.0f,
        "a non-rotating atmosphere must permit pole-to-pole overturning");

    const float warmSurfacePressure = climateatmosphere::thermalSurfacePressureAnomalyHpa(
        10.0f, 1000.0f, 288.0f, 0.12f);
    const float coldSurfacePressure = climateatmosphere::thermalSurfacePressureAnomalyHpa(
        -10.0f, 1000.0f, 288.0f, 0.12f);
    expect(
        warmSurfacePressure < 0.0f && coldSurfacePressure > 0.0f &&
            std::abs(warmSurfacePressure + coldSurfacePressure) < 0.0001f,
        "warm columns must form thermal lows and cold columns thermal highs");

    const float warmModePressure = climateatmosphere::thermalModePressureAnomalyHpa(
        10.0f, 1.225f, 100000.0f, 50000.0f);
    expect(
        std::abs(warmModePressure + 24.37f) < 0.05f,
        "the thermal-mode pressure must equal the hypsometric geopotential anomaly");

    const float equatorialPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(0.0f, 0.0f, 25.0f, 10.0f);
    const float subtropicalPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(25.0f, 0.0f, 25.0f, 10.0f);
    const float subpolarPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(57.5f, 0.0f, 25.0f, 10.0f);
    const float polarPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(90.0f, 0.0f, 25.0f, 10.0f);
    expect(
        equatorialPressure < 0.0f && subtropicalPressure > 0.0f &&
            subpolarPressure < 0.0f && polarPressure > 0.0f,
        "axisymmetric overturning must produce equatorial and subpolar lows with subtropical and polar highs");

    constexpr int responseColumns = 72;
    constexpr int responseRows = 37;
    const auto responseIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * responseColumns + x);
    };
    std::vector<float> thermalImpulse(
        static_cast<size_t>(responseColumns * responseRows),
        0.0f);
    const int impulseColumn = responseColumns / 2;
    const int equatorRow = responseRows / 2;
    thermalImpulse[responseIndex(impulseColumn, equatorRow)] = 1.0f;
    const auto progradeResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        1.0f);
    const auto retrogradeResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        -1.0f);
    expect(
        progradeResponse[responseIndex(impulseColumn - 1, equatorRow)] >
            progradeResponse[responseIndex(impulseColumn + 1, equatorRow)] &&
        retrogradeResponse[responseIndex(impulseColumn + 1, equatorRow)] >
            retrogradeResponse[responseIndex(impulseColumn - 1, equatorRow)],
        "the tropical response must extend westward relative to planetary rotation");
    const int extratropicalRow = 10;
    thermalImpulse.assign(static_cast<size_t>(responseColumns * responseRows), 0.0f);
    thermalImpulse[responseIndex(impulseColumn, extratropicalRow)] = 1.0f;
    const auto extratropicalResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        1.0f);
    expect(
        std::abs(
            extratropicalResponse[responseIndex(impulseColumn - 1, extratropicalRow)] -
            extratropicalResponse[responseIndex(impulseColumn + 1, extratropicalRow)]) < 1.0e-6f,
        "the extratropical thermal response must be zonally symmetric");
    double equatorialResponseSum = 0.0;
    for (int x = 0; x < responseColumns; x++)
        equatorialResponseSum += progradeResponse[responseIndex(x, equatorRow)];
    expect(
        equatorialResponseSum > 0.0 && equatorialResponseSum < 1.0,
        "the nonlocal response must redistribute an impulse without creating a new extremum");

    constexpr int topographicColumns = 24;
    constexpr int topographicRows = topographicColumns / 2;
    const auto topographicIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * topographicColumns + x);
    };
    std::vector<float> ridgeTerrain(
        static_cast<size_t>(topographicColumns * topographicRows),
        0.0f);
    std::vector<float> ridgeEastWind(ridgeTerrain.size(), 10.0f);
    std::vector<float> ridgeSouthWind(ridgeTerrain.size(), 0.0f);
    for (int y = 0; y < topographicRows; y++)
        ridgeTerrain[topographicIndex(10, y)] = 2000.0f;
    const auto ridgeForcing = climateatmosphere::mechanicalTopographicPressureForcingHpa(
        topographicColumns,
        topographicRows,
        ridgeTerrain,
        ridgeEastWind,
        ridgeSouthWind,
        3.0f,
        2000.0f,
        2.0f,
        1.0f,
        8.0f,
        10.0f,
        30.0f);
    expect(
        ridgeForcing[topographicIndex(7, 2)] > 0.0f &&
            ridgeForcing[topographicIndex(13, 2)] < 0.0f,
        "westerly flow across a ridge must produce a windward high and lee trough");
    ridgeEastWind.assign(ridgeTerrain.size(), -10.0f);
    const auto reversedRidgeForcing =
        climateatmosphere::mechanicalTopographicPressureForcingHpa(
            topographicColumns,
            topographicRows,
            ridgeTerrain,
            ridgeEastWind,
            ridgeSouthWind,
            3.0f,
            2000.0f,
            2.0f,
            1.0f,
            8.0f,
            10.0f,
            30.0f);
    expect(
        reversedRidgeForcing[topographicIndex(13, 2)] > 0.0f &&
            reversedRidgeForcing[topographicIndex(7, 2)] < 0.0f,
        "topographic pressure forcing must reverse when the background wind reverses");

    const auto tropicalRidge = climateatmosphere::mechanicalTopographicPressureForcingHpa(
        topographicColumns, topographicRows, ridgeTerrain, ridgeEastWind, ridgeSouthWind,
        3.0f, 2000.0f, 2.0f, 1.0f, 8.0f, 0.0f, 0.0f);
    expect(tropicalRidge[topographicIndex(13, topographicRows / 2)] > 0.0f &&
        tropicalRidge[topographicIndex(7, topographicRows / 2)] < 0.0f,
        "trade winds must retain windward and lee pressure responses to tropical mountains");
    const auto flatTropics = climateatmosphere::mechanicalTopographicPressureForcingHpa(
        topographicColumns, topographicRows, std::vector<float>(ridgeTerrain.size(), 0.0f),
        ridgeEastWind, ridgeSouthWind, 3.0f, 2000.0f, 2.0f, 1.0f, 8.0f, 0.0f, 0.0f);
    expect(std::all_of(flatTropics.begin(), flatTropics.end(), [](float v) { return v == 0.0f; }),
        "removing the tropical gate must not invent terrain forcing over flat surfaces");

    const auto spacing = climateatmosphere::cellSpacingMetres(0.0f, 2048, 1024, 6371000.0f);
    expect(
        spacing.zonalMetres > 19000.0f && spacing.zonalMetres < 20000.0f &&
            spacing.meridionalMetres > 19000.0f && spacing.meridionalMetres < 20000.0f,
        "Earth benchmark cells must be about 19.5 kilometres at the equator");

    const auto unrotated = climateatmosphere::steadyRayleighCoriolisWind(
        0.001f, 0.002f, 0.0f, 1000.0f, earthRotation);
    expect(
        std::abs(unrotated.eastMetresPerSecond - 1.0f) < 0.001f &&
            std::abs(unrotated.southMetresPerSecond + 2.0f) < 0.001f,
        "Rayleigh flow at the equator must follow the pressure-gradient acceleration");

    const auto northern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, 45.0f, 86400.0f, earthRotation);
    const auto southern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, -45.0f, 86400.0f, earthRotation);
    expect(
        northern.eastMetresPerSecond > 0.0f && southern.eastMetresPerSecond < 0.0f,
        "the same meridional height force must produce opposite zonal flow across the equator");

    const auto quadraticEquatorial = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.001f, 0.0f, 0.0f, 0.001f, 100.0f, earthRotation);
    const auto quadraticEquatorialStronger = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.004f, 0.0f, 0.0f, 0.001f, 100.0f, earthRotation);
    expect(
        std::abs(quadraticEquatorial.eastMetresPerSecond - 10.0f) < 0.001f &&
            std::abs(quadraticEquatorial.southMetresPerSecond) < 0.001f &&
            std::abs(quadraticEquatorialStronger.eastMetresPerSecond - 20.0f) < 0.001f,
        "quadratic drag must keep equatorial flow finite and scale speed with the square root of force");

    const auto quadraticNorthern = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.0f, 0.001f, 45.0f, 0.0013f, 300.0f, earthRotation);
    const auto quadraticSouthern = climateatmosphere::steadyQuadraticDragCoriolisWind(
        0.0f, 0.001f, -45.0f, 0.0013f, 300.0f, earthRotation);
    expect(
        quadraticNorthern.eastMetresPerSecond > 0.0f &&
            quadraticSouthern.eastMetresPerSecond < 0.0f,
        "quadratic surface drag must preserve the Coriolis reversal across the equator");

    constexpr int waveColumns = 64;
    constexpr int waveRows = waveColumns / 2;
    constexpr int waveCentreX = 16;
    constexpr int waveCentreY = waveRows / 2;
    const auto waveIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y) * waveColumns + x;
    };
    std::vector<float> waveForcing(waveColumns * waveRows, 0.0f);
    std::vector<float> waveDragTime(waveForcing.size(), 43200.0f);
    for (int y = 0; y < waveRows; y++)
    {
        for (int x = 0; x < waveColumns; x++)
        {
            int distanceX = x - waveCentreX;
            if (distanceX > waveColumns / 2)
                distanceX -= waveColumns;
            if (distanceX < -waveColumns / 2)
                distanceX += waveColumns;
            const int distanceY = y - waveCentreY;
            waveForcing[waveIndex(x, y)] = std::exp(
                -static_cast<float>(distanceX * distanceX) / 25.0f -
                static_cast<float>(distanceY * distanceY) / 9.0f);
        }
    }
    const auto progradeWave = climateatmosphere::solveSteadyStationaryWavePressure(
        waveColumns,
        waveRows,
        waveForcing,
        waveDragTime,
        48.0f,
        2.2f * 86400.0f,
        1.225f,
        6371000.0f,
        earthRotation,
        1.0f,
        true,
        200,
        1.0e-4f);
    const auto retrogradeWave = climateatmosphere::solveSteadyStationaryWavePressure(
        waveColumns,
        waveRows,
        waveForcing,
        waveDragTime,
        48.0f,
        2.2f * 86400.0f,
        1.225f,
        6371000.0f,
        earthRotation,
        -1.0f,
        true,
        200,
        1.0e-4f);
    const auto fullPressureWave = climateatmosphere::solveSteadyStationaryWavePressure(
        waveColumns,
        waveRows,
        waveForcing,
        waveDragTime,
        48.0f,
        2.2f * 86400.0f,
        1.225f,
        6371000.0f,
        earthRotation,
        1.0f,
        false,
        300,
        1.0e-4f);
    double waveMagnitude = 0.0;
    double waveAsymmetry = 0.0;
    double reversedWaveDifference = 0.0;
    double maximumRowMean = 0.0;
    double fullPressureAreaTotal = 0.0;
    double fullPressureAreaWeight = 0.0;
    double fullPressureRowMeanMagnitude = 0.0;
    for (int y = 0; y < waveRows; y++)
    {
        double rowMean = 0.0;
        double fullPressureRowMean = 0.0;
        for (int x = 0; x < waveColumns; x++)
        {
            const int mirroredX = (2 * waveCentreX - x + waveColumns) % waveColumns;
            const float value = progradeWave.pressureAnomalyHpa[waveIndex(x, y)];
            waveMagnitude += std::abs(value);
            waveAsymmetry += std::abs(
                value - progradeWave.pressureAnomalyHpa[waveIndex(mirroredX, y)]);
            reversedWaveDifference += std::abs(
                value - retrogradeWave.pressureAnomalyHpa[waveIndex(mirroredX, y)]);
            rowMean += value;
            fullPressureRowMean +=
                fullPressureWave.pressureAnomalyHpa[waveIndex(x, y)];
        }
        maximumRowMean = std::max(
            maximumRowMean,
            std::abs(rowMean / static_cast<double>(waveColumns)));
        fullPressureRowMean /= static_cast<double>(waveColumns);
        fullPressureRowMeanMagnitude += std::abs(fullPressureRowMean);
        const double latitudeRadians =
            (90.0 - 180.0 * (static_cast<double>(y) + 0.5) /
                    static_cast<double>(waveRows)) *
            3.14159265358979323846 / 180.0;
        const double areaWeight = std::max(0.0, std::cos(latitudeRadians));
        fullPressureAreaTotal += fullPressureRowMean * areaWeight;
        fullPressureAreaWeight += areaWeight;
    }
    if (!progradeWave.converged || !retrogradeWave.converged)
    {
        std::cerr
            << "stationary-wave diagnostics prograde_iterations=" << progradeWave.iterations
            << " prograde_residual=" << progradeWave.relativeResidual
            << " retrograde_iterations=" << retrogradeWave.iterations
            << " retrograde_residual=" << retrogradeWave.relativeResidual << '\n';
    }
    expect(
        progradeWave.converged && retrogradeWave.converged &&
            progradeWave.relativeResidual < 1.0e-4f &&
            retrogradeWave.relativeResidual < 1.0e-4f,
        "stationary-wave pressure solves must converge to their requested residual");
    expect(
        waveMagnitude > 0.0 && waveAsymmetry / waveMagnitude > 0.01 &&
            reversedWaveDifference / waveMagnitude < 0.001,
        "rotation must create a longitudinally asymmetric response that mirrors when rotation reverses");
    expect(
        maximumRowMean < 1.0e-5,
        "the stationary-wave response must preserve each row's zonal mean");
    expect(
        fullPressureWave.converged && fullPressureRowMeanMagnitude > 0.01 &&
            std::abs(fullPressureAreaTotal / fullPressureAreaWeight) < 1.0e-5,
        "the full pressure solve must retain zonal structure while conserving global mean pressure");
    expect(!progradeWave.residualHistory.empty() && progradeWave.restartCycles > 0 &&
            progradeWave.residualHistory.back() <= progradeWave.residualHistory.front(),
        "stationary solves must retain restart-cycle physical residual histories");

    constexpr int modeColumns = 16;
    constexpr int modeRows = 8;
    constexpr std::size_t modeCellCount = modeColumns * modeRows;
    std::vector<float> absorbed(modeCellCount, 220.0f);
    std::vector<float> outgoing(modeCellCount, 220.0f);
    std::vector<float> sensible(modeCellCount, 0.0f);
    std::vector<float> condensation(modeCellCount, 0.0f);
    for (int y = 0; y < modeRows; y++)
    {
        for (int x = 0; x < modeColumns; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * modeColumns + x;
            absorbed[cell] += 20.0f * std::cos(2.0f * pi * x / modeColumns);
            sensible[cell] = 4.0f * std::sin(2.0f * pi * x / modeColumns);
            condensation[cell] = x < modeColumns / 2 ? 10.0f : 0.0f;
        }
    }
    const auto heating = climateatmosphere::diagnoseDiabaticHeating(
        modeColumns,
        modeRows,
        absorbed,
        outgoing,
        sensible,
        condensation,
        86400.0f,
        0.35f,
        1.0f);
    expect(heating.areaWeightedLatentHeatingWm2 > 140.0 &&
            heating.areaWeightedLatentHeatingWm2 < 150.0,
        "latent heating must convert millimetres per day to watts per square metre");
    expect(heating.maximumAbsoluteRowMeanWm2 < 1.0e-5f,
        "stationary heating projection must have a controlled zero zonal mean");
    const auto noLatentHeating = climateatmosphere::diagnoseDiabaticHeating(
        modeColumns,
        modeRows,
        absorbed,
        outgoing,
        sensible,
        condensation,
        86400.0f,
        0.35f,
        0.0f);
    expect(std::abs(noLatentHeating.areaWeightedLatentHeatingWm2) < 1.0e-9,
        "the latent projection switch must prevent double-counting hydrology energy");

    const auto earthParameters = climateatmosphere::diagnoseStationaryParameters(
        0.01f, 10000.0f, 9.80665f, 6371000.0f, earthRotation,
        64, 32, 0.10f, 1200000.0f);
    const auto slowRotationParameters = climateatmosphere::diagnoseStationaryParameters(
        0.01f, 10000.0f, 9.80665f, 6371000.0f, earthRotation * 0.5f,
        64, 32, 0.10f, 1200000.0f);
    expect(earthParameters.equivalentDepthMetres > 0.0f &&
            earthParameters.maximumZonalWavenumber <= modeColumns * 2 &&
            slowRotationParameters.adjustmentLengthMetres >
                earthParameters.adjustmentLengthMetres,
        "stationary parameters must derive from stratification, rotation, and resolution");

    std::vector<float> zonalPressure(modeRows, 0.0f);
    std::vector<float> orographic(modeCellCount, 0.0f);
    for (int y = 0; y < modeRows; y++)
        zonalPressure[y] = 3.0f * std::cos(pi * (y + 0.5f) / modeRows);
    climateatmosphere::ModeSeparatedCirculationConfig modeConfig;
    modeConfig.maximumIterations = 500;
    const auto separated = climateatmosphere::solveModeSeparatedCirculation(
        modeColumns,
        modeRows,
        zonalPressure,
        heating.stationaryProjectedHeatingWm2,
        orographic,
        modeConfig);
    double surfaceMagnitude = 0.0;
    double upperMagnitude = 0.0;
    double maximumZonalTransfer = 0.0;
    for (int y = 0; y < modeRows; y++)
    {
        double rowMean = 0.0;
        for (int x = 0; x < modeColumns; x++)
        {
            const std::size_t cell = static_cast<std::size_t>(y) * modeColumns + x;
            rowMean += separated.surfacePressureAnomalyHpa[cell];
            surfaceMagnitude += std::abs(separated.surfaceEastWindMps[cell]) +
                std::abs(separated.surfaceSouthWindMps[cell]);
            upperMagnitude += std::abs(separated.upperEastWindMps[cell]) +
                std::abs(separated.upperSouthWindMps[cell]);
        }
        maximumZonalTransfer = std::max(
            maximumZonalTransfer,
            std::abs(rowMean / modeColumns - zonalPressure[y]));
    }
    expect(separated.surfaceStationarySolver.converged &&
            separated.upperStationarySolver.converged &&
            maximumZonalTransfer < 1.0e-4,
        "surface and upper stationary modes must converge without collapsing the zonal mode");
    expect(surfaceMagnitude > 0.0 && upperMagnitude > 0.0,
        "separately closed surface and upper modes must both respond to heating");
    modeConfig.enabled.stationary = false;
    modeConfig.enabled.upper = false;
    const auto zonalOnly = climateatmosphere::solveModeSeparatedCirculation(
        modeColumns,
        modeRows,
        zonalPressure,
        heating.stationaryProjectedHeatingWm2,
        orographic,
        modeConfig);
    expect(std::all_of(
            zonalOnly.upperHeightAnomalyMetres.begin(),
            zonalOnly.upperHeightAnomalyMetres.end(),
            [](float value) { return value == 0.0f; }),
        "upper-only and stationary-only responses must remain isolatable");

    const std::vector<float> zeroHeating(modeCellCount, 0.0f);
    const std::vector<float> zeroZonal(modeRows, 0.0f);
    modeConfig = {};
    const auto zeroResponse = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zeroZonal, zeroHeating, zeroHeating, modeConfig);
    expect(zeroResponse.surfaceStationarySolver.converged && zeroResponse.upperStationarySolver.converged &&
        zeroResponse.areaWeightedKineticEnergyJm2 == 0.0, "zero forcing must produce zero flow without noise");
    std::vector<float> thermalPressure(modeCellCount);
    for (int y = 0; y < modeRows; ++y)
        for (int x = 0; x < modeColumns; ++x)
            thermalPressure[y * modeColumns + x] = climateatmosphere::thermalModePressureAnomalyHpa(
                5.0f * std::cos(2.0f * pi * x / modeColumns), 1.225f, 100000.0f, 70000.0f);
    const auto thermalResponse = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zeroZonal, zeroHeating, thermalPressure, modeConfig);
    double thermalAlignment = 0.0;
    double thermalRowMean = 0.0;
    for (int y = 0; y < modeRows; ++y)
    {
        double rowMean = 0.0;
        for (int x = 0; x < modeColumns; ++x)
        {
            const int cell = y * modeColumns + x;
            rowMean += thermalResponse.surfacePressureAnomalyHpa[cell];
            thermalAlignment += thermalResponse.surfacePressureAnomalyHpa[cell] * thermalPressure[cell];
        }
        thermalRowMean = std::max(thermalRowMean, std::abs(rowMean / modeColumns));
    }
    expect(thermalResponse.surfaceStationarySolver.converged && thermalAlignment > 0.0 &&
        thermalResponse.areaWeightedKineticEnergyJm2 > 0.0 && thermalRowMean < 1.0e-5,
        "surface thermal gradients must drive pressure and wind without diabatic heating or zonal mass transfer");
    // The two numerical models should approach the same forced, damped mode.
    // Feeding the already adjusted steady pressure into the evolving model
    // instead of its source attenuates this mode a second time.
    constexpr int comparisonColumns = 32, comparisonRows = 16;
    std::vector<float> waveSource(comparisonColumns * comparisonRows);
    for (int y = 0; y < comparisonRows; ++y)
        for (int x = 0; x < comparisonColumns; ++x)
            waveSource[y * comparisonColumns + x] = 4.0f *
                std::sin(pi * (y + 0.5f) / comparisonRows) * std::cos(2.0f * pi * x / comparisonColumns);
    const auto steadyMode = climateatmosphere::solveSteadyStationaryWavePressure(
        comparisonColumns, comparisonRows, waveSource, std::vector<float>(waveSource.size(), 43200.0f),
        48.0f, 86400.0f, 1.225f, 6371000.0f, 0.0f, 1.0f, true, 1000, 1.0e-5f);
    climateweather::ShallowWaterConfig evolvingConfig;
    evolvingConfig.rotationRatePerSecond = 0.0f;
    evolvingConfig.lowerDragTimeSeconds = 43200.0f;
    evolvingConfig.heightRelaxationTimeSeconds = 86400.0f;
    evolvingConfig.baroclinicCoupling = 0.0f;
    const float pressurePerMetre = 1.225f * evolvingConfig.gravityMetresPerSecondSquared / 100.0f;
    evolvingConfig.lowerMeanDepthMetres = 48.0f / pressurePerMetre;
    auto evolvingMode = climateweather::makeState(comparisonColumns, comparisonRows, 1, 73);
    climateweather::ShallowWaterForcing evolvingForcing;
    evolvingForcing.equilibriumHeightMetres = {steadyMode.equilibriumPressureAnomalyHpa};
    for (float& value : evolvingForcing.equilibriumHeightMetres[0]) value /= pressurePerMetre;
    bool evolvingStable = true;
    for (int step = 0; step < 64; ++step)
    {
        const auto d = climateweather::advance(evolvingMode, evolvingConfig, evolvingForcing, 21600.0f);
        evolvingStable = evolvingStable && d.finite && d.bounded;
    }
    double modeDifference = 0.0, modeScale = 0.0;
    for (int y = 0; y < comparisonRows; ++y)
        for (int x = 0; x < comparisonColumns; ++x)
        {
            const int cell = y * comparisonColumns + x;
            const double weight = std::sin(pi * (y + 0.5) / comparisonRows);
            const double expected = steadyMode.pressureAnomalyHpa[cell];
            modeDifference += weight * std::pow(pressurePerMetre *
                evolvingMode.layers[0].heightAnomalyMetres[cell] - expected, 2);
            modeScale += weight * expected * expected;
        }
    expect(steadyMode.converged && evolvingStable && std::sqrt(modeDifference / modeScale) < 0.25,
        "steady and evolving pressure modes must share unadjusted forcing without double attenuation");
    for (int mask = 0; mask < 16; ++mask)
    {
        modeConfig.enabled = {(mask & 1) != 0, (mask & 2) != 0, (mask & 4) != 0, (mask & 8) != 0};
        modeConfig.zonalUpperHeightMetres = zonalPressure;
        const auto isolated = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
            zonalPressure, heating.stationaryProjectedHeatingWm2, orographic, modeConfig);
        const auto zero = [](const auto& values) { return std::all_of(values.begin(), values.end(), [](float v) { return v == 0.0f; }); };
        expect((modeConfig.enabled.surface || (zero(isolated.surfaceEastWindMps) && zero(isolated.surfacePressureAnomalyHpa))) &&
            (modeConfig.enabled.upper || (zero(isolated.upperEastWindMps) && zero(isolated.upperHeightAnomalyMetres))),
            "all 16 mode-isolation combinations must respect disabled layers");
    }
    modeConfig = {};
    modeConfig.maximumIterations = 1;
    modeConfig.stationarySolver = climateatmosphere::StationarySolver::LegacyJacobi;
    modeConfig.relativeTolerance = 1.0e-12f;
    const auto failed = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, heating.stationaryProjectedHeatingWm2, orographic, modeConfig);
    expect(!failed.surfaceStationarySolver.converged && failed.surfacePressureAnomalyHpa[0] == zonalPressure[0],
        "a failed stationary solve must fall back to the independently closed zonal mode");
    std::vector<float> mountain(modeCellCount), jets(modeCellCount, 20.0f), calm(modeCellCount, 0.0f);
    for (int y = 0; y < modeRows; ++y)
        for (int x = 0; x < modeColumns; ++x)
            mountain[y * modeColumns + x] = 200.0f * std::cos(2.0f * pi * x / modeColumns);
    modeConfig = {};
    const auto lowerMountain = climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
        mountain, jets, 0.01f, 3000.0f, 86400.0f, modeConfig);
    const auto upperMountain = climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
        mountain, jets, 0.01f, 8000.0f, 86400.0f, modeConfig);
    const auto calmMountain = climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
        mountain, calm, 0.01f, 5000.0f, 86400.0f, modeConfig);
    expect(lowerMountain != upperMountain && std::all_of(upperMountain.begin(), upperMountain.end(),
        [](float v) { return std::isfinite(v); }) && std::all_of(calmMountain.begin(), calmMountain.end(),
        [](float v) { return v == 0.0f; }), "mountain waves must propagate/damp vertically and vanish without incident wind");
    {
        const auto effective = climateatmosphere::detail::effectiveMountainWaveHeightMetres;
        expect(effective(100.0f, 20.0f, 0.01f) == 100.0f &&
               effective(100.0f, -20.0f, 0.01f) == 100.0f,
            "weak mountain perturbations must retain their height for either wind direction");
        expect(effective(5000.0f, 10.0f, 0.01f) == effective(10000.0f, 10.0f, 0.01f) &&
               effective(5000.0f, 10.0f, 0.02f) < effective(5000.0f, 10.0f, 0.01f),
            "blocked launch height must saturate and decrease as stability increases");
        expect(effective(1000.0f, 0.0f, 0.01f) == 0.0f &&
               effective(-100.0f, 20.0f, 0.01f) == 0.0f &&
               effective(1000.0f, 20.0f, 0.0f) == 0.0f,
            "calm, submerged and unstratified inputs must not launch a limited mountain wave");
        std::vector<float> bounded(modeCellCount), twiceTall(modeCellCount);
        for (int cell = 0; cell < modeCellCount; ++cell)
        {
            const float peak = cell % modeColumns < modeColumns / 4 ? 5000.0f : 0.0f;
            bounded[cell] = effective(peak, jets[cell], 0.01f);
            twiceTall[cell] = effective(2.0f * peak, jets[cell], 0.01f);
        }
        expect(climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
            bounded, jets, 0.01f, 5000.0f, 86400.0f, modeConfig) ==
            climateatmosphere::upperOrographicHeightForcing(modeColumns, modeRows,
            twiceTall, jets, 0.01f, 5000.0f, 86400.0f, modeConfig),
            "a taller blocked obstacle must not indefinitely amplify the linear upper response");
    }
    expect(climateatmosphere::diagnoseBruntVaisalaFrequency(280.0f, 0.006f, 9.80665f) > 0.0f &&
        climateatmosphere::diagnoseBruntVaisalaFrequency(280.0f, 0.012f, 9.80665f) == 0.0f,
        "stratification must distinguish stable and convectively unstable lapse rates");
    const auto refined = climateatmosphere::diagnoseStationaryParameters(0.01f, 10000.0f, 9.80665f,
        6371000.0f, earthRotation, 128, 64, 0.1f, 1200000.0f);
    const auto twiceRefined = climateatmosphere::diagnoseStationaryParameters(0.01f, 10000.0f, 9.80665f,
        6371000.0f, earthRotation, 256, 128, 0.1f, 1200000.0f);
    expect(refined.maximumZonalWavenumber == twiceRefined.maximumZonalWavenumber,
        "refinement must not invent physical forcing bandwidth");
    const auto largeRestart = climateatmosphere::solveSteadyStationaryWavePressure(
        modeColumns, modeRows, heating.stationaryProjectedHeatingWm2, std::vector<float>(modeCellCount, 43200.0f),
        48.0f, 190080.0f, 1.225f, 6371000.0f, earthRotation, 1.0f, true, 500, 1.0e-4f, 128);
    expect(largeRestart.converged && largeRestart.relativeResidual <= 1.0e-4f,
        "configurable larger GMRES restart windows must retain physical residual acceptance");
    modeConfig = {};
    modeConfig.enabled = {true, false, true, false};
    climateatmosphere::ColumnHeatingInput column;
    column.incomingSolarWm2 = 340.0;
    column.sensibleHeatingWm2 = 20.0;
    column.condensationMm = {1.0, 4.0};
    column.reevaporationMm = 2.0;
    column.surfaceEvaporationMm = 3.0;
    const auto heatColumn = climateatmosphere::diagnoseColumnHeating(column);
    expect(std::abs(heatColumn.closureResidualWm2) < 1.0e-10 &&
        heatColumn.latentWm2[0] < 0.0 && heatColumn.latentWm2[1] > 0.0,
        "grey radiation and phase changes must close energy and locate re-evaporative cooling below condensation");
    column.longwaveOpticalDepth = {0.0, 0.0};
    column.shortwaveOpticalDepth = {0.0, 0.0};
    const auto transparent = climateatmosphere::diagnoseColumnHeating(column);
    expect(transparent.radiativeWm2[0] == 0.0 && transparent.radiativeWm2[1] == 0.0,
        "transparent air must not receive surface radiative heating");
    const auto uniformDrag = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, zeroHeating, zeroHeating, modeConfig);
    modeConfig.surfaceDragCoefficients.assign(modeCellCount, modeConfig.surfaceDragCoefficient);
    for (int y = 0; y < modeRows; ++y)
        for (int x = modeColumns / 2; x < modeColumns; ++x)
            modeConfig.surfaceDragCoefficients[y * modeColumns + x] *= 8.0f;
    const auto roughLand = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, zeroHeating, zeroHeating, modeConfig);
    const int equatorialRow = modeRows / 2;
    const int oceanCell = equatorialRow * modeColumns;
    const int landCell = oceanCell + modeColumns / 2;
    const auto speed = [](const auto& flow, int cell) {
        return std::hypot(flow.surfaceEastWindMps[cell], flow.surfaceSouthWindMps[cell]); };
    expect(speed(roughLand, landCell) < speed(roughLand, oceanCell) &&
        speed(roughLand, oceanCell) == speed(uniformDrag, oceanCell),
        "the final quadratic wind solver must use local drag without changing ocean cells");
    modeConfig.adjustZonalSurfacePressure = true;
    modeConfig.surfaceDragTimesSeconds.assign(modeCellCount, modeConfig.surfaceDragTimeSeconds);
    for (int y = 0; y < modeRows; ++y)
        for (int x = modeColumns / 2; x < modeColumns; ++x)
            modeConfig.surfaceDragTimesSeconds[y * modeColumns + x] /= 8.0f;
    const auto adjusted = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, zeroHeating, zeroHeating, modeConfig);
    expect(adjusted.surfaceStationarySolver.converged &&
        adjusted.surfaceStationarySolver.returnedRelativeResidual <= modeConfig.relativeTolerance,
        "zonal pressure participates in the accepted full pressure/divergence solve");
    expect(adjusted.surfacePressureAnomalyHpa == adjusted.surfaceStationarySolver.pressureAnomalyHpa,
        "adjusted zonal pressure must not be added a second time after the solve");
    double regionalPressure = 0.0;
    for (int y = 0; y < modeRows; ++y)
        regionalPressure = std::max(regionalPressure, static_cast<double>(std::abs(
            adjusted.surfacePressureAnomalyHpa[y * modeColumns] -
            adjusted.surfacePressureAnomalyHpa[y * modeColumns + modeColumns / 2])));
    expect(regionalPressure > 1.0e-4,
        "a rough continent must modify pressure adjustment even with zonally uniform forcing");
    modeConfig = {};
    modeConfig.enabled = {true, false, false, true};
    modeConfig.zonalUpperHeightMetres.resize(modeRows);
    for (int y = 0; y < modeRows; ++y)
        modeConfig.zonalUpperHeightMetres[y] = 40.0f * zonalPressure[y];
    modeConfig.adjustZonalUpperHeight = true;
    const auto adjustedUpper = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, zeroHeating, zeroHeating, modeConfig);
    expect(adjustedUpper.upperStationarySolver.converged &&
        adjustedUpper.upperStationarySolver.returnedRelativeResidual <= modeConfig.relativeTolerance,
        "the full upper height mode must satisfy its pressure/divergence equation");
    double upperReconstructionError = 0.0;
    const float adjustedPressurePerMetre = modeConfig.airDensityKgM3 * modeConfig.gravityMetresPerSecondSquared / 100.0f;
    for (int cell = 0; cell < modeCellCount; ++cell)
        upperReconstructionError = std::max(upperReconstructionError, static_cast<double>(std::abs(
            adjustedUpper.upperHeightAnomalyMetres[cell] -
            adjustedUpper.upperStationarySolver.pressureAnomalyHpa[cell] / adjustedPressurePerMetre)));
    expect(upperReconstructionError < 1.0e-5,
        "adjusted upper forcing must not be added again as a diagnostic height");
    modeConfig = {};
    modeConfig.interlayerMomentumCoupling = 0.0f;
    modeConfig.upperMaximumZonalWavenumber = 0;
    const auto upperFiltered = climateatmosphere::solveModeSeparatedCirculation(modeColumns, modeRows,
        zonalPressure, heating.stationaryProjectedHeatingWm2, orographic, modeConfig);
    expect(upperFiltered.surfacePressureAnomalyHpa == separated.surfacePressureAnomalyHpa &&
        std::all_of(upperFiltered.upperHeightAnomalyMetres.begin(), upperFiltered.upperHeightAnomalyMetres.end(),
            [](float v) { return std::abs(v) < 1.0e-4f; }),
        "upper bandwidth must filter upper forcing independently of surface pressure");
    return failures == 0 ? 0 : 1;
}
