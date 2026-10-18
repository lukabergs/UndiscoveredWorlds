#pragma once

// The climate module's world: terrain and planet inputs (read only) plus the
// state the climate processes own: seasonal and legacy climate fields, sea ice,
// biomes, subtropical highs, the weather state and the solver context. planet
// holds one climateworld and forwards its climate accessors to it; a standalone
// caller owns one and binds its own terrain. Accessor names follow planet so the
// climate code reads the same (docs/architecture/world-model.md).

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

#include "simulations/climate/climate_context.hpp"
#include "core/climate_seasons.hpp"
#include "simulations/climate/climate_weather.hpp"
#include "core/world/field_access.hpp"
#include "core/world/transient_state.hpp"
#include "core/world/world_grid.hpp"

#ifndef ARRAYWIDTH
#define ARRAYWIDTH 2048
#define ARRAYHEIGHT 1024
#endif


// Terrain and planet inputs of the climate processes. Grids hold one value per
// cell of the world grid in its id order (on lat-lon row-major, y * (width + 1)
// + x over 0..width, 0..height) and are owned by the caller.
struct climateinputs
{
    int width = 0;                    // largest x index
    int height = 0;                   // largest y index
    long seed = 0;
    bool rotation = true;
    float tilt = 22.5f;
    float eccentricity = 0.0167f;
    int perihelion = 0;
    float gravity = 1.0f;
    float lunar = 1.0f;
    float tempdecrease = 6.5f;
    int northpolaradjust = 0;
    int southpolaradjust = 0;
    int averagetemp = 14;
    float waterpickup = 1.0f;
    int maxelevation = 0;
    int sealevel = 0;
    int landtotal = 0;
    int seatotal = 0;
    const float* nom = nullptr;   // base elevation, metres above the datum
    const float* mountainheights = nullptr;
    const float* extraelev = nullptr;
    const float* craterrims = nullptr;
    const float* cratercentres = nullptr;
    const float* volcano = nullptr;
    const float* oceanridgeheight = nullptr;
    const float* lakesurface = nullptr;
    const float* riftlakesurface = nullptr;
    const short* special = nullptr;
};

// The climate fields were int16 arrays, some of them scaled. They are stored as
// float in physical units; the integer accessors still return the former
// encodings for code that has not moved to the float accessors, the .uww save
// format and the renderers. These are the factors of the scaled encodings.
namespace climatelegacyscale
{
inline constexpr float verticalvelocity = 100.0f;        // hPa/day
inline constexpr float current = 100.0f;                 // m/s, stored as cm/s
inline constexpr float maritimeinfluence = 1000.0f;      // 0..1
inline constexpr float maritimethermalanomaly = 10.0f;   // C
inline constexpr float convergence = 100.0f;             // mm/month
inline constexpr float orographicmotion = 4.0f;          // m of parcel displacement
}

// One level of the atmosphere column (the 10 m surface or a pressure level):
// seasonal fields on the world grid in its id order (row-major on lat-lon), in
// physical units. The current two-mode model fills 10 m, 850 hPa and 500 hPa
// after each wind solve; a multi-level solver would fill more levels.
struct atmospherelevel
{
    std::string name;                 // "10m", "850hPa", "500hPa"
    float pressurehpa = 0.0f;         // nominal pressure; 0 for the 10 m surface level
    std::array<std::vector<float>, CLIMATESEASONCOUNT> eastwind;            // m/s
    std::array<std::vector<float>, CLIMATESEASONCOUNT> northwind;           // m/s
    std::array<std::vector<float>, CLIMATESEASONCOUNT> geopotentialheight;  // m; empty when not modelled
    std::array<std::vector<float>, CLIMATESEASONCOUNT> verticalvelocity;    // hPa/day, ascent positive; empty when not modelled
    std::array<std::vector<std::uint8_t>, CLIMATESEASONCOUNT> available;    // 1 above ground and solved
};

// What the last climate solve recorded on the solvers' own grids: per-season
// moisture and heat processes, ocean state and weather statistics. Diagnostics
// exports, validation and run bundles read it. Transient: a copy of the world
// does not carry it, as the copy did not run the solve.
struct climatediagnostics
{
    std::array<climatehydrology::SeasonalProcessFields, CLIMATESEASONCOUNT> process;
    bool hasprocess = false;
    std::array<climateocean::OceanState, CLIMATESEASONCOUNT> ocean;
    int oceancolumns = 0;
    int oceanrows = 0;
    bool hasocean = false;
    std::array<std::array<climateweather::WeatherStatistics, CLIMATESEASONCOUNT>, 2> weather;
    int weathercolumns = 0;
    int weatherrows = 0;
    bool hasweather = false;
};

class climateworld
{
    friend class planet;
    friend struct planetfieldaccess;
    friend struct climateworldfieldaccess;

public:
    // Binds the terrain and planet inputs; called before climate processes run.
    void bindinputs(const climateinputs& inputs) const
    {
        const int width = inputs_.width, height = inputs_.height;
        inputs_ = inputs;
        inputs_.width = width;      // the grid size belongs to the state (resize)
        inputs_.height = height;
    }
    const climateinputs& inputs() const { return inputs_; }
    // Sizes the seasonal fields for the world grid and resets the solver state:
    // a lat-lon grid of largest indices width, height, or any cell grid (on the
    // cube width() and height() are -1 and only the CellId accessors apply).
    void resize(int width, int height);
    void resize(const mesh::CellGrid& grid);

    climatecontext::Context& climatesimulation() { noteclimatesolverread(); return climatecontext_; }
    const climatecontext::Context& climatesimulation() const { noteclimatesolverread(); return climatecontext_; }

    // Inputs.
    int width() const { return inputs_.width; }
    int height() const { return inputs_.height; }
    long seed() const { UW_FIELD_READ(seed); return inputs_.seed; }
    bool rotation() const { UW_FIELD_READ(prograde_rotation); return inputs_.rotation; }
    float tilt() const { UW_FIELD_READ(axial_tilt); return inputs_.tilt; }
    float eccentricity() const { UW_FIELD_READ(eccentricity); return inputs_.eccentricity; }
    int perihelion() const { UW_FIELD_READ(perihelion); return inputs_.perihelion; }
    float gravity() const { UW_FIELD_READ(gravity); return inputs_.gravity; }
    float lunar() const { UW_FIELD_READ(lunar_strength); return inputs_.lunar; }
    float tempdecrease() const { UW_FIELD_READ(lapse_rate); return inputs_.tempdecrease; }
    int northpolaradjust() const { UW_FIELD_READ(north_polar_adjust); return inputs_.northpolaradjust; }
    int southpolaradjust() const { UW_FIELD_READ(south_polar_adjust); return inputs_.southpolaradjust; }
    int averagetemp() const { UW_FIELD_READ(average_temperature); return inputs_.averagetemp; }
    float waterpickup() const { UW_FIELD_READ(water_pickup); return inputs_.waterpickup; }
    int maxelevation() const { UW_FIELD_READ(max_elevation); return inputs_.maxelevation; }
    int sealevel() const { UW_FIELD_READ(sea_level); return inputs_.sealevel; }
    int landtotal() const { UW_FIELD_READ(land_cells); return inputs_.landtotal; }
    int seatotal() const { UW_FIELD_READ(sea_cells); return inputs_.seatotal; }
    int nom(int x, int y) const;    // whole metres
    int nom(mesh::CellId c) const;
    int map(int x, int y) const;    // whole metres
    int map(mesh::CellId c) const;
    float nomfloat(int x, int y) const;
    float nomfloat(mesh::CellId c) const;
    float mapfloat(int x, int y) const;
    float mapfloat(mesh::CellId c) const;
    bool sea(int x, int y) const;
    bool sea(mesh::CellId c) const;
    bool seawrap(int x, int y) const;
    bool basebelowsea(int x, int y) const;   // the base elevation in whole metres is at or below sea level (planet::basebelowsea) and the cell has no lake surface
    bool basebelowsea(mesh::CellId c) const;
    bool outline(int x, int y) const;
    int mountainheight(int x, int y) const;
    int mountainheight(mesh::CellId c) const;
    float mountainheightfloat(int x, int y) const;
    float mountainheightfloat(mesh::CellId c) const;
    int craterrim(int x, int y) const;
    int craterrim(mesh::CellId c) const;
    float craterrimfloat(int x, int y) const;
    float craterrimfloat(mesh::CellId c) const;
    int volcano(int x, int y) const;
    int volcano(mesh::CellId c) const;
    float volcanofloat(int x, int y) const;
    float volcanofloat(mesh::CellId c) const;
    int lakesurface(int x, int y) const;
    int lakesurface(mesh::CellId c) const;
    float lakesurfacefloat(int x, int y) const;
    float lakesurfacefloat(mesh::CellId c) const;
    int riftlakesurface(int x, int y) const;
    int riftlakesurface(mesh::CellId c) const;
    float riftlakesurfacefloat(int x, int y) const;
    float riftlakesurfacefloat(mesh::CellId c) const;
    int special(int x, int y) const;
    int special(mesh::CellId c) const;
    int truelake(int x, int y) const;
    int truelake(mesh::CellId c) const;

    // Climate state.
    int northpolartemp() const;
    void setnorthpolartemp(int amount);
    int southpolartemp() const;
    void setsouthpolartemp(int amount);
    int eqtemp() const;
    void seteqtemp(int amount);
    // The former int16 encoding of a whole field, cell by cell in storage order:
    // for the .uww save format and the renderers.
    static std::vector<short> legacyfield(const std::vector<float>& values, float scale = 1.0f);
    static std::vector<short> legacyfield(const worldgrid<float>& values);
    std::vector<short> legacyjantempmap() const;
    std::vector<short> legacyjultempmap() const;
    std::vector<short> legacyjanrainmap() const;
    std::vector<short> legacyjulrainmap() const;
    const short* rawseaicemap() const;
    const short* rawclimatemap() const;
    const short* rawbiomemap() const;
    int jantemp(int x, int y) const;
    int jantemp(mesh::CellId c) const;
    void setjantemp(int x, int y, int amount);
    void setjantemp(mesh::CellId c, int amount);
    int jultemp(int x, int y) const;
    int jultemp(mesh::CellId c) const;
    void setjultemp(int x, int y, int amount);
    void setjultemp(mesh::CellId c, int amount);
    float jantempfloat(int x, int y) const;
    float jantempfloat(mesh::CellId c) const;
    void setjantempfloat(int x, int y, float amount);
    void setjantempfloat(mesh::CellId c, float amount);
    float jultempfloat(int x, int y) const;
    float jultempfloat(mesh::CellId c) const;
    void setjultempfloat(int x, int y, float amount);
    void setjultempfloat(mesh::CellId c, float amount);
    int seasonaltemp(int season, int x, int y) const;
    int seasonaltemp(int season, mesh::CellId c) const;
    void setseasonaltemp(int season, int x, int y, int amount);
    void setseasonaltemp(int season, mesh::CellId c, int amount);
    float seasonaltempfloat(int season, int x, int y) const;
    float seasonaltempfloat(int season, mesh::CellId c) const;
    void setseasonaltempfloat(int season, int x, int y, float amount);
    void setseasonaltempfloat(int season, mesh::CellId c, float amount);
    int janrain(int x, int y) const;
    int janrain(mesh::CellId c) const;
    void setjanrain(int x, int y, int amount);
    void setjanrain(mesh::CellId c, int amount);
    int julrain(int x, int y) const;
    int julrain(mesh::CellId c) const;
    void setjulrain(int x, int y, int amount);
    void setjulrain(mesh::CellId c, int amount);
    float janrainfloat(int x, int y) const;
    float janrainfloat(mesh::CellId c) const;
    void setjanrainfloat(int x, int y, float amount);
    void setjanrainfloat(mesh::CellId c, float amount);
    float julrainfloat(int x, int y) const;
    float julrainfloat(mesh::CellId c) const;
    void setjulrainfloat(int x, int y, float amount);
    void setjulrainfloat(mesh::CellId c, float amount);
    int seasonalrain(int season, int x, int y) const;
    int seasonalrain(int season, mesh::CellId c) const;
    void setseasonalrain(int season, int x, int y, int amount);
    void setseasonalrain(int season, mesh::CellId c, int amount);
    float seasonalrainfloat(int season, int x, int y) const;
    float seasonalrainfloat(int season, mesh::CellId c) const;
    void setseasonalrainfloat(int season, int x, int y, float amount);
    void setseasonalrainfloat(int season, mesh::CellId c, float amount);
    int seasonalpressure(int season, int x, int y) const;
    int seasonalpressure(int season, mesh::CellId c) const;
    void setseasonalpressure(int season, int x, int y, int amount);
    void setseasonalpressure(int season, mesh::CellId c, int amount);
    float seasonalpressurefloat(int season, int x, int y) const;
    float seasonalpressurefloat(int season, mesh::CellId c) const;
    void setseasonalpressurefloat(int season, int x, int y, float amount);
    void setseasonalpressurefloat(int season, mesh::CellId c, float amount);
    int seasonaluwind(int season, int x, int y) const;
    int seasonaluwind(int season, mesh::CellId c) const;
    void setseasonaluwind(int season, int x, int y, int amount);
    void setseasonaluwind(int season, mesh::CellId c, int amount);
    float seasonaluwindfloat(int season, int x, int y) const;
    float seasonaluwindfloat(int season, mesh::CellId c) const;
    void setseasonaluwindfloat(int season, int x, int y, float amount);
    void setseasonaluwindfloat(int season, mesh::CellId c, float amount);
    int seasonalvwind(int season, int x, int y) const;
    int seasonalvwind(int season, mesh::CellId c) const;
    void setseasonalvwind(int season, int x, int y, int amount);
    void setseasonalvwind(int season, mesh::CellId c, int amount);
    float seasonalvwindfloat(int season, int x, int y) const;
    float seasonalvwindfloat(int season, mesh::CellId c) const;
    void setseasonalvwindfloat(int season, int x, int y, float amount);
    void setseasonalvwindfloat(int season, mesh::CellId c, float amount);
    int seasonalupperheight(int season, int x, int y) const;
    int seasonalupperheight(int season, mesh::CellId c) const;
    void setseasonalupperheight(int season, int x, int y, int amount);
    void setseasonalupperheight(int season, mesh::CellId c, int amount);
    float seasonalupperheightfloat(int season, int x, int y) const;
    float seasonalupperheightfloat(int season, mesh::CellId c) const;
    void setseasonalupperheightfloat(int season, int x, int y, float amount);
    void setseasonalupperheightfloat(int season, mesh::CellId c, float amount);
    int seasonalupperuwind(int season, int x, int y) const;
    int seasonalupperuwind(int season, mesh::CellId c) const;
    void setseasonalupperuwind(int season, int x, int y, int amount);
    void setseasonalupperuwind(int season, mesh::CellId c, int amount);
    float seasonalupperuwindfloat(int season, int x, int y) const;
    float seasonalupperuwindfloat(int season, mesh::CellId c) const;
    void setseasonalupperuwindfloat(int season, int x, int y, float amount);
    void setseasonalupperuwindfloat(int season, mesh::CellId c, float amount);
    int seasonaluppervwind(int season, int x, int y) const;
    int seasonaluppervwind(int season, mesh::CellId c) const;
    void setseasonaluppervwind(int season, int x, int y, int amount);
    void setseasonaluppervwind(int season, mesh::CellId c, int amount);
    float seasonaluppervwindfloat(int season, int x, int y) const;
    float seasonaluppervwindfloat(int season, mesh::CellId c) const;
    void setseasonaluppervwindfloat(int season, int x, int y, float amount);
    void setseasonaluppervwindfloat(int season, mesh::CellId c, float amount);
    int seasonalverticalvelocity(int season, int x, int y) const;
    int seasonalverticalvelocity(int season, mesh::CellId c) const;
    void setseasonalverticalvelocity(int season, int x, int y, int amount);
    void setseasonalverticalvelocity(int season, mesh::CellId c, int amount);
    float seasonalverticalvelocityfloat(int season, int x, int y) const;
    float seasonalverticalvelocityfloat(int season, mesh::CellId c) const;
    void setseasonalverticalvelocityfloat(int season, int x, int y, float amount);
    void setseasonalverticalvelocityfloat(int season, mesh::CellId c, float amount);
    int seasonalcurrentu(int season, int x, int y) const;
    int seasonalcurrentu(int season, mesh::CellId c) const;
    void setseasonalcurrentu(int season, int x, int y, int amount);
    void setseasonalcurrentu(int season, mesh::CellId c, int amount);
    float seasonalcurrentufloat(int season, int x, int y) const;
    float seasonalcurrentufloat(int season, mesh::CellId c) const;
    void setseasonalcurrentufloat(int season, int x, int y, float amount);
    void setseasonalcurrentufloat(int season, mesh::CellId c, float amount);
    int seasonalcurrentv(int season, int x, int y) const;
    int seasonalcurrentv(int season, mesh::CellId c) const;
    void setseasonalcurrentv(int season, int x, int y, int amount);
    void setseasonalcurrentv(int season, mesh::CellId c, int amount);
    float seasonalcurrentvfloat(int season, int x, int y) const;
    float seasonalcurrentvfloat(int season, mesh::CellId c) const;
    void setseasonalcurrentvfloat(int season, int x, int y, float amount);
    void setseasonalcurrentvfloat(int season, mesh::CellId c, float amount);
    int seasonalsst(int season, int x, int y) const;
    int seasonalsst(int season, mesh::CellId c) const;
    void setseasonalsst(int season, int x, int y, int amount);
    void setseasonalsst(int season, mesh::CellId c, int amount);
    float seasonalsstfloat(int season, int x, int y) const;
    float seasonalsstfloat(int season, mesh::CellId c) const;
    void setseasonalsstfloat(int season, int x, int y, float amount);
    void setseasonalsstfloat(int season, mesh::CellId c, float amount);
    int seasonalevaporation(int season, int x, int y) const;
    int seasonalevaporation(int season, mesh::CellId c) const;
    void setseasonalevaporation(int season, int x, int y, int amount);
    void setseasonalevaporation(int season, mesh::CellId c, int amount);
    float seasonalevaporationfloat(int season, int x, int y) const;
    float seasonalevaporationfloat(int season, mesh::CellId c) const;
    void setseasonalevaporationfloat(int season, int x, int y, float amount);
    void setseasonalevaporationfloat(int season, mesh::CellId c, float amount);
    // Potential evapotranspiration of the rainfall solver (mm/month; erosion plan C2).
    float seasonalpotentialevaporationfloat(int season, int x, int y) const;
    float seasonalpotentialevaporationfloat(int season, mesh::CellId c) const;
    void setseasonalpotentialevaporationfloat(int season, int x, int y, float amount);
    void setseasonalpotentialevaporationfloat(int season, mesh::CellId c, float amount);
    int seasonalmaritimeinfluence(int season, int x, int y) const;
    int seasonalmaritimeinfluence(int season, mesh::CellId c) const;
    void setseasonalmaritimeinfluence(int season, int x, int y, int amount);
    void setseasonalmaritimeinfluence(int season, mesh::CellId c, int amount);
    float seasonalmaritimeinfluencefloat(int season, int x, int y) const;
    float seasonalmaritimeinfluencefloat(int season, mesh::CellId c) const;
    void setseasonalmaritimeinfluencefloat(int season, int x, int y, float amount);
    void setseasonalmaritimeinfluencefloat(int season, mesh::CellId c, float amount);
    int seasonalmaritimethermalanomaly(int season, int x, int y) const;
    int seasonalmaritimethermalanomaly(int season, mesh::CellId c) const;
    void setseasonalmaritimethermalanomaly(int season, int x, int y, int amount);
    void setseasonalmaritimethermalanomaly(int season, mesh::CellId c, int amount);
    float seasonalmaritimethermalanomalyfloat(int season, int x, int y) const;
    float seasonalmaritimethermalanomalyfloat(int season, mesh::CellId c) const;
    void setseasonalmaritimethermalanomalyfloat(int season, int x, int y, float amount);
    void setseasonalmaritimethermalanomalyfloat(int season, mesh::CellId c, float amount);
    int seasonalmaritimefetch(int season, int x, int y) const;
    int seasonalmaritimefetch(int season, mesh::CellId c) const;
    void setseasonalmaritimefetch(int season, int x, int y, int amount);
    void setseasonalmaritimefetch(int season, mesh::CellId c, int amount);
    int seasonalmoisture(int season, int x, int y) const;
    int seasonalmoisture(int season, mesh::CellId c) const;
    void setseasonalmoisture(int season, int x, int y, int amount);
    void setseasonalmoisture(int season, mesh::CellId c, int amount);
    float seasonalmoisturefloat(int season, int x, int y) const;
    float seasonalmoisturefloat(int season, mesh::CellId c) const;
    void setseasonalmoisturefloat(int season, int x, int y, float amount);
    void setseasonalmoisturefloat(int season, mesh::CellId c, float amount);
    int seasonalconvergence(int season, int x, int y) const;
    int seasonalconvergence(int season, mesh::CellId c) const;
    void setseasonalconvergence(int season, int x, int y, int amount);
    void setseasonalconvergence(int season, mesh::CellId c, int amount);
    float seasonalconvergencefloat(int season, int x, int y) const;
    float seasonalconvergencefloat(int season, mesh::CellId c) const;
    void setseasonalconvergencefloat(int season, int x, int y, float amount);
    void setseasonalconvergencefloat(int season, mesh::CellId c, float amount);
    int seasonaluplift(int season, int x, int y) const;
    int seasonaluplift(int season, mesh::CellId c) const;
    void setseasonaluplift(int season, int x, int y, int amount);
    void setseasonaluplift(int season, mesh::CellId c, int amount);
    float seasonalupliftfloat(int season, int x, int y) const;
    float seasonalupliftfloat(int season, mesh::CellId c) const;
    void setseasonalupliftfloat(int season, int x, int y, float amount);
    void setseasonalupliftfloat(int season, mesh::CellId c, float amount);
    int seasonalsubsidence(int season, int x, int y) const;
    int seasonalsubsidence(int season, mesh::CellId c) const;
    void setseasonalsubsidence(int season, int x, int y, int amount);
    void setseasonalsubsidence(int season, mesh::CellId c, int amount);
    float seasonalsubsidencefloat(int season, int x, int y) const;
    float seasonalsubsidencefloat(int season, mesh::CellId c) const;
    void setseasonalsubsidencefloat(int season, int x, int y, float amount);
    void setseasonalsubsidencefloat(int season, mesh::CellId c, float amount);
    int janmountainrain(int x, int y) const;
    int janmountainrain(mesh::CellId c) const;
    void setjanmountainrain(int x, int y, int amount);
    void setjanmountainrain(mesh::CellId c, int amount);
    int julmountainrain(int x, int y) const;
    int julmountainrain(mesh::CellId c) const;
    void setjulmountainrain(int x, int y, int amount);
    void setjulmountainrain(mesh::CellId c, int amount);
    int wintermountainrain(int x, int y) const;
    int wintermountainrain(mesh::CellId c) const;
    void setwintermountainrain(int x, int y, int amount);
    void setwintermountainrain(mesh::CellId c, int amount);
    int summermountainrain(int x, int y) const;
    int summermountainrain(mesh::CellId c) const;
    void setsummermountainrain(int x, int y, int amount);
    void setsummermountainrain(mesh::CellId c, int amount);
    void setwintermountainrainfloat(int x, int y, float amount);
    void setwintermountainrainfloat(mesh::CellId c, float amount);
    void setsummermountainrainfloat(int x, int y, float amount);
    void setsummermountainrainfloat(mesh::CellId c, float amount);
    int janmountainraindir(int x, int y) const;
    int janmountainraindir(mesh::CellId c) const;
    void setjanmountainraindir(int x, int y, int amount);
    void setjanmountainraindir(mesh::CellId c, int amount);
    int julmountainraindir(int x, int y) const;
    int julmountainraindir(mesh::CellId c) const;
    void setjulmountainraindir(int x, int y, int amount);
    void setjulmountainraindir(mesh::CellId c, int amount);
    int wintermountainraindir(int x, int y) const;
    int wintermountainraindir(mesh::CellId c) const;
    int summermountainraindir(int x, int y) const;
    int summermountainraindir(mesh::CellId c) const;
    void setwintermountainraindir(int x, int y, int amount);
    void setwintermountainraindir(mesh::CellId c, int amount);
    void setsummermountainraindir(int x, int y, int amount);
    void setsummermountainraindir(mesh::CellId c, int amount);
    int climate(int x, int y) const;
    int climate(mesh::CellId c) const;
    void setclimate(int x, int y, int amount);
    void setclimate(mesh::CellId c, int amount);
    int biome(int x, int y) const;
    int biome(mesh::CellId c) const;
    void setbiome(int x, int y, int amount);
    void setbiome(mesh::CellId c, int amount);
    int seaice(int x, int y) const;
    int seaice(mesh::CellId c) const;
    void setseaice(int x, int y, int amount);
    void setseaice(mesh::CellId c, int amount);
    int wind(int x, int y) const;
    int wind(mesh::CellId c) const;
    void setwind(int x, int y, int amount);
    void setwind(mesh::CellId c, int amount);
    int horse(int x, int y) const;
    void sethorse(int x, int y, int amount);
    int mintemp(int x, int y) const;
    int mintemp(mesh::CellId c) const;
    void setmintemp(int x, int y, int amount);
    void setmintemp(mesh::CellId c, int amount);
    int maxtemp(int x, int y) const;
    int maxtemp(mesh::CellId c) const;
    void setmaxtemp(int x, int y, int amount);
    void setmaxtemp(mesh::CellId c, int amount);
    int winterrain(int x, int y) const;
    int winterrain(mesh::CellId c) const;
    void setwinterrain(int x, int y, int amount);
    void setwinterrain(mesh::CellId c, int amount);
    int summerrain(int x, int y) const;
    int summerrain(mesh::CellId c) const;
    void setsummerrain(int x, int y, int amount);
    void setsummerrain(mesh::CellId c, int amount);
    // The same selections on the stored floats.
    float mintempfloat(int x, int y) const;
    float mintempfloat(mesh::CellId c) const;
    void setmintempfloat(int x, int y, float amount);
    void setmintempfloat(mesh::CellId c, float amount);
    float maxtempfloat(int x, int y) const;
    float maxtempfloat(mesh::CellId c) const;
    void setmaxtempfloat(int x, int y, float amount);
    void setmaxtempfloat(mesh::CellId c, float amount);
    float winterrainfloat(int x, int y) const;
    float winterrainfloat(mesh::CellId c) const;
    void setwinterrainfloat(int x, int y, float amount);
    void setwinterrainfloat(mesh::CellId c, float amount);
    float summerrainfloat(int x, int y) const;
    float summerrainfloat(mesh::CellId c) const;
    void setsummerrainfloat(int x, int y, float amount);
    void setsummerrainfloat(mesh::CellId c, float amount);

    // Annual and transition-season views of the seasonal fields.
    int avetemp(int x, int y) const;
    int avetemp(mesh::CellId c) const;
    float avetempfloat(int x, int y) const;
    float avetempfloat(mesh::CellId c) const;
    int averain(int x, int y) const;
    int averain(mesh::CellId c) const;
    float averainfloat(int x, int y) const;
    float averainfloat(mesh::CellId c) const;
    int aprtemp(int x, int y) const;
    int aprtemp(mesh::CellId c) const;
    int octtemp(int x, int y) const;
    int octtemp(mesh::CellId c) const;
    int aprrain(int x, int y) const;
    int aprrain(mesh::CellId c) const;
    int octrain(int x, int y) const;
    int octrain(mesh::CellId c) const;

    // Atmosphere column. updatelevels() copies the model's level fields after a
    // wind solve; consumers read levels() instead of model-specific fields.
    const std::vector<atmospherelevel>& levels() const { UW_FIELD_READ(level_wind_east); return levels_; }
    void updatelevels();

    // Diagnostics of the last solve; null when none was recorded.
    const climatediagnostics* diagnostics() const { return diagnostics_.get(); }
    // Season 0 starts a new annual record; layer 0 is the surface, 1 the upper level.
    void recordprocessdiagnostics(int season, const climatehydrology::SeasonalProcessFields& fields);
    void recordoceandiagnostics(int season, int columns, int rows, const climateocean::OceanState& state);
    void recordweatherdiagnostics(int season, int columns, int rows,
        const climateweather::WeatherStatistics& statistics, int layer);

    // Seasonal state from the legacy January/July fields when the transition seasons are empty.
    void syncseasonalclimatefromlegacy();
    const std::vector<std::uint8_t>& weatheranomalystate() const;
    void setweatheranomalystate(std::vector<std::uint8_t> state);
    bool advanceweather(float elapsedSeconds, int horizontalCells = 0);

private:
    mutable climateinputs inputs_;
    std::vector<atmospherelevel> levels_;

    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalconvergencemaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalcurrentumaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalcurrentvmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalevaporationmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalpotentialevaporationmaps;
    std::array<std::vector<short>, CLIMATESEASONCOUNT> seasonalmaritimefetchmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalmaritimeinfluencemaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalmaritimethermalanomalymaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalmoisturemaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalpressuremaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalrainmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalsstmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalsubsidencemaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonaltempmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalupliftmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalupperheightmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalupperuwindmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonaluppervwindmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonaluwindmaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalverticalvelocitymaps;
    std::array<std::vector<float>, CLIMATESEASONCOUNT> seasonalvwindmaps;
    worldgrid<float> jantempmap;
    worldgrid<float> jultempmap;
    worldgrid<short> climatemap;
    worldgrid<short> biomemap;
    worldgrid<float> janrainmap;
    worldgrid<float> julrainmap;
    worldgrid<float> janmountainrainmap;
    worldgrid<float> julmountainrainmap;
    worldgrid<short> janmountainraindirmap;
    worldgrid<short> julmountainraindirmap;
    worldgrid<short> seaicemap;
    worldgrid<int> windmap;
    worldgrid<int> horselats;
    std::vector<std::uint8_t> weatheranomalystatebytes;
    climatecontext::Context climatecontext_;
    transientstate<climatediagnostics> diagnostics_;
    int itsnorthpolartemp = 0;
    int itssouthpolartemp = 0;
    int itseqtemp = 0;

    // A value in its former int16 encoding: what the integer accessors return.
    static int legacyvalue(float value, float scale = 1.0f)
    {
        return static_cast<int>(std::lround(std::clamp(value * scale, -32768.0f, 32767.0f)));
    }

    int seasonalclimateindex(int x, int y) const { return y * (inputs_.width + 1) + x; }
    bool validseasonindex(int season) const { return season >= 0 && season < CLIMATESEASONCOUNT; }
    // The cell id of a lat-lon (x, y) inside the grid (the world grid's id order).
    mesh::CellId cellindex(int x, int y) const { return static_cast<mesh::CellId>(y * (inputs_.width + 1) + x); }
    std::size_t inputindex(int x, int y) const { return cellindex(x, y); }
    void resizecells(int width, int height, const mesh::CellGrid* grid);
    bool inside(int x, int y) const { return !(y < 0 || y > inputs_.height || x < 0 || x > inputs_.width); }
    // Same wrapping as planet::wrapx, which steps by the largest index.
    int wrapx(int x) const
    {
        while (x > inputs_.width)
            x = x - inputs_.width;
        while (x < 0)
            x = x + inputs_.width;
        return x;
    }
    int clipy(int y) const { return y < 0 ? 0 : (y > inputs_.height ? inputs_.height : y); }
    void noteclimatesolverread() const;
};

inline void climateworld::noteclimatesolverread() const
{
    UW_FIELD_READ(solver_monthly_temperature);
    UW_FIELD_READ(solver_monthly_precipitation);
    UW_FIELD_READ(solver_surface_wind_east);
    UW_FIELD_READ(solver_surface_wind_south);
    UW_FIELD_READ(solver_upper_wind_east);
    UW_FIELD_READ(solver_upper_wind_south);
    UW_FIELD_READ(solver_low_wind_east);
    UW_FIELD_READ(solver_low_wind_south);
    UW_FIELD_READ(solver_sea_surface_temperature);
    UW_FIELD_READ(solver_sea_ice);
}

inline int climateworld::nom(mesh::CellId c) const
{
    UW_FIELD_READ(base_elevation);
    return wholemetres(inputs_.nom[c]);
}

inline int climateworld::nom(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    return inside(x, y) ? wholemetres(inputs_.nom[inputindex(x, y)]) : 0;
}

inline float climateworld::nomfloat(mesh::CellId c) const
{
    UW_FIELD_READ(base_elevation);
    return inputs_.nom[c];
}

inline float climateworld::nomfloat(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    return inside(x, y) ? inputs_.nom[inputindex(x, y)] : 0.0f;
}

// The base elevation plus what map() adds to it; the land/sea test stays on whole metres.
inline float climateworld::mapfloat(mesh::CellId c) const
{
    UW_FIELD_READ(volcano);
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(ocean_ridge_height);
    UW_FIELD_READ(mountain_height);
    UW_FIELD_READ(extra_elevation);
    UW_FIELD_READ(crater_rims);
    UW_FIELD_READ(crater_centres);
    const std::size_t cell = c;
    const float base = inputs_.nom[cell];
    const float thisvolcano = std::fabs(inputs_.volcano[cell]);

    if (wholemetres(base) <= inputs_.sealevel)
        return base + inputs_.oceanridgeheight[cell] + thisvolcano;

    return base + (std::max)(thisvolcano, inputs_.mountainheights[cell]) + inputs_.extraelev[cell] + inputs_.craterrims[cell] + inputs_.cratercentres[cell];
}

inline float climateworld::mapfloat(int x, int y) const
{
    if (!inside(x, y))
        return 0.0f;

    return mapfloat(cellindex(x, y));
}

// Whole metres: the base and what is added to it are rounded separately, as planet::map does.
inline int climateworld::map(mesh::CellId c) const
{
    const float base = inputs_.nom[c];
    return wholemetres(base) + wholemetres(mapfloat(c) - base);
}

inline int climateworld::map(int x, int y) const
{
    if (!inside(x, y))
        return 0;

    return map(cellindex(x, y));
}

inline bool climateworld::sea(mesh::CellId c) const
{
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    UW_FIELD_READ(volcano);
    const std::size_t cell = c;
    if (wholemetres(inputs_.nom[cell]) <= inputs_.sealevel && inputs_.lakesurface[cell] == 0)
        return 1;

    if (inputs_.volcano[cell] == 0)
        return 0;

    const float thisvolcano = std::fabs(inputs_.volcano[cell]);

    if (wholemetres(inputs_.nom[cell] + thisvolcano) <= inputs_.sealevel)
        return 1;
    else
        return 0;
}

inline bool climateworld::sea(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    UW_FIELD_READ(volcano);
    if (!inside(x, y))
        return 0;

    return sea(cellindex(x, y));
}

inline bool climateworld::seawrap(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    x = wrapx(x);
    y = clipy(y);

    const std::size_t cell = inputindex(x, y);
    if (wholemetres(inputs_.nom[cell]) <= inputs_.sealevel && inputs_.lakesurface[cell] == 0)
        return 1;
    else
        return 0;
}

// A cell with a lake surface (a lake or salt pan in an enclosed below-sea basin,
// erosion plan H3) is not below the sea, as in sea().
inline bool climateworld::basebelowsea(mesh::CellId c) const
{
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    if (nom(c) > inputs_.sealevel)
        return false;
    return inputs_.lakesurface[c] == 0;
}

inline bool climateworld::basebelowsea(int x, int y) const
{
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    if (nom(x, y) > inputs_.sealevel)
        return false;
    return !inside(x, y) || inputs_.lakesurface[inputindex(x, y)] == 0;
}

inline bool climateworld::outline(int x, int y) const
{
    if (y < 1 || y > inputs_.height - 1)
        return 0;

    if (sea(x, y) == 0)
    {
        if (seawrap(x - 1, y) == 1)
            return 1;

        if (seawrap(x, y - 1) == 1)
            return 1;

        if (seawrap(x + 1, y) == 1)
            return 1;

        if (seawrap(x, y + 1) == 1)
            return 1;
    }

    return 0;
}

inline int climateworld::mountainheight(mesh::CellId c) const
{
    UW_FIELD_READ(mountain_height);
    return wholemetres(inputs_.mountainheights[c]);
}

inline int climateworld::mountainheight(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    return inside(x, y) ? wholemetres(inputs_.mountainheights[inputindex(x, y)]) : 0;
}

inline float climateworld::mountainheightfloat(mesh::CellId c) const
{
    UW_FIELD_READ(mountain_height);
    return inputs_.mountainheights[c];
}

inline float climateworld::mountainheightfloat(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    return inside(x, y) ? inputs_.mountainheights[inputindex(x, y)] : 0.0f;
}

inline int climateworld::craterrim(mesh::CellId c) const
{
    UW_FIELD_READ(crater_rims);
    return wholemetres(inputs_.craterrims[c]);
}

inline int climateworld::craterrim(int x, int y) const
{
    UW_FIELD_READ(crater_rims);
    return inside(x, y) ? wholemetres(inputs_.craterrims[inputindex(x, y)]) : 0;
}

inline float climateworld::craterrimfloat(mesh::CellId c) const
{
    UW_FIELD_READ(crater_rims);
    return inputs_.craterrims[c];
}

inline float climateworld::craterrimfloat(int x, int y) const
{
    UW_FIELD_READ(crater_rims);
    return inside(x, y) ? inputs_.craterrims[inputindex(x, y)] : 0.0f;
}

inline int climateworld::volcano(mesh::CellId c) const
{
    UW_FIELD_READ(volcano);
    return wholemetres(inputs_.volcano[c]);
}

inline int climateworld::volcano(int x, int y) const
{
    UW_FIELD_READ(volcano);
    return inside(x, y) ? wholemetres(inputs_.volcano[inputindex(x, y)]) : 0;
}

inline float climateworld::volcanofloat(mesh::CellId c) const
{
    UW_FIELD_READ(volcano);
    return inputs_.volcano[c];
}

inline float climateworld::volcanofloat(int x, int y) const
{
    UW_FIELD_READ(volcano);
    return inside(x, y) ? inputs_.volcano[inputindex(x, y)] : 0.0f;
}

inline int climateworld::lakesurface(mesh::CellId c) const
{
    UW_FIELD_READ(lake_surface);
    return wholemetres(inputs_.lakesurface[c]);
}

inline int climateworld::lakesurface(int x, int y) const
{
    UW_FIELD_READ(lake_surface);
    return inside(x, y) ? wholemetres(inputs_.lakesurface[inputindex(x, y)]) : 0;
}

inline float climateworld::lakesurfacefloat(mesh::CellId c) const
{
    UW_FIELD_READ(lake_surface);
    return inputs_.lakesurface[c];
}

inline float climateworld::lakesurfacefloat(int x, int y) const
{
    UW_FIELD_READ(lake_surface);
    return inside(x, y) ? inputs_.lakesurface[inputindex(x, y)] : 0.0f;
}

inline int climateworld::riftlakesurface(mesh::CellId c) const
{
    UW_FIELD_READ(rift_lake_surface);
    return wholemetres(inputs_.riftlakesurface[c]);
}

inline int climateworld::riftlakesurface(int x, int y) const
{
    UW_FIELD_READ(rift_lake_surface);
    return inside(x, y) ? wholemetres(inputs_.riftlakesurface[inputindex(x, y)]) : 0;
}

inline float climateworld::riftlakesurfacefloat(mesh::CellId c) const
{
    UW_FIELD_READ(rift_lake_surface);
    return inputs_.riftlakesurface[c];
}

inline float climateworld::riftlakesurfacefloat(int x, int y) const
{
    UW_FIELD_READ(rift_lake_surface);
    return inside(x, y) ? inputs_.riftlakesurface[inputindex(x, y)] : 0.0f;
}

inline int climateworld::special(mesh::CellId c) const
{
    UW_FIELD_READ(special_features);
    return static_cast<int>(inputs_.special[c]);
}

inline int climateworld::special(int x, int y) const
{
    UW_FIELD_READ(special_features);
    return inside(x, y) ? static_cast<int>(inputs_.special[inputindex(x, y)]) : 0;
}

inline int climateworld::truelake(mesh::CellId c) const
{
    if (lakesurface(c) != 0 && special(c) < 110)
        return (1);
    else
        return(0);
}

inline int climateworld::truelake(int x, int y) const
{
    if (!inside(x, y))
        return 0;

    return truelake(cellindex(x, y));
}

inline int climateworld::northpolartemp() const { UW_FIELD_READ(north_polar_temperature); return itsnorthpolartemp; }

inline void climateworld::setnorthpolartemp(int amount) { itsnorthpolartemp = amount; }

inline int climateworld::southpolartemp() const { UW_FIELD_READ(south_polar_temperature); return itssouthpolartemp; }

inline void climateworld::setsouthpolartemp(int amount) { itssouthpolartemp = amount; }

inline int climateworld::eqtemp() const { UW_FIELD_READ(equatorial_temperature); return itseqtemp; }

inline void climateworld::seteqtemp(int amount) { itseqtemp = amount; }

inline std::vector<short> climateworld::legacyjantempmap() const { UW_FIELD_READ(temperature_january); return legacyfield(jantempmap); }

inline std::vector<short> climateworld::legacyjultempmap() const { UW_FIELD_READ(temperature_july); return legacyfield(jultempmap); }

inline std::vector<short> climateworld::legacyjanrainmap() const { UW_FIELD_READ(precipitation_january); return legacyfield(janrainmap); }

inline std::vector<short> climateworld::legacyjulrainmap() const { UW_FIELD_READ(precipitation_july); return legacyfield(julrainmap); }

inline const short* climateworld::rawseaicemap() const { UW_FIELD_READ(sea_ice); return seaicemap.data(); }

inline const short* climateworld::rawclimatemap() const { UW_FIELD_READ(climate_class); return climatemap.data(); }

inline const short* climateworld::rawbiomemap() const { UW_FIELD_READ(biome); return biomemap.data(); }

inline int climateworld::jantemp(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    return legacyvalue(jantempmap.at(c));
}

inline int climateworld::jantemp(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return jantemp(cellindex(x, y));
}

inline void climateworld::setjantemp(mesh::CellId c, int amount)
{
    jantempmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setjantemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjantemp(cellindex(x, y), amount);
}

inline int climateworld::jultemp(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_july);
    return legacyvalue(jultempmap.at(c));
}

inline int climateworld::jultemp(int x, int y) const
{
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return jultemp(cellindex(x, y));
}

inline void climateworld::setjultemp(mesh::CellId c, int amount)
{
    jultempmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setjultemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjultemp(cellindex(x, y), amount);
}

inline int climateworld::seasonaltemp(int season, mesh::CellId c) const
{
    UW_FIELD_READ(temperature);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonaltempmaps[season][c]);
}

inline int climateworld::seasonaltemp(int season, int x, int y) const
{
    UW_FIELD_READ(temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonaltemp(season, cellindex(x, y));
}

inline void climateworld::setseasonaltemp(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonaltempmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonaltemp(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaltemp(season, cellindex(x, y), amount);
}

inline int climateworld::janrain(mesh::CellId c) const
{
    UW_FIELD_READ(precipitation_january);
    return legacyvalue(janrainmap.at(c));
}

inline int climateworld::janrain(int x, int y) const
{
    UW_FIELD_READ(precipitation_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return janrain(cellindex(x, y));
}

inline void climateworld::setjanrain(mesh::CellId c, int amount)
{
    janrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setjanrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjanrain(cellindex(x, y), amount);
}

inline int climateworld::julrain(mesh::CellId c) const
{
    UW_FIELD_READ(precipitation_july);
    return legacyvalue(julrainmap.at(c));
}

inline int climateworld::julrain(int x, int y) const
{
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return julrain(cellindex(x, y));
}

inline void climateworld::setjulrain(mesh::CellId c, int amount)
{
    julrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setjulrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjulrain(cellindex(x, y), amount);
}

inline int climateworld::seasonalrain(int season, mesh::CellId c) const
{
    UW_FIELD_READ(precipitation);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalrainmaps[season][c]);
}

inline int climateworld::seasonalrain(int season, int x, int y) const
{
    UW_FIELD_READ(precipitation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalrain(season, cellindex(x, y));
}

inline void climateworld::setseasonalrain(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    const int index = c;
    seasonalrainmaps[season][index] = static_cast<float>(amount);
}

inline void climateworld::setseasonalrain(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalrain(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalrainfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(precipitation);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalrainmaps[season][c];
}

inline float climateworld::seasonalrainfloat(int season, int x, int y) const
{
    UW_FIELD_READ(precipitation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalrainfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalrainfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    const int index = c;
    const float physicalamount = (std::max)(0.0f, amount);
    seasonalrainmaps[season][index] = physicalamount;
}

inline void climateworld::setseasonalrainfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalrainfloat(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalpressure(int season, mesh::CellId c) const
{
    UW_FIELD_READ(pressure);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalpressuremaps[season][c]);
}

inline int climateworld::seasonalpressure(int season, int x, int y) const
{
    UW_FIELD_READ(pressure);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalpressure(season, cellindex(x, y));
}

inline void climateworld::setseasonalpressure(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalpressuremaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalpressure(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalpressure(season, cellindex(x, y), amount);
}

inline int climateworld::seasonaluwind(int season, mesh::CellId c) const
{
    UW_FIELD_READ(wind_east);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonaluwindmaps[season][c]);
}

inline int climateworld::seasonaluwind(int season, int x, int y) const
{
    UW_FIELD_READ(wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonaluwind(season, cellindex(x, y));
}

inline void climateworld::setseasonaluwind(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonaluwindmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonaluwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaluwind(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalvwind(int season, mesh::CellId c) const
{
    UW_FIELD_READ(wind_south);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalvwindmaps[season][c]);
}

inline int climateworld::seasonalvwind(int season, int x, int y) const
{
    UW_FIELD_READ(wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalvwind(season, cellindex(x, y));
}

inline void climateworld::setseasonalvwind(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalvwindmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalvwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalvwind(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalupperheight(int season, mesh::CellId c) const
{
    UW_FIELD_READ(upper_height);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalupperheightmaps[season][c]);
}

inline int climateworld::seasonalupperheight(int season, int x, int y) const
{
    UW_FIELD_READ(upper_height);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalupperheight(season, cellindex(x, y));
}

inline void climateworld::setseasonalupperheight(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalupperheightmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalupperheight(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalupperheight(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalupperuwind(int season, mesh::CellId c) const
{
    UW_FIELD_READ(upper_wind_east);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalupperuwindmaps[season][c]);
}

inline int climateworld::seasonalupperuwind(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalupperuwind(season, cellindex(x, y));
}

inline void climateworld::setseasonalupperuwind(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalupperuwindmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalupperuwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalupperuwind(season, cellindex(x, y), amount);
}

inline int climateworld::seasonaluppervwind(int season, mesh::CellId c) const
{
    UW_FIELD_READ(upper_wind_south);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonaluppervwindmaps[season][c]);
}

inline int climateworld::seasonaluppervwind(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonaluppervwind(season, cellindex(x, y));
}

inline void climateworld::setseasonaluppervwind(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonaluppervwindmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonaluppervwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaluppervwind(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalverticalvelocity(int season, mesh::CellId c) const
{
    UW_FIELD_READ(vertical_velocity);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalverticalvelocitymaps[season][c], climatelegacyscale::verticalvelocity);
}

inline int climateworld::seasonalverticalvelocity(int season, int x, int y) const
{
    UW_FIELD_READ(vertical_velocity);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalverticalvelocity(season, cellindex(x, y));
}

inline void climateworld::setseasonalverticalvelocity(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalverticalvelocitymaps[season][c] = static_cast<float>(amount) / climatelegacyscale::verticalvelocity;
}

inline void climateworld::setseasonalverticalvelocity(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalverticalvelocity(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalcurrentu(int season, mesh::CellId c) const
{
    UW_FIELD_READ(ocean_current_east);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalcurrentumaps[season][c], climatelegacyscale::current);
}

inline int climateworld::seasonalcurrentu(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalcurrentu(season, cellindex(x, y));
}

inline void climateworld::setseasonalcurrentu(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalcurrentumaps[season][c] = static_cast<float>(amount) / climatelegacyscale::current;
}

inline void climateworld::setseasonalcurrentu(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalcurrentu(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalcurrentv(int season, mesh::CellId c) const
{
    UW_FIELD_READ(ocean_current_south);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalcurrentvmaps[season][c], climatelegacyscale::current);
}

inline int climateworld::seasonalcurrentv(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalcurrentv(season, cellindex(x, y));
}

inline void climateworld::setseasonalcurrentv(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalcurrentvmaps[season][c] = static_cast<float>(amount) / climatelegacyscale::current;
}

inline void climateworld::setseasonalcurrentv(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalcurrentv(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalsst(int season, mesh::CellId c) const
{
    UW_FIELD_READ(sea_surface_temperature);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalsstmaps[season][c]);
}

inline int climateworld::seasonalsst(int season, int x, int y) const
{
    UW_FIELD_READ(sea_surface_temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalsst(season, cellindex(x, y));
}

inline void climateworld::setseasonalsst(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalsstmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalsst(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalsst(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalevaporation(int season, mesh::CellId c) const
{
    UW_FIELD_READ(evaporation);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalevaporationmaps[season][c]);
}

inline int climateworld::seasonalevaporation(int season, int x, int y) const
{
    UW_FIELD_READ(evaporation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalevaporation(season, cellindex(x, y));
}

inline void climateworld::setseasonalevaporation(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalevaporationmaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalevaporation(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalevaporation(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalmaritimeinfluence(int season, mesh::CellId c) const
{
    UW_FIELD_READ(maritime_influence);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalmaritimeinfluencemaps[season][c], climatelegacyscale::maritimeinfluence);
}

inline int climateworld::seasonalmaritimeinfluence(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_influence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalmaritimeinfluence(season, cellindex(x, y));
}

inline void climateworld::setseasonalmaritimeinfluence(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmaritimeinfluencemaps[season][c] = static_cast<float>(amount) / climatelegacyscale::maritimeinfluence;
}

inline void climateworld::setseasonalmaritimeinfluence(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmaritimeinfluence(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalmaritimethermalanomaly(int season, mesh::CellId c) const
{
    UW_FIELD_READ(maritime_thermal_anomaly);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalmaritimethermalanomalymaps[season][c], climatelegacyscale::maritimethermalanomaly);
}

inline int climateworld::seasonalmaritimethermalanomaly(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_thermal_anomaly);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalmaritimethermalanomaly(season, cellindex(x, y));
}

inline void climateworld::setseasonalmaritimethermalanomaly(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmaritimethermalanomalymaps[season][c] = static_cast<float>(amount) / climatelegacyscale::maritimethermalanomaly;
}

inline void climateworld::setseasonalmaritimethermalanomaly(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmaritimethermalanomaly(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalmaritimefetch(int season, mesh::CellId c) const
{
    UW_FIELD_READ(maritime_fetch);
    if (!validseasonindex(season))
        return 0;
    return static_cast<int>(seasonalmaritimefetchmaps[season][c]);
}

inline int climateworld::seasonalmaritimefetch(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_fetch);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalmaritimefetch(season, cellindex(x, y));
}

inline void climateworld::setseasonalmaritimefetch(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmaritimefetchmaps[season][c] = static_cast<short>(amount);
}

inline void climateworld::setseasonalmaritimefetch(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmaritimefetch(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalmoisture(int season, mesh::CellId c) const
{
    UW_FIELD_READ(moisture);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalmoisturemaps[season][c]);
}

inline int climateworld::seasonalmoisture(int season, int x, int y) const
{
    UW_FIELD_READ(moisture);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalmoisture(season, cellindex(x, y));
}

inline void climateworld::setseasonalmoisture(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmoisturemaps[season][c] = static_cast<float>(amount);
}

inline void climateworld::setseasonalmoisture(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmoisture(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalconvergence(int season, mesh::CellId c) const
{
    UW_FIELD_READ(convergence);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalconvergencemaps[season][c], climatelegacyscale::convergence);
}

inline int climateworld::seasonalconvergence(int season, int x, int y) const
{
    UW_FIELD_READ(convergence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalconvergence(season, cellindex(x, y));
}

inline void climateworld::setseasonalconvergence(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalconvergencemaps[season][c] = static_cast<float>(amount) / climatelegacyscale::convergence;
}

inline void climateworld::setseasonalconvergence(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalconvergence(season, cellindex(x, y), amount);
}

inline int climateworld::seasonaluplift(int season, mesh::CellId c) const
{
    UW_FIELD_READ(orographic_uplift);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalupliftmaps[season][c], climatelegacyscale::orographicmotion);
}

inline int climateworld::seasonaluplift(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_uplift);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonaluplift(season, cellindex(x, y));
}

inline void climateworld::setseasonaluplift(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalupliftmaps[season][c] = static_cast<float>(amount) / climatelegacyscale::orographicmotion;
}

inline void climateworld::setseasonaluplift(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaluplift(season, cellindex(x, y), amount);
}

inline int climateworld::seasonalsubsidence(int season, mesh::CellId c) const
{
    UW_FIELD_READ(orographic_subsidence);
    if (!validseasonindex(season))
        return 0;
    return legacyvalue(seasonalsubsidencemaps[season][c], climatelegacyscale::orographicmotion);
}

inline int climateworld::seasonalsubsidence(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_subsidence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seasonalsubsidence(season, cellindex(x, y));
}

inline void climateworld::setseasonalsubsidence(int season, mesh::CellId c, int amount)
{
    if (!validseasonindex(season))
        return;
    seasonalsubsidencemaps[season][c] = static_cast<float>(amount) / climatelegacyscale::orographicmotion;
}

inline void climateworld::setseasonalsubsidence(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalsubsidence(season, cellindex(x, y), amount);
}

// Float accessors: the stored values in physical units.

inline float climateworld::seasonaltempfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(temperature);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonaltempmaps[season][c];
}

inline float climateworld::seasonaltempfloat(int season, int x, int y) const
{
    UW_FIELD_READ(temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonaltempfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonaltempfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonaltempmaps[season][c] = amount;
}

inline void climateworld::setseasonaltempfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaltempfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalpressurefloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(pressure);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalpressuremaps[season][c];
}

inline float climateworld::seasonalpressurefloat(int season, int x, int y) const
{
    UW_FIELD_READ(pressure);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalpressurefloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalpressurefloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalpressuremaps[season][c] = amount;
}

inline void climateworld::setseasonalpressurefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalpressurefloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonaluwindfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(wind_east);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonaluwindmaps[season][c];
}

inline float climateworld::seasonaluwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonaluwindfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonaluwindfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonaluwindmaps[season][c] = amount;
}

inline void climateworld::setseasonaluwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaluwindfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalvwindfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(wind_south);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalvwindmaps[season][c];
}

inline float climateworld::seasonalvwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalvwindfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalvwindfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalvwindmaps[season][c] = amount;
}

inline void climateworld::setseasonalvwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalvwindfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalupperheightfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(upper_height);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalupperheightmaps[season][c];
}

inline float climateworld::seasonalupperheightfloat(int season, int x, int y) const
{
    UW_FIELD_READ(upper_height);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalupperheightfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalupperheightfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalupperheightmaps[season][c] = amount;
}

inline void climateworld::setseasonalupperheightfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalupperheightfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalupperuwindfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(upper_wind_east);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalupperuwindmaps[season][c];
}

inline float climateworld::seasonalupperuwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalupperuwindfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalupperuwindfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalupperuwindmaps[season][c] = amount;
}

inline void climateworld::setseasonalupperuwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalupperuwindfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonaluppervwindfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(upper_wind_south);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonaluppervwindmaps[season][c];
}

inline float climateworld::seasonaluppervwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonaluppervwindfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonaluppervwindfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonaluppervwindmaps[season][c] = amount;
}

inline void climateworld::setseasonaluppervwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonaluppervwindfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalverticalvelocityfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(vertical_velocity);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalverticalvelocitymaps[season][c];
}

inline float climateworld::seasonalverticalvelocityfloat(int season, int x, int y) const
{
    UW_FIELD_READ(vertical_velocity);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalverticalvelocityfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalverticalvelocityfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalverticalvelocitymaps[season][c] = amount;
}

inline void climateworld::setseasonalverticalvelocityfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalverticalvelocityfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalcurrentufloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(ocean_current_east);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalcurrentumaps[season][c];
}

inline float climateworld::seasonalcurrentufloat(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalcurrentufloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalcurrentufloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalcurrentumaps[season][c] = amount;
}

inline void climateworld::setseasonalcurrentufloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalcurrentufloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalcurrentvfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(ocean_current_south);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalcurrentvmaps[season][c];
}

inline float climateworld::seasonalcurrentvfloat(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalcurrentvfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalcurrentvfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalcurrentvmaps[season][c] = amount;
}

inline void climateworld::setseasonalcurrentvfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalcurrentvfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalsstfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(sea_surface_temperature);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalsstmaps[season][c];
}

inline float climateworld::seasonalsstfloat(int season, int x, int y) const
{
    UW_FIELD_READ(sea_surface_temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalsstfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalsstfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalsstmaps[season][c] = amount;
}

inline void climateworld::setseasonalsstfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalsstfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalevaporationfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(evaporation);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalevaporationmaps[season][c];
}

inline float climateworld::seasonalevaporationfloat(int season, int x, int y) const
{
    UW_FIELD_READ(evaporation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalevaporationfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalevaporationfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalevaporationmaps[season][c] = amount;
}

inline void climateworld::setseasonalevaporationfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalevaporationfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalpotentialevaporationfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(potential_evaporation);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalpotentialevaporationmaps[season][c];
}

inline float climateworld::seasonalpotentialevaporationfloat(int season, int x, int y) const
{
    UW_FIELD_READ(potential_evaporation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalpotentialevaporationfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalpotentialevaporationfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalpotentialevaporationmaps[season][c] = amount;
}

inline void climateworld::setseasonalpotentialevaporationfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalpotentialevaporationfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalmaritimeinfluencefloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(maritime_influence);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalmaritimeinfluencemaps[season][c];
}

inline float climateworld::seasonalmaritimeinfluencefloat(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_influence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalmaritimeinfluencefloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalmaritimeinfluencefloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmaritimeinfluencemaps[season][c] = amount;
}

inline void climateworld::setseasonalmaritimeinfluencefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmaritimeinfluencefloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalmaritimethermalanomalyfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(maritime_thermal_anomaly);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalmaritimethermalanomalymaps[season][c];
}

inline float climateworld::seasonalmaritimethermalanomalyfloat(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_thermal_anomaly);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalmaritimethermalanomalyfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalmaritimethermalanomalyfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmaritimethermalanomalymaps[season][c] = amount;
}

inline void climateworld::setseasonalmaritimethermalanomalyfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmaritimethermalanomalyfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalmoisturefloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(moisture);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalmoisturemaps[season][c];
}

inline float climateworld::seasonalmoisturefloat(int season, int x, int y) const
{
    UW_FIELD_READ(moisture);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalmoisturefloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalmoisturefloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalmoisturemaps[season][c] = amount;
}

inline void climateworld::setseasonalmoisturefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalmoisturefloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalconvergencefloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(convergence);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalconvergencemaps[season][c];
}

inline float climateworld::seasonalconvergencefloat(int season, int x, int y) const
{
    UW_FIELD_READ(convergence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalconvergencefloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalconvergencefloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalconvergencemaps[season][c] = amount;
}

inline void climateworld::setseasonalconvergencefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalconvergencefloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalupliftfloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(orographic_uplift);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalupliftmaps[season][c];
}

inline float climateworld::seasonalupliftfloat(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_uplift);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalupliftfloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalupliftfloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalupliftmaps[season][c] = amount;
}

inline void climateworld::setseasonalupliftfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalupliftfloat(season, cellindex(x, y), amount);
}

inline float climateworld::seasonalsubsidencefloat(int season, mesh::CellId c) const
{
    UW_FIELD_READ(orographic_subsidence);
    if (!validseasonindex(season))
        return 0.0f;
    return seasonalsubsidencemaps[season][c];
}

inline float climateworld::seasonalsubsidencefloat(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_subsidence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalsubsidencefloat(season, cellindex(x, y));
}

inline void climateworld::setseasonalsubsidencefloat(int season, mesh::CellId c, float amount)
{
    if (!validseasonindex(season))
        return;
    seasonalsubsidencemaps[season][c] = amount;
}

inline void climateworld::setseasonalsubsidencefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseasonalsubsidencefloat(season, cellindex(x, y), amount);
}

inline float climateworld::jantempfloat(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    return jantempmap.at(c);
}

inline float climateworld::jantempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jantempfloat(cellindex(x, y));
}

inline void climateworld::setjantempfloat(mesh::CellId c, float amount)
{
    jantempmap.at(c) = amount;
}

inline void climateworld::setjantempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjantempfloat(cellindex(x, y), amount);
}

inline float climateworld::jultempfloat(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_july);
    return jultempmap.at(c);
}

inline float climateworld::jultempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jultempfloat(cellindex(x, y));
}

inline void climateworld::setjultempfloat(mesh::CellId c, float amount)
{
    jultempmap.at(c) = amount;
}

inline void climateworld::setjultempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjultempfloat(cellindex(x, y), amount);
}

inline float climateworld::janrainfloat(mesh::CellId c) const
{
    UW_FIELD_READ(precipitation_january);
    return janrainmap.at(c);
}

inline float climateworld::janrainfloat(int x, int y) const
{
    UW_FIELD_READ(precipitation_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return janrainfloat(cellindex(x, y));
}

inline void climateworld::setjanrainfloat(mesh::CellId c, float amount)
{
    janrainmap.at(c) = amount;
}

inline void climateworld::setjanrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjanrainfloat(cellindex(x, y), amount);
}

inline float climateworld::julrainfloat(mesh::CellId c) const
{
    UW_FIELD_READ(precipitation_july);
    return julrainmap.at(c);
}

inline float climateworld::julrainfloat(int x, int y) const
{
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return julrainfloat(cellindex(x, y));
}

inline void climateworld::setjulrainfloat(mesh::CellId c, float amount)
{
    julrainmap.at(c) = amount;
}

inline void climateworld::setjulrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjulrainfloat(cellindex(x, y), amount);
}

inline int climateworld::janmountainrain(mesh::CellId c) const
{
    UW_FIELD_READ(mountain_precipitation_january);
    return legacyvalue(janmountainrainmap.at(c));
}

inline int climateworld::janmountainrain(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return janmountainrain(cellindex(x, y));
}

inline void climateworld::setjanmountainrain(mesh::CellId c, int amount)
{
    janmountainrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setjanmountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjanmountainrain(cellindex(x, y), amount);
}

inline int climateworld::julmountainrain(mesh::CellId c) const
{
    UW_FIELD_READ(mountain_precipitation_july);
    return legacyvalue(julmountainrainmap.at(c));
}

inline int climateworld::julmountainrain(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return julmountainrain(cellindex(x, y));
}

inline void climateworld::setjulmountainrain(mesh::CellId c, int amount)
{
    julmountainrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setjulmountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjulmountainrain(cellindex(x, y), amount);
}

inline int climateworld::wintermountainrain(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_january);
    UW_FIELD_READ(mountain_precipitation_july);
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        return legacyvalue(janmountainrainmap.at(c));
    else
        return legacyvalue(julmountainrainmap.at(c));
}

inline int climateworld::wintermountainrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_january);
    UW_FIELD_READ(mountain_precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return wintermountainrain(cellindex(x, y));
}

inline void climateworld::setwintermountainrain(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        janmountainrainmap.at(c) = static_cast<float>(amount);
    else
        julmountainrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setwintermountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setwintermountainrain(cellindex(x, y), amount);
}

inline int climateworld::summermountainrain(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_january);
    UW_FIELD_READ(mountain_precipitation_july);
    if (legacyvalue(jantempmap.at(c)) >= legacyvalue(jultempmap.at(c)))
        return legacyvalue(janmountainrainmap.at(c));
    else
        return legacyvalue(julmountainrainmap.at(c));
}

inline int climateworld::summermountainrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_january);
    UW_FIELD_READ(mountain_precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return summermountainrain(cellindex(x, y));
}

inline void climateworld::setsummermountainrain(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) >= legacyvalue(jultempmap.at(c)))
        janmountainrainmap.at(c) = static_cast<float>(amount);
    else
        julmountainrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setsummermountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setsummermountainrain(cellindex(x, y), amount);
}

inline void climateworld::setwintermountainrainfloat(mesh::CellId c, float amount)
{
    (jantempmap.at(c) < jultempmap.at(c) ? janmountainrainmap.at(c) : julmountainrainmap.at(c)) = amount;
}

inline void climateworld::setwintermountainrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setwintermountainrainfloat(cellindex(x, y), amount);
}

inline void climateworld::setsummermountainrainfloat(mesh::CellId c, float amount)
{
    (jantempmap.at(c) >= jultempmap.at(c) ? janmountainrainmap.at(c) : julmountainrainmap.at(c)) = amount;
}

inline void climateworld::setsummermountainrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setsummermountainrainfloat(cellindex(x, y), amount);
}

inline int climateworld::janmountainraindir(mesh::CellId c) const
{
    UW_FIELD_READ(mountain_precipitation_direction_january);
    return (int)janmountainraindirmap.at(c);
}

inline int climateworld::janmountainraindir(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_direction_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return janmountainraindir(cellindex(x, y));
}

inline void climateworld::setjanmountainraindir(mesh::CellId c, int amount)
{
    janmountainraindirmap.at(c) = (short)amount;
}

inline void climateworld::setjanmountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjanmountainraindir(cellindex(x, y), amount);
}

inline int climateworld::julmountainraindir(mesh::CellId c) const
{
    UW_FIELD_READ(mountain_precipitation_direction_july);
    return (int)julmountainraindirmap.at(c);
}

inline int climateworld::julmountainraindir(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return julmountainraindir(cellindex(x, y));
}

inline void climateworld::setjulmountainraindir(mesh::CellId c, int amount)
{
    julmountainraindirmap.at(c) = (short)amount;
}

inline void climateworld::setjulmountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setjulmountainraindir(cellindex(x, y), amount);
}

inline int climateworld::wintermountainraindir(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_direction_january);
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        return (int)janmountainraindirmap.at(c);
    else
        return (int)julmountainraindirmap.at(c);
}

inline int climateworld::wintermountainraindir(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_direction_january);
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return wintermountainraindir(cellindex(x, y));
}

inline int climateworld::summermountainraindir(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_direction_january);
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (legacyvalue(jantempmap.at(c)) >= legacyvalue(jultempmap.at(c)))
        return (int)janmountainraindirmap.at(c);
    else
        return (int)julmountainraindirmap.at(c);
}

inline int climateworld::summermountainraindir(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_direction_january);
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return summermountainraindir(cellindex(x, y));
}

inline void climateworld::setwintermountainraindir(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        janmountainraindirmap.at(c) = (short)amount;
    else
        julmountainraindirmap.at(c) = (short)amount;
}

inline void climateworld::setwintermountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setwintermountainraindir(cellindex(x, y), amount);
}

inline void climateworld::setsummermountainraindir(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) >= legacyvalue(jultempmap.at(c)))
        janmountainraindirmap.at(c) = (short)amount;
    else
        julmountainraindirmap.at(c) = (short)amount;
}

inline void climateworld::setsummermountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setsummermountainraindir(cellindex(x, y), amount);
}

inline int climateworld::climate(mesh::CellId c) const
{
    UW_FIELD_READ(climate_class);
    return (int)climatemap.at(c);
}

inline int climateworld::climate(int x, int y) const
{
    UW_FIELD_READ(climate_class);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return climate(cellindex(x, y));
}

inline void climateworld::setclimate(mesh::CellId c, int amount)
{
    climatemap.at(c) = (short)amount;
}

inline void climateworld::setclimate(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setclimate(cellindex(x, y), amount);
}

inline int climateworld::biome(mesh::CellId c) const
{
    UW_FIELD_READ(biome);
    return (int)biomemap.at(c);
}

inline int climateworld::biome(int x, int y) const
{
    UW_FIELD_READ(biome);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return biome(cellindex(x, y));
}

inline void climateworld::setbiome(mesh::CellId c, int amount)
{
    biomemap.at(c) = (short)amount;
}

inline void climateworld::setbiome(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setbiome(cellindex(x, y), amount);
}

inline int climateworld::seaice(mesh::CellId c) const
{
    UW_FIELD_READ(sea_ice);
    return (int)seaicemap.at(c);
}

inline int climateworld::seaice(int x, int y) const
{
    UW_FIELD_READ(sea_ice);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return seaice(cellindex(x, y));
}

inline void climateworld::setseaice(mesh::CellId c, int amount)
{
    seaicemap.at(c) = (short)amount;
}

inline void climateworld::setseaice(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setseaice(cellindex(x, y), amount);
}

inline int climateworld::wind(mesh::CellId c) const
{
    UW_FIELD_READ(prevailing_wind);
    return windmap.at(c);
}

inline int climateworld::wind(int x, int y) const
{
    UW_FIELD_READ(prevailing_wind);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return wind(cellindex(x, y));
}

inline void climateworld::setwind(mesh::CellId c, int amount)
{
    windmap.at(c) = amount;
}

inline void climateworld::setwind(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setwind(cellindex(x, y), amount);
}

inline int climateworld::horse(int x, int y) const
{
    UW_FIELD_READ(subtropical_high_latitudes);
    if (x<0 || x>inputs_.width)
        return 0;

    return horselats[x][y];
}

inline void climateworld::sethorse(int x, int y, int amount)
{
    if (x<0 || x>inputs_.width)
        return;

    horselats[x][y] = amount;
}

inline int climateworld::mintemp(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        return legacyvalue(jantempmap.at(c));
    else
        return legacyvalue(jultempmap.at(c));
}

inline int climateworld::mintemp(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return mintemp(cellindex(x, y));
}

inline void climateworld::setmintemp(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        jantempmap.at(c) = static_cast<float>(amount);
    else
        jultempmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setmintemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setmintemp(cellindex(x, y), amount);
}

inline int climateworld::maxtemp(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (legacyvalue(jantempmap.at(c)) > legacyvalue(jultempmap.at(c)))
        return legacyvalue(jantempmap.at(c));
    else
        return legacyvalue(jultempmap.at(c));
}

inline int climateworld::maxtemp(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return maxtemp(cellindex(x, y));
}

inline void climateworld::setmaxtemp(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) > legacyvalue(jultempmap.at(c)))
        jantempmap.at(c) = static_cast<float>(amount);
    else
        jultempmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setmaxtemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setmaxtemp(cellindex(x, y), amount);
}

inline int climateworld::winterrain(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        return legacyvalue(janrainmap.at(c));
    else
        return legacyvalue(julrainmap.at(c));
}

inline int climateworld::winterrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return winterrain(cellindex(x, y));
}

inline void climateworld::setwinterrain(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) < legacyvalue(jultempmap.at(c)))
        janrainmap.at(c) = static_cast<float>(amount);
    else
        julrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setwinterrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setwinterrain(cellindex(x, y), amount);
}

inline int climateworld::summerrain(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (legacyvalue(jantempmap.at(c)) >= legacyvalue(jultempmap.at(c)))
        return legacyvalue(janrainmap.at(c));
    else
        return legacyvalue(julrainmap.at(c));
}

inline int climateworld::summerrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return summerrain(cellindex(x, y));
}

inline void climateworld::setsummerrain(mesh::CellId c, int amount)
{
    if (legacyvalue(jantempmap.at(c)) >= legacyvalue(jultempmap.at(c)))
        janrainmap.at(c) = static_cast<float>(amount);
    else
        julrainmap.at(c) = static_cast<float>(amount);
}

inline void climateworld::setsummerrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setsummerrain(cellindex(x, y), amount);
}

inline float climateworld::mintempfloat(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    return jantempmap.at(c) < jultempmap.at(c) ? jantempmap.at(c) : jultempmap.at(c);
}

inline float climateworld::mintempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return mintempfloat(cellindex(x, y));
}

inline void climateworld::setmintempfloat(mesh::CellId c, float amount)
{
    (jantempmap.at(c) < jultempmap.at(c) ? jantempmap.at(c) : jultempmap.at(c)) = amount;
}

inline void climateworld::setmintempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setmintempfloat(cellindex(x, y), amount);
}

inline float climateworld::maxtempfloat(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    return jantempmap.at(c) > jultempmap.at(c) ? jantempmap.at(c) : jultempmap.at(c);
}

inline float climateworld::maxtempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return maxtempfloat(cellindex(x, y));
}

inline void climateworld::setmaxtempfloat(mesh::CellId c, float amount)
{
    (jantempmap.at(c) > jultempmap.at(c) ? jantempmap.at(c) : jultempmap.at(c)) = amount;
}

inline void climateworld::setmaxtempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setmaxtempfloat(cellindex(x, y), amount);
}

inline float climateworld::winterrainfloat(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    return jantempmap.at(c) < jultempmap.at(c) ? janrainmap.at(c) : julrainmap.at(c);
}

inline float climateworld::winterrainfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return winterrainfloat(cellindex(x, y));
}

inline void climateworld::setwinterrainfloat(mesh::CellId c, float amount)
{
    (jantempmap.at(c) < jultempmap.at(c) ? janrainmap.at(c) : julrainmap.at(c)) = amount;
}

inline void climateworld::setwinterrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setwinterrainfloat(cellindex(x, y), amount);
}

inline float climateworld::summerrainfloat(mesh::CellId c) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    return jantempmap.at(c) >= jultempmap.at(c) ? janrainmap.at(c) : julrainmap.at(c);
}

inline float climateworld::summerrainfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return summerrainfloat(cellindex(x, y));
}

inline void climateworld::setsummerrainfloat(mesh::CellId c, float amount)
{
    (jantempmap.at(c) >= jultempmap.at(c) ? janrainmap.at(c) : julrainmap.at(c)) = amount;
}

inline void climateworld::setsummerrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    setsummerrainfloat(cellindex(x, y), amount);
}

inline float climateworld::avetempfloat(mesh::CellId c) const
{
    float total = 0.0f;
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += seasonaltempfloat(season, c);
    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

inline float climateworld::avetempfloat(int x, int y) const
{
    float total = 0.0f;
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += seasonaltempfloat(season, x, y);
    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

inline int climateworld::avetemp(mesh::CellId c) const
{
    int total = 0;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total = total + seasonaltemp(season, c);

    return total / CLIMATESEASONCOUNT;
}

inline int climateworld::avetemp(int x, int y) const
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return avetemp(cellindex(x, y));
}

inline int climateworld::averain(mesh::CellId c) const
{
    int total = 0;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total = total + seasonalrain(season, c);

    return total / CLIMATESEASONCOUNT;
}

inline int climateworld::averain(int x, int y) const
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return averain(cellindex(x, y));
}

inline float climateworld::averainfloat(mesh::CellId c) const
{
    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += seasonalrainfloat(season, c);

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

inline float climateworld::averainfloat(int x, int y) const
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return averainfloat(cellindex(x, y));
}

inline int climateworld::aprtemp(mesh::CellId c) const
{
    return seasonaltemp(seasonapril, c);
}

inline int climateworld::aprtemp(int x, int y) const
{
    return seasonaltemp(seasonapril, x, y);
}

inline int climateworld::octtemp(mesh::CellId c) const
{
    return seasonaltemp(seasonoctober, c);
}

inline int climateworld::octtemp(int x, int y) const
{
    return seasonaltemp(seasonoctober, x, y);
}

inline int climateworld::aprrain(mesh::CellId c) const
{
    return seasonalrain(seasonapril, c);
}

inline int climateworld::aprrain(int x, int y) const
{
    return seasonalrain(seasonapril, x, y);
}

inline int climateworld::octrain(mesh::CellId c) const
{
    return seasonalrain(seasonoctober, c);
}

inline int climateworld::octrain(int x, int y) const
{
    return seasonalrain(seasonoctober, x, y);
}
