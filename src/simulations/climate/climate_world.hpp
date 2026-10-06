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


// Terrain and planet inputs of the climate processes. Grids are x-major with a
// stride of `rows` cells over 0..width, 0..height and are owned by the caller.
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
    int rows = ARRAYHEIGHT;
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
// seasonal fields on the world grid, row-major (y * (width + 1) + x), in
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
    // Sizes the seasonal fields for the world grid and resets the solver state.
    void resize(int width, int height);

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
    int map(int x, int y) const;    // whole metres
    float nomfloat(int x, int y) const;
    float mapfloat(int x, int y) const;
    bool sea(int x, int y) const;
    bool seawrap(int x, int y) const;
    bool basebelowsea(int x, int y) const;   // the base elevation in whole metres is at or below sea level (planet::basebelowsea)
    bool outline(int x, int y) const;
    int mountainheight(int x, int y) const;
    float mountainheightfloat(int x, int y) const;
    int craterrim(int x, int y) const;
    float craterrimfloat(int x, int y) const;
    int volcano(int x, int y) const;
    float volcanofloat(int x, int y) const;
    int lakesurface(int x, int y) const;
    float lakesurfacefloat(int x, int y) const;
    int riftlakesurface(int x, int y) const;
    float riftlakesurfacefloat(int x, int y) const;
    int special(int x, int y) const;
    int truelake(int x, int y) const;

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
    void setjantemp(int x, int y, int amount);
    int jultemp(int x, int y) const;
    void setjultemp(int x, int y, int amount);
    float jantempfloat(int x, int y) const;
    void setjantempfloat(int x, int y, float amount);
    float jultempfloat(int x, int y) const;
    void setjultempfloat(int x, int y, float amount);
    int seasonaltemp(int season, int x, int y) const;
    void setseasonaltemp(int season, int x, int y, int amount);
    float seasonaltempfloat(int season, int x, int y) const;
    void setseasonaltempfloat(int season, int x, int y, float amount);
    int janrain(int x, int y) const;
    void setjanrain(int x, int y, int amount);
    int julrain(int x, int y) const;
    void setjulrain(int x, int y, int amount);
    float janrainfloat(int x, int y) const;
    void setjanrainfloat(int x, int y, float amount);
    float julrainfloat(int x, int y) const;
    void setjulrainfloat(int x, int y, float amount);
    int seasonalrain(int season, int x, int y) const;
    void setseasonalrain(int season, int x, int y, int amount);
    float seasonalrainfloat(int season, int x, int y) const;
    void setseasonalrainfloat(int season, int x, int y, float amount);
    int seasonalpressure(int season, int x, int y) const;
    void setseasonalpressure(int season, int x, int y, int amount);
    float seasonalpressurefloat(int season, int x, int y) const;
    void setseasonalpressurefloat(int season, int x, int y, float amount);
    int seasonaluwind(int season, int x, int y) const;
    void setseasonaluwind(int season, int x, int y, int amount);
    float seasonaluwindfloat(int season, int x, int y) const;
    void setseasonaluwindfloat(int season, int x, int y, float amount);
    int seasonalvwind(int season, int x, int y) const;
    void setseasonalvwind(int season, int x, int y, int amount);
    float seasonalvwindfloat(int season, int x, int y) const;
    void setseasonalvwindfloat(int season, int x, int y, float amount);
    int seasonalupperheight(int season, int x, int y) const;
    void setseasonalupperheight(int season, int x, int y, int amount);
    float seasonalupperheightfloat(int season, int x, int y) const;
    void setseasonalupperheightfloat(int season, int x, int y, float amount);
    int seasonalupperuwind(int season, int x, int y) const;
    void setseasonalupperuwind(int season, int x, int y, int amount);
    float seasonalupperuwindfloat(int season, int x, int y) const;
    void setseasonalupperuwindfloat(int season, int x, int y, float amount);
    int seasonaluppervwind(int season, int x, int y) const;
    void setseasonaluppervwind(int season, int x, int y, int amount);
    float seasonaluppervwindfloat(int season, int x, int y) const;
    void setseasonaluppervwindfloat(int season, int x, int y, float amount);
    int seasonalverticalvelocity(int season, int x, int y) const;
    void setseasonalverticalvelocity(int season, int x, int y, int amount);
    float seasonalverticalvelocityfloat(int season, int x, int y) const;
    void setseasonalverticalvelocityfloat(int season, int x, int y, float amount);
    int seasonalcurrentu(int season, int x, int y) const;
    void setseasonalcurrentu(int season, int x, int y, int amount);
    float seasonalcurrentufloat(int season, int x, int y) const;
    void setseasonalcurrentufloat(int season, int x, int y, float amount);
    int seasonalcurrentv(int season, int x, int y) const;
    void setseasonalcurrentv(int season, int x, int y, int amount);
    float seasonalcurrentvfloat(int season, int x, int y) const;
    void setseasonalcurrentvfloat(int season, int x, int y, float amount);
    int seasonalsst(int season, int x, int y) const;
    void setseasonalsst(int season, int x, int y, int amount);
    float seasonalsstfloat(int season, int x, int y) const;
    void setseasonalsstfloat(int season, int x, int y, float amount);
    int seasonalevaporation(int season, int x, int y) const;
    void setseasonalevaporation(int season, int x, int y, int amount);
    float seasonalevaporationfloat(int season, int x, int y) const;
    void setseasonalevaporationfloat(int season, int x, int y, float amount);
    int seasonalmaritimeinfluence(int season, int x, int y) const;
    void setseasonalmaritimeinfluence(int season, int x, int y, int amount);
    float seasonalmaritimeinfluencefloat(int season, int x, int y) const;
    void setseasonalmaritimeinfluencefloat(int season, int x, int y, float amount);
    int seasonalmaritimethermalanomaly(int season, int x, int y) const;
    void setseasonalmaritimethermalanomaly(int season, int x, int y, int amount);
    float seasonalmaritimethermalanomalyfloat(int season, int x, int y) const;
    void setseasonalmaritimethermalanomalyfloat(int season, int x, int y, float amount);
    int seasonalmaritimefetch(int season, int x, int y) const;
    void setseasonalmaritimefetch(int season, int x, int y, int amount);
    int seasonalmoisture(int season, int x, int y) const;
    void setseasonalmoisture(int season, int x, int y, int amount);
    float seasonalmoisturefloat(int season, int x, int y) const;
    void setseasonalmoisturefloat(int season, int x, int y, float amount);
    int seasonalconvergence(int season, int x, int y) const;
    void setseasonalconvergence(int season, int x, int y, int amount);
    float seasonalconvergencefloat(int season, int x, int y) const;
    void setseasonalconvergencefloat(int season, int x, int y, float amount);
    int seasonaluplift(int season, int x, int y) const;
    void setseasonaluplift(int season, int x, int y, int amount);
    float seasonalupliftfloat(int season, int x, int y) const;
    void setseasonalupliftfloat(int season, int x, int y, float amount);
    int seasonalsubsidence(int season, int x, int y) const;
    void setseasonalsubsidence(int season, int x, int y, int amount);
    float seasonalsubsidencefloat(int season, int x, int y) const;
    void setseasonalsubsidencefloat(int season, int x, int y, float amount);
    int janmountainrain(int x, int y) const;
    void setjanmountainrain(int x, int y, int amount);
    int julmountainrain(int x, int y) const;
    void setjulmountainrain(int x, int y, int amount);
    int wintermountainrain(int x, int y) const;
    void setwintermountainrain(int x, int y, int amount);
    int summermountainrain(int x, int y) const;
    void setsummermountainrain(int x, int y, int amount);
    void setwintermountainrainfloat(int x, int y, float amount);
    void setsummermountainrainfloat(int x, int y, float amount);
    int janmountainraindir(int x, int y) const;
    void setjanmountainraindir(int x, int y, int amount);
    int julmountainraindir(int x, int y) const;
    void setjulmountainraindir(int x, int y, int amount);
    int wintermountainraindir(int x, int y) const;
    int summermountainraindir(int x, int y) const;
    void setwintermountainraindir(int x, int y, int amount);
    void setsummermountainraindir(int x, int y, int amount);
    int climate(int x, int y) const;
    void setclimate(int x, int y, int amount);
    int biome(int x, int y) const;
    void setbiome(int x, int y, int amount);
    int seaice(int x, int y) const;
    void setseaice(int x, int y, int amount);
    int wind(int x, int y) const;
    void setwind(int x, int y, int amount);
    int horse(int x, int y) const;
    void sethorse(int x, int y, int amount);
    int mintemp(int x, int y) const;
    void setmintemp(int x, int y, int amount);
    int maxtemp(int x, int y) const;
    void setmaxtemp(int x, int y, int amount);
    int winterrain(int x, int y) const;
    void setwinterrain(int x, int y, int amount);
    int summerrain(int x, int y) const;
    void setsummerrain(int x, int y, int amount);
    // The same selections on the stored floats.
    float mintempfloat(int x, int y) const;
    void setmintempfloat(int x, int y, float amount);
    float maxtempfloat(int x, int y) const;
    void setmaxtempfloat(int x, int y, float amount);
    float winterrainfloat(int x, int y) const;
    void setwinterrainfloat(int x, int y, float amount);
    float summerrainfloat(int x, int y) const;
    void setsummerrainfloat(int x, int y, float amount);

    // Annual and transition-season views of the seasonal fields.
    int avetemp(int x, int y) const;
    float avetempfloat(int x, int y) const;
    int averain(int x, int y) const;
    float averainfloat(int x, int y) const;
    int aprtemp(int x, int y) const;
    int octtemp(int x, int y) const;
    int aprrain(int x, int y) const;
    int octrain(int x, int y) const;

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
    std::size_t inputindex(int x, int y) const { return static_cast<std::size_t>(x) * inputs_.rows + y; }
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

inline int climateworld::nom(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    return inside(x, y) ? wholemetres(inputs_.nom[inputindex(x, y)]) : 0;
}

inline float climateworld::nomfloat(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    return inside(x, y) ? inputs_.nom[inputindex(x, y)] : 0.0f;
}

// The base elevation plus what map() adds to it; the land/sea test stays on whole metres.
inline float climateworld::mapfloat(int x, int y) const
{
    if (!inside(x, y))
        return 0.0f;
    UW_FIELD_READ(volcano);
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(ocean_ridge_height);
    UW_FIELD_READ(mountain_height);
    UW_FIELD_READ(extra_elevation);
    UW_FIELD_READ(crater_rims);
    UW_FIELD_READ(crater_centres);
    const std::size_t cell = inputindex(x, y);
    const float base = inputs_.nom[cell];
    const float thisvolcano = std::fabs(inputs_.volcano[cell]);

    if (wholemetres(base) <= inputs_.sealevel)
        return base + inputs_.oceanridgeheight[cell] + thisvolcano;

    return base + (std::max)(thisvolcano, inputs_.mountainheights[cell]) + inputs_.extraelev[cell] + inputs_.craterrims[cell] + inputs_.cratercentres[cell];
}

// Whole metres: the base and what is added to it are rounded separately, as planet::map does.
inline int climateworld::map(int x, int y) const
{
    if (!inside(x, y))
        return 0;

    const float base = inputs_.nom[inputindex(x, y)];
    return wholemetres(base) + wholemetres(mapfloat(x, y) - base);
}

inline bool climateworld::sea(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    UW_FIELD_READ(volcano);
    if (!inside(x, y))
        return 0;

    const std::size_t cell = inputindex(x, y);
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

inline bool climateworld::basebelowsea(int x, int y) const
{
    UW_FIELD_READ(sea_level);
    return nom(x, y) <= inputs_.sealevel;
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

inline int climateworld::mountainheight(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    return inside(x, y) ? wholemetres(inputs_.mountainheights[inputindex(x, y)]) : 0;
}

inline float climateworld::mountainheightfloat(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    return inside(x, y) ? inputs_.mountainheights[inputindex(x, y)] : 0.0f;
}

inline int climateworld::craterrim(int x, int y) const
{
    UW_FIELD_READ(crater_rims);
    return inside(x, y) ? wholemetres(inputs_.craterrims[inputindex(x, y)]) : 0;
}

inline float climateworld::craterrimfloat(int x, int y) const
{
    UW_FIELD_READ(crater_rims);
    return inside(x, y) ? inputs_.craterrims[inputindex(x, y)] : 0.0f;
}

inline int climateworld::volcano(int x, int y) const
{
    UW_FIELD_READ(volcano);
    return inside(x, y) ? wholemetres(inputs_.volcano[inputindex(x, y)]) : 0;
}

inline float climateworld::volcanofloat(int x, int y) const
{
    UW_FIELD_READ(volcano);
    return inside(x, y) ? inputs_.volcano[inputindex(x, y)] : 0.0f;
}

inline int climateworld::lakesurface(int x, int y) const
{
    UW_FIELD_READ(lake_surface);
    return inside(x, y) ? wholemetres(inputs_.lakesurface[inputindex(x, y)]) : 0;
}

inline float climateworld::lakesurfacefloat(int x, int y) const
{
    UW_FIELD_READ(lake_surface);
    return inside(x, y) ? inputs_.lakesurface[inputindex(x, y)] : 0.0f;
}

inline int climateworld::riftlakesurface(int x, int y) const
{
    UW_FIELD_READ(rift_lake_surface);
    return inside(x, y) ? wholemetres(inputs_.riftlakesurface[inputindex(x, y)]) : 0;
}

inline float climateworld::riftlakesurfacefloat(int x, int y) const
{
    UW_FIELD_READ(rift_lake_surface);
    return inside(x, y) ? inputs_.riftlakesurface[inputindex(x, y)] : 0.0f;
}

inline int climateworld::special(int x, int y) const
{
    UW_FIELD_READ(special_features);
    return inside(x, y) ? static_cast<int>(inputs_.special[inputindex(x, y)]) : 0;
}

inline int climateworld::truelake(int x, int y) const
{
    if (!inside(x, y))
        return 0;

    if (lakesurface(x, y) != 0 && special(x, y) < 110)
        return (1);
    else
        return(0);
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

inline const short* climateworld::rawseaicemap() const { UW_FIELD_READ(sea_ice); return &seaicemap[0][0]; }

inline const short* climateworld::rawclimatemap() const { UW_FIELD_READ(climate_class); return &climatemap[0][0]; }

inline const short* climateworld::rawbiomemap() const { UW_FIELD_READ(biome); return &biomemap[0][0]; }

inline int climateworld::jantemp(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(jantempmap[x][y]);
}

inline void climateworld::setjantemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    jantempmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::jultemp(int x, int y) const
{
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(jultempmap[x][y]);
}

inline void climateworld::setjultemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    jultempmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::seasonaltemp(int season, int x, int y) const
{
    UW_FIELD_READ(temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonaltempmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonaltemp(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonaltempmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::janrain(int x, int y) const
{
    UW_FIELD_READ(precipitation_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(janrainmap[x][y]);
}

inline void climateworld::setjanrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    janrainmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::julrain(int x, int y) const
{
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(julrainmap[x][y]);
}

inline void climateworld::setjulrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    julrainmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::seasonalrain(int season, int x, int y) const
{
    UW_FIELD_READ(precipitation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalrainmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalrain(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    const int index = seasonalclimateindex(x, y);
    seasonalrainmaps[season][index] = static_cast<float>(amount);
}

inline float climateworld::seasonalrainfloat(int season, int x, int y) const
{
    UW_FIELD_READ(precipitation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalrainmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalrainfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    const int index = seasonalclimateindex(x, y);
    const float physicalamount = (std::max)(0.0f, amount);
    seasonalrainmaps[season][index] = physicalamount;
}

inline int climateworld::seasonalpressure(int season, int x, int y) const
{
    UW_FIELD_READ(pressure);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalpressuremaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalpressure(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalpressuremaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonaluwind(int season, int x, int y) const
{
    UW_FIELD_READ(wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonaluwindmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonaluwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonaluwindmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalvwind(int season, int x, int y) const
{
    UW_FIELD_READ(wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalvwindmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalvwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalvwindmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalupperheight(int season, int x, int y) const
{
    UW_FIELD_READ(upper_height);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalupperheightmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalupperheight(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalupperheightmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalupperuwind(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalupperuwindmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalupperuwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalupperuwindmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonaluppervwind(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonaluppervwindmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonaluppervwind(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonaluppervwindmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalverticalvelocity(int season, int x, int y) const
{
    UW_FIELD_READ(vertical_velocity);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalverticalvelocitymaps[season][seasonalclimateindex(x, y)], climatelegacyscale::verticalvelocity);
}

inline void climateworld::setseasonalverticalvelocity(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalverticalvelocitymaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::verticalvelocity;
}

inline int climateworld::seasonalcurrentu(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalcurrentumaps[season][seasonalclimateindex(x, y)], climatelegacyscale::current);
}

inline void climateworld::setseasonalcurrentu(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalcurrentumaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::current;
}

inline int climateworld::seasonalcurrentv(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalcurrentvmaps[season][seasonalclimateindex(x, y)], climatelegacyscale::current);
}

inline void climateworld::setseasonalcurrentv(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalcurrentvmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::current;
}

inline int climateworld::seasonalsst(int season, int x, int y) const
{
    UW_FIELD_READ(sea_surface_temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalsstmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalsst(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalsstmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalevaporation(int season, int x, int y) const
{
    UW_FIELD_READ(evaporation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalevaporationmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalevaporation(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalevaporationmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalmaritimeinfluence(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_influence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalmaritimeinfluencemaps[season][seasonalclimateindex(x, y)], climatelegacyscale::maritimeinfluence);
}

inline void climateworld::setseasonalmaritimeinfluence(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmaritimeinfluencemaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::maritimeinfluence;
}

inline int climateworld::seasonalmaritimethermalanomaly(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_thermal_anomaly);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalmaritimethermalanomalymaps[season][seasonalclimateindex(x, y)], climatelegacyscale::maritimethermalanomaly);
}

inline void climateworld::setseasonalmaritimethermalanomaly(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmaritimethermalanomalymaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::maritimethermalanomaly;
}

inline int climateworld::seasonalmaritimefetch(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_fetch);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return static_cast<int>(seasonalmaritimefetchmaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalmaritimefetch(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmaritimefetchmaps[season][seasonalclimateindex(x, y)] = static_cast<short>(amount);
}

inline int climateworld::seasonalmoisture(int season, int x, int y) const
{
    UW_FIELD_READ(moisture);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalmoisturemaps[season][seasonalclimateindex(x, y)]);
}

inline void climateworld::setseasonalmoisture(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmoisturemaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount);
}

inline int climateworld::seasonalconvergence(int season, int x, int y) const
{
    UW_FIELD_READ(convergence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalconvergencemaps[season][seasonalclimateindex(x, y)], climatelegacyscale::convergence);
}

inline void climateworld::setseasonalconvergence(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalconvergencemaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::convergence;
}

inline int climateworld::seasonaluplift(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_uplift);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalupliftmaps[season][seasonalclimateindex(x, y)], climatelegacyscale::orographicmotion);
}

inline void climateworld::setseasonaluplift(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalupliftmaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::orographicmotion;
}

inline int climateworld::seasonalsubsidence(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_subsidence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(seasonalsubsidencemaps[season][seasonalclimateindex(x, y)], climatelegacyscale::orographicmotion);
}

inline void climateworld::setseasonalsubsidence(int season, int x, int y, int amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalsubsidencemaps[season][seasonalclimateindex(x, y)] = static_cast<float>(amount) / climatelegacyscale::orographicmotion;
}

// Float accessors: the stored values in physical units.

inline float climateworld::seasonaltempfloat(int season, int x, int y) const
{
    UW_FIELD_READ(temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonaltempmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonaltempfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonaltempmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalpressurefloat(int season, int x, int y) const
{
    UW_FIELD_READ(pressure);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalpressuremaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalpressurefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalpressuremaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonaluwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonaluwindmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonaluwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonaluwindmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalvwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalvwindmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalvwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalvwindmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalupperheightfloat(int season, int x, int y) const
{
    UW_FIELD_READ(upper_height);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalupperheightmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalupperheightfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalupperheightmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalupperuwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalupperuwindmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalupperuwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalupperuwindmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonaluppervwindfloat(int season, int x, int y) const
{
    UW_FIELD_READ(upper_wind_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonaluppervwindmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonaluppervwindfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonaluppervwindmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalverticalvelocityfloat(int season, int x, int y) const
{
    UW_FIELD_READ(vertical_velocity);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalverticalvelocitymaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalverticalvelocityfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalverticalvelocitymaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalcurrentufloat(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_east);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalcurrentumaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalcurrentufloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalcurrentumaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalcurrentvfloat(int season, int x, int y) const
{
    UW_FIELD_READ(ocean_current_south);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalcurrentvmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalcurrentvfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalcurrentvmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalsstfloat(int season, int x, int y) const
{
    UW_FIELD_READ(sea_surface_temperature);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalsstmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalsstfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalsstmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalevaporationfloat(int season, int x, int y) const
{
    UW_FIELD_READ(evaporation);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalevaporationmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalevaporationfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalevaporationmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalmaritimeinfluencefloat(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_influence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalmaritimeinfluencemaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalmaritimeinfluencefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmaritimeinfluencemaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalmaritimethermalanomalyfloat(int season, int x, int y) const
{
    UW_FIELD_READ(maritime_thermal_anomaly);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalmaritimethermalanomalymaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalmaritimethermalanomalyfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmaritimethermalanomalymaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalmoisturefloat(int season, int x, int y) const
{
    UW_FIELD_READ(moisture);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalmoisturemaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalmoisturefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalmoisturemaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalconvergencefloat(int season, int x, int y) const
{
    UW_FIELD_READ(convergence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalconvergencemaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalconvergencefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalconvergencemaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalupliftfloat(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_uplift);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalupliftmaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalupliftfloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalupliftmaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::seasonalsubsidencefloat(int season, int x, int y) const
{
    UW_FIELD_READ(orographic_subsidence);
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return seasonalsubsidencemaps[season][seasonalclimateindex(x, y)];
}

inline void climateworld::setseasonalsubsidencefloat(int season, int x, int y, float amount)
{
    if (!validseasonindex(season) || y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seasonalsubsidencemaps[season][seasonalclimateindex(x, y)] = amount;
}

inline float climateworld::jantempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jantempmap[x][y];
}

inline void climateworld::setjantempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    jantempmap[x][y] = amount;
}

inline float climateworld::jultempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jultempmap[x][y];
}

inline void climateworld::setjultempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    jultempmap[x][y] = amount;
}

inline float climateworld::janrainfloat(int x, int y) const
{
    UW_FIELD_READ(precipitation_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return janrainmap[x][y];
}

inline void climateworld::setjanrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    janrainmap[x][y] = amount;
}

inline float climateworld::julrainfloat(int x, int y) const
{
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return julrainmap[x][y];
}

inline void climateworld::setjulrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    julrainmap[x][y] = amount;
}

inline int climateworld::janmountainrain(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(janmountainrainmap[x][y]);
}

inline void climateworld::setjanmountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    janmountainrainmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::julmountainrain(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return legacyvalue(julmountainrainmap[x][y]);
}

inline void climateworld::setjulmountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    julmountainrainmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::wintermountainrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_january);
    UW_FIELD_READ(mountain_precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        return legacyvalue(janmountainrainmap[x][y]);
    else
        return legacyvalue(julmountainrainmap[x][y]);
}

inline void climateworld::setwintermountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        janmountainrainmap[x][y] = static_cast<float>(amount);
    else
        julmountainrainmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::summermountainrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_january);
    UW_FIELD_READ(mountain_precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) >= legacyvalue(jultempmap[x][y]))
        return legacyvalue(janmountainrainmap[x][y]);
    else
        return legacyvalue(julmountainrainmap[x][y]);
}

inline void climateworld::setsummermountainrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) >= legacyvalue(jultempmap[x][y]))
        janmountainrainmap[x][y] = static_cast<float>(amount);
    else
        julmountainrainmap[x][y] = static_cast<float>(amount);
}

inline void climateworld::setwintermountainrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    (jantempmap[x][y] < jultempmap[x][y] ? janmountainrainmap[x][y] : julmountainrainmap[x][y]) = amount;
}

inline void climateworld::setsummermountainrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    (jantempmap[x][y] >= jultempmap[x][y] ? janmountainrainmap[x][y] : julmountainrainmap[x][y]) = amount;
}

inline int climateworld::janmountainraindir(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_direction_january);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return (int)janmountainraindirmap[x][y];
}

inline void climateworld::setjanmountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    janmountainraindirmap[x][y] = (short)amount;
}

inline int climateworld::julmountainraindir(int x, int y) const
{
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return (int)julmountainraindirmap[x][y];
}

inline void climateworld::setjulmountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    julmountainraindirmap[x][y] = (short)amount;
}

inline int climateworld::wintermountainraindir(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_direction_january);
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        return (int)janmountainraindirmap[x][y];
    else
        return (int)julmountainraindirmap[x][y];
}

inline int climateworld::summermountainraindir(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(mountain_precipitation_direction_january);
    UW_FIELD_READ(mountain_precipitation_direction_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) >= legacyvalue(jultempmap[x][y]))
        return (int)janmountainraindirmap[x][y];
    else
        return (int)julmountainraindirmap[x][y];
}

inline void climateworld::setwintermountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        janmountainraindirmap[x][y] = (short)amount;
    else
        julmountainraindirmap[x][y] = (short)amount;
}

inline void climateworld::setsummermountainraindir(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) >= legacyvalue(jultempmap[x][y]))
        janmountainraindirmap[x][y] = (short)amount;
    else
        julmountainraindirmap[x][y] = (short)amount;
}

inline int climateworld::climate(int x, int y) const
{
    UW_FIELD_READ(climate_class);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return (int)climatemap[x][y];
}

inline void climateworld::setclimate(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    climatemap[x][y] = (short)amount;
}

inline int climateworld::biome(int x, int y) const
{
    UW_FIELD_READ(biome);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return (int)biomemap[x][y];
}

inline void climateworld::setbiome(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    biomemap[x][y] = (short)amount;
}

inline int climateworld::seaice(int x, int y) const
{
    UW_FIELD_READ(sea_ice);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return (int)seaicemap[x][y];
}

inline void climateworld::setseaice(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    seaicemap[x][y] = (short)amount;
}

inline int climateworld::wind(int x, int y) const
{
    UW_FIELD_READ(prevailing_wind);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    return windmap[x][y];
}

inline void climateworld::setwind(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    windmap[x][y] = amount;
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

inline int climateworld::mintemp(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        return legacyvalue(jantempmap[x][y]);
    else
        return legacyvalue(jultempmap[x][y]);
}

inline void climateworld::setmintemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        jantempmap[x][y] = static_cast<float>(amount);
    else
        jultempmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::maxtemp(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) > legacyvalue(jultempmap[x][y]))
        return legacyvalue(jantempmap[x][y]);
    else
        return legacyvalue(jultempmap[x][y]);
}

inline void climateworld::setmaxtemp(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) > legacyvalue(jultempmap[x][y]))
        jantempmap[x][y] = static_cast<float>(amount);
    else
        jultempmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::winterrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        return legacyvalue(janrainmap[x][y]);
    else
        return legacyvalue(julrainmap[x][y]);
}

inline void climateworld::setwinterrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) < legacyvalue(jultempmap[x][y]))
        janrainmap[x][y] = static_cast<float>(amount);
    else
        julrainmap[x][y] = static_cast<float>(amount);
}

inline int climateworld::summerrain(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    if (legacyvalue(jantempmap[x][y]) >= legacyvalue(jultempmap[x][y]))
        return legacyvalue(janrainmap[x][y]);
    else
        return legacyvalue(julrainmap[x][y]);
}

inline void climateworld::setsummerrain(int x, int y, int amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    if (legacyvalue(jantempmap[x][y]) >= legacyvalue(jultempmap[x][y]))
        janrainmap[x][y] = static_cast<float>(amount);
    else
        julrainmap[x][y] = static_cast<float>(amount);
}

inline float climateworld::mintempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jantempmap[x][y] < jultempmap[x][y] ? jantempmap[x][y] : jultempmap[x][y];
}

inline void climateworld::setmintempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    (jantempmap[x][y] < jultempmap[x][y] ? jantempmap[x][y] : jultempmap[x][y]) = amount;
}

inline float climateworld::maxtempfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jantempmap[x][y] > jultempmap[x][y] ? jantempmap[x][y] : jultempmap[x][y];
}

inline void climateworld::setmaxtempfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    (jantempmap[x][y] > jultempmap[x][y] ? jantempmap[x][y] : jultempmap[x][y]) = amount;
}

inline float climateworld::winterrainfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jantempmap[x][y] < jultempmap[x][y] ? janrainmap[x][y] : julrainmap[x][y];
}

inline void climateworld::setwinterrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    (jantempmap[x][y] < jultempmap[x][y] ? janrainmap[x][y] : julrainmap[x][y]) = amount;
}

inline float climateworld::summerrainfloat(int x, int y) const
{
    UW_FIELD_READ(temperature_january);
    UW_FIELD_READ(temperature_july);
    UW_FIELD_READ(precipitation_january);
    UW_FIELD_READ(precipitation_july);
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    return jantempmap[x][y] >= jultempmap[x][y] ? janrainmap[x][y] : julrainmap[x][y];
}

inline void climateworld::setsummerrainfloat(int x, int y, float amount)
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return;

    (jantempmap[x][y] >= jultempmap[x][y] ? janrainmap[x][y] : julrainmap[x][y]) = amount;
}

inline float climateworld::avetempfloat(int x, int y) const
{
    float total = 0.0f;
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += seasonaltempfloat(season, x, y);
    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

inline int climateworld::avetemp(int x, int y) const
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    int total = 0;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total = total + seasonaltemp(season, x, y);

    return total / CLIMATESEASONCOUNT;
}

inline int climateworld::averain(int x, int y) const
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0;

    int total = 0;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total = total + seasonalrain(season, x, y);

    return total / CLIMATESEASONCOUNT;
}

inline float climateworld::averainfloat(int x, int y) const
{
    if (y<0 || y>inputs_.height || x<0 || x>inputs_.width)
        return 0.0f;

    float total = 0.0f;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        total += seasonalrainfloat(season, x, y);

    return total / static_cast<float>(CLIMATESEASONCOUNT);
}

inline int climateworld::aprtemp(int x, int y) const
{
    return seasonaltemp(seasonapril, x, y);
}

inline int climateworld::octtemp(int x, int y) const
{
    return seasonaltemp(seasonoctober, x, y);
}

inline int climateworld::aprrain(int x, int y) const
{
    return seasonalrain(seasonapril, x, y);
}

inline int climateworld::octrain(int x, int y) const
{
    return seasonalrain(seasonoctober, x, y);
}
