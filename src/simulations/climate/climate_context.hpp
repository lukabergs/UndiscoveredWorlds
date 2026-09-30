#pragma once

#include "climate_hydrology.hpp"
#include "climate_atmosphere.hpp"
#include "climate_ocean_dynamics.hpp"
#include <algorithm>
#include <array>
#include <vector>

namespace climatecontext
{
enum class TransportMethod { Ffsl, Mpdata };
using FloatGrid = std::vector<std::vector<float>>;

template<typename T>
class RowMajorGrid
{
public:
    template<typename Value> struct Column
    {
        Value* values;
        int columns, x;
        Value& operator[](int y) const { return values[static_cast<std::size_t>(y) * columns + x]; }
    };
    RowMajorGrid() = default;
    RowMajorGrid(int columns, int rows, T value = {})
        : columns_(columns), values_(static_cast<std::size_t>(columns) * rows, value) {}
    Column<T> operator[](int x) { return {values_.data(), columns_, x}; }
    Column<const T> operator[](int x) const { return {values_.data(), columns_, x}; }
    void fill(T value) { std::fill(values_.begin(), values_.end(), value); }
    std::vector<T>& values() { return values_; }
    const std::vector<T>& values() const { return values_; }
private:
    int columns_ = 0;
    std::vector<T> values_;
};

struct Resolution
{
    int atmosphereColumns = 64;
    int oceanColumns = 64;
    int hydrologyColumns = 128;
    int weatherColumns = 32;
};

struct WindFields
{
    int width = -1, height = -1;
    std::array<FloatGrid, 4> surfacewindu, surfacewindv, upperwindu, upperwindv;
    std::array<FloatGrid, 4> lowlevelwindu, lowlevelwindv, lowlevelavailable;
    std::array<bool, 4> hasLowLevelProfile{};
    // Surface-minus-carrier drainage increment, in the same east/south basis.
    std::array<FloatGrid, 4> drainagewindu, drainagewindv;
    std::array<bool, 4> hasDrainageProfile{};
    std::array<bool, 4> populated{};
    void reset(int w, int h)
    {
        width = w; height = h; populated.fill(false); hasDrainageProfile.fill(false);
        hasLowLevelProfile.fill(false);
    }
    bool completeFor(int w, int h) const
    {
        return width == w && height == h &&
            std::all_of(populated.begin(), populated.end(), [](bool value) { return value; });
    }
};

struct HydrologyStorage
{
    int columns = 0, rows = 0;
    RowMajorGrid<float> boundary, free, soil, snow;
    RowMajorGrid<double> landSurfaceC, landDeepC;
};

struct OceanStorage
{
    int columns = 0, rows = 0;
    std::vector<float> sstC, iceMetres, deepTemperatureC;
    std::vector<float> storageDepthMetres, reservoirTemperatureC;
    double annualEnthalpyDriftK = 0.0;
    bool periodic = false;
};

struct MonthlyClimate
{
    int columns = 0, rows = 0;
    std::array<std::vector<float>, 12> temperatureC, rainfallMm;
    // Reference seasonal values allow later terrain/import adjustments to be
    // applied without losing the generated monthly shape.
    std::array<std::vector<float>, 4> referenceTemperatureC, referenceRainfallMm;
    template<std::size_t N> bool complete(const std::array<std::vector<float>, N>& fields) const
    {
        const auto cells = static_cast<std::size_t>(columns) * rows;
        return columns > 0 && rows > 0 && std::all_of(fields.begin(), fields.end(),
            [cells](const auto& field) { return field.size() == cells; });
    }
    bool hasTemperature() const { return complete(temperatureC) && complete(referenceTemperatureC); }
    bool hasRainfall() const { return complete(rainfallMm) && complete(referenceRainfallMm); }
};

// Owned by one world and copied by value with workbench snapshots. No process-
// global world pointers, save-format fields, or presentation callbacks.
struct Context
{
    Resolution resolution;
    bool calibrateMeanTemperature = true;
    // Baseline generation solves seasonal circulation; transient weather is
    // an explicit comparison mode, independent of the moisture integrator.
    bool transientWeather = false;
    // Experimental until the coupled rainfall/upper-flow regression is resolved.
    bool mixedLayerTropics = false;
    TransportMethod transportMethod = TransportMethod::Ffsl;
    int oceanHeatStepHours = 24;
    bool reuseOceanCirculation = true;
    bool krylovOceanCirculation = true;
    bool gatherOceanHeatFluxes = true;
    int oceanHeatWorkers = 0;
    climateocean::CoastalScheme oceanCoastalScheme = climateocean::CoastalScheme::FaceForcing;
    std::string oceanCaptureDirectory;
    climateatmosphere::StationarySolver atmosphereSolver = climateatmosphere::StationarySolver::Zonal;
    std::string atmosphereCaptureDirectory;
    std::array<std::vector<float>, 4> radiativeZonalTemperature;
    std::array<climateatmosphere::BoundaryLayerState, 4> boundaryLayer;
    int couplingIterationLimit = 0; // Zero uses the tuning default.
    MonthlyClimate monthly;
    WindFields winds;
    HydrologyStorage storage;
    OceanStorage oceanStorage;
    int hydrologyYears = 0;
    double landAnnualDriftK = 0.0, landEnergyResidualJm2 = 0.0;
    std::array<climatehydrology::SeasonalProcessFields, 4> processes;
    std::array<FloatGrid, 4> oceanSst, oceanSkin, oceanIce;
    std::array<bool, 4> oceanAccepted{}, atmosphereAccepted{};
    bool hasOcean = false, hasProcesses = false;

    void resetDynamics()
    {
        winds = {}; storage = {}; processes = {}; oceanStorage = {};
        hydrologyYears = 0;
        radiativeZonalTemperature = {};
        boundaryLayer = {};
        landAnnualDriftK = landEnergyResidualJm2 = 0.0;
        oceanSst = {}; oceanSkin = {}; oceanIce = {};
        oceanAccepted.fill(false); atmosphereAccepted.fill(false);
        hasOcean = false; hasProcesses = false;
    }
};
}
