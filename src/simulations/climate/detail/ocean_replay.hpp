#pragma once
#include "climate_ocean_dynamics.hpp"
#include <atomic>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace climateocean::detail
{
// Diagnostic-only, versioned text fixture. No struct padding or save-format
// dependency; max_digits10 roundtrips every float forcing/config value.
struct OceanReplay
{
    int columns = 0, rows = 0, season = 0;
    OceanConfig config;
    OceanForcing forcing;

    template<class C, class F> static void configFields(C& c, F&& io)
    {
        io(c.planetRadiusMetres, c.rotationRatePerSecond, c.rotationDirection, c.airDensityKgM3,
            c.waterDensityKgM3, c.dragCoefficient, c.barotropicDragPerSecond, c.linearBottomDragMps,
            c.minimumCoriolisPerSecond, c.mixedLayerDepthMetres, c.heatDiffusivityM2S,
            c.surfaceHeatExchangeWm2K, c.waterHeatCapacityJkgK, c.freezingTemperatureC, c.iceDensityKgM3,
            c.latentHeatFusionJkg, c.iceConductivityWmK, c.oceanTimeStepSeconds, c.streamfunctionIterations,
            c.heatStepsPerIteration, c.couplingIterations, c.underRelaxation, c.convergenceTolerance,
            c.sstWindFeedbackMpsPerK, c.deepWaterTemperatureContrastK, c.streamfunctionTolerance,
            c.maximumCurrentMps, c.oneWay, c.implicitZonalHeatDiffusion, c.reuseSeasonalCirculation,
            c.krylovCirculation, c.parallelZonalHeatDiffusion, c.gatherHeatFluxes, c.heatWorkers);
    }
    template<class F, class IO> static void forcingFields(F& f, IO&& io)
    {
        io(f.landMask); io(f.bathymetryMetres); io(f.eastWindMps); io(f.southWindMps);
        io(f.atmosphericTemperatureC); io(f.initialSstC); io(f.initialIceThicknessMetres);
        io(f.surfaceHeatFluxWm2); io(f.surfaceHeatFluxReferenceTemperatureC); io(f.deepWaterTemperatureC);
    }
    template<class C, class F> static void storageConfigFields(C& c, F&& io)
    {
        io(c.variableHeatStorage, c.minimumStorageDepthMetres, c.maximumStorageDepthMetres,
            c.storageColumnDepthMetres, c.storageStratificationPerSecond2, c.storageWindMixingEfficiency,
            c.storageMixingMemoryDays, c.storageAdjustmentDays);
    }
    std::filesystem::path write(const std::string& directory) const
    {
        static std::atomic<unsigned> sequence{0};
        std::filesystem::create_directories(directory);
        const auto path = std::filesystem::path(directory) / ("ocean-" + std::to_string(sequence++) + ".txt");
        if (std::filesystem::exists(path)) throw std::runtime_error("Ocean capture already exists: " + path.string());
        std::ofstream out(path);
        const bool storage = config.variableHeatStorage || !forcing.initialStorageDepthMetres.empty() ||
            !forcing.initialReservoirTemperatureC.empty();
        out << std::setprecision(std::numeric_limits<float>::max_digits10)
            << (storage ? "UW_OCEAN_V2 " : "UW_OCEAN_V1 ") << columns << ' ' << rows << ' ' << season << '\n';
        configFields(config, [&](const auto&... v) { ((out << v << ' '), ...); });
        out << static_cast<int>(config.coastalScheme) << '\n';
        if (storage)
        {
            storageConfigFields(config, [&](const auto&... v) { ((out << v << ' '), ...); });
            out << '\n';
        }
        const auto writeField = [&](const auto& field)
        {
            out << field.size() << '\n';
            for (const auto value : field) out << +value << '\n';
        };
        forcingFields(forcing, writeField);
        if (storage)
        {
            writeField(forcing.initialStorageDepthMetres);
            writeField(forcing.initialReservoirTemperatureC);
        }
        if (!out) throw std::runtime_error("Ocean capture write failed: " + path.string());
        return path;
    }
    static OceanReplay read(const std::string& path)
    {
        OceanReplay f;
        std::ifstream in(path);
        std::string version;
        in >> version >> f.columns >> f.rows >> f.season;
        if (!in || (version != "UW_OCEAN_V1" && version != "UW_OCEAN_V2") || f.columns < 8 || f.columns > 2048 ||
            f.columns % 4 != 0 || f.rows != f.columns / 2 || f.season < 0 || f.season > 3)
            throw std::runtime_error("Invalid ocean replay header");
        configFields(f.config, [&](auto&... v) { ((in >> v), ...); });
        int scheme = -1;
        in >> scheme;
        if (!in || scheme < 0 || scheme > 2) throw std::runtime_error("Invalid ocean replay config");
        f.config.coastalScheme = static_cast<CoastalScheme>(scheme);
        if (version == "UW_OCEAN_V2") storageConfigFields(f.config, [&](auto&... v) { ((in >> v), ...); });
        const auto readField = [&](auto& field)
        {
            std::size_t size = 0;
            in >> size;
            if (!in || (size != 0 && size != static_cast<std::size_t>(f.columns) * f.rows))
                throw std::runtime_error("Invalid ocean replay field size");
            field.resize(size);
            for (auto& value : field)
            {
                if constexpr (std::is_same_v<std::decay_t<decltype(value)>, std::uint8_t>)
                {
                    int mask = -1; in >> mask;
                    if (mask != 0 && mask != 1) throw std::runtime_error("Invalid ocean replay mask");
                    value = static_cast<std::uint8_t>(mask);
                }
                else in >> value;
            }
        };
        forcingFields(f.forcing, readField);
        if (version == "UW_OCEAN_V2")
        {
            readField(f.forcing.initialStorageDepthMetres);
            readField(f.forcing.initialReservoirTemperatureC);
        }
        const auto count = static_cast<std::size_t>(f.columns) * f.rows;
        if (!in || f.forcing.landMask.size() != count || f.forcing.bathymetryMetres.size() != count ||
            f.forcing.eastWindMps.size() != count || f.forcing.southWindMps.size() != count ||
            f.forcing.atmosphericTemperatureC.size() != count || f.forcing.initialSstC.size() != count)
            throw std::runtime_error("Truncated ocean replay");
        in >> std::ws;
        if (!in.eof()) throw std::runtime_error("Trailing ocean replay data");
        return f;
    }
};
}
