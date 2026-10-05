#include "simulations/climate/climate_world_fields.hpp"

#include <algorithm>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>

#include "simulations/climate/climate_world.hpp"
#include "world/generation_tuning.hpp"

using worldmodel::FieldData;
using worldmodel::FieldId;
using worldmodel::TableId;

namespace
{
// An input grid seen through its x-major pointer, indexed view[x][y].
template<typename T>
struct inputview
{
    const T* values;
    int stride;
    const T* operator[](int x) const { return values + static_cast<std::size_t>(x) * stride; }
};
}

// Friend of climateworld: maps catalogue ids onto the climate state.
struct climateworldfieldaccess
{
    // Climate state stored on the world grid, for either direction (field_encoding.hpp).
    template<typename World, typename Operation>
    static bool state(World& world, FieldId id, Operation& op)
    {
        switch (id)
        {
        case FieldId::north_polar_temperature: return op.scalar(id, world.itsnorthpolartemp);
        case FieldId::south_polar_temperature: return op.scalar(id, world.itssouthpolartemp);
        case FieldId::equatorial_temperature: return op.scalar(id, world.itseqtemp);
        case FieldId::sea_ice: return op.grid(id, world.seaicemap);
        case FieldId::temperature_january: return op.grid(id, world.jantempmap);
        case FieldId::temperature_july: return op.grid(id, world.jultempmap);
        case FieldId::precipitation_january: return op.grid(id, world.janrainmap);
        case FieldId::precipitation_july: return op.grid(id, world.julrainmap);
        case FieldId::mountain_precipitation_january: return op.grid(id, world.janmountainrainmap);
        case FieldId::mountain_precipitation_july: return op.grid(id, world.julmountainrainmap);
        case FieldId::mountain_precipitation_direction_january: return op.grid(id, world.janmountainraindirmap);
        case FieldId::mountain_precipitation_direction_july: return op.grid(id, world.julmountainraindirmap);
        case FieldId::prevailing_wind: return op.grid(id, world.windmap);
        case FieldId::climate_class: return op.grid(id, world.climatemap);
        case FieldId::biome: return op.grid(id, world.biomemap);
        case FieldId::subtropical_high_latitudes: return op.zonalprofile(id, world.horselats);
        case FieldId::temperature: return op.seasonal(id, world.seasonaltempmaps);
        case FieldId::precipitation: return op.seasonal(id, world.seasonalrainmaps);
        case FieldId::precipitation_exact: return op.seasonal(id, world.seasonalrainfloatmaps);
        case FieldId::pressure: return op.seasonal(id, world.seasonalpressuremaps);
        case FieldId::wind_east: return op.seasonal(id, world.seasonaluwindmaps);
        case FieldId::wind_south: return op.seasonal(id, world.seasonalvwindmaps);
        case FieldId::upper_height: return op.seasonal(id, world.seasonalupperheightmaps);
        case FieldId::upper_wind_east: return op.seasonal(id, world.seasonalupperuwindmaps);
        case FieldId::upper_wind_south: return op.seasonal(id, world.seasonaluppervwindmaps);
        case FieldId::vertical_velocity: return op.seasonal(id, world.seasonalverticalvelocitymaps);
        case FieldId::evaporation: return op.seasonal(id, world.seasonalevaporationmaps);
        case FieldId::maritime_influence: return op.seasonal(id, world.seasonalmaritimeinfluencemaps);
        case FieldId::maritime_thermal_anomaly: return op.seasonal(id, world.seasonalmaritimethermalanomalymaps);
        case FieldId::maritime_fetch: return op.seasonal(id, world.seasonalmaritimefetchmaps);
        case FieldId::moisture: return op.seasonal(id, world.seasonalmoisturemaps);
        case FieldId::convergence: return op.seasonal(id, world.seasonalconvergencemaps);
        case FieldId::orographic_uplift: return op.seasonal(id, world.seasonalupliftmaps);
        case FieldId::orographic_subsidence: return op.seasonal(id, world.seasonalsubsidencemaps);
        case FieldId::ocean_current_east: return op.seasonal(id, world.seasonalcurrentumaps);
        case FieldId::ocean_current_south: return op.seasonal(id, world.seasonalcurrentvmaps);
        case FieldId::sea_surface_temperature: return op.seasonal(id, world.seasonalsstmaps);
        default: return false;
        }
    }

    static bool read(const climateworld& world, FieldId id, FieldData& data, bool includeinputs)
    {
        namespace encode = worldmodel::encoding;
        const std::size_t columns = static_cast<std::size_t>(world.inputs_.width) + 1;
        const std::size_t rows = static_cast<std::size_t>(world.inputs_.height) + 1;
        const auto& context = world.climatecontext_;

        encode::ReadOperation op;
        op.columns = columns;
        op.rows = rows;
        if (state(world, id, op))
        {
            data = std::move(op.data);
            return true;
        }

        switch (id)
        {
        case FieldId::solver_monthly_temperature:
            data = encode::monthly(id, static_cast<std::size_t>(context.monthly.columns), static_cast<std::size_t>(context.monthly.rows), context.monthly.temperatureC);
            return true;
        case FieldId::solver_monthly_precipitation:
            data = encode::monthly(id, static_cast<std::size_t>(context.monthly.columns), static_cast<std::size_t>(context.monthly.rows), context.monthly.rainfallMm);
            return true;
        case FieldId::solver_surface_wind_east: data = encode::solverseasonal(id, context.winds.surfacewindu); return true;
        case FieldId::solver_surface_wind_south: data = encode::solverseasonal(id, context.winds.surfacewindv); return true;
        case FieldId::solver_upper_wind_east: data = encode::solverseasonal(id, context.winds.upperwindu); return true;
        case FieldId::solver_upper_wind_south: data = encode::solverseasonal(id, context.winds.upperwindv); return true;
        case FieldId::solver_low_wind_east: data = encode::solverseasonal(id, context.winds.lowlevelwindu); return true;
        case FieldId::solver_low_wind_south: data = encode::solverseasonal(id, context.winds.lowlevelwindv); return true;
        case FieldId::solver_sea_surface_temperature: data = encode::solverseasonal(id, context.oceanSst); return true;
        case FieldId::solver_sea_ice: data = encode::solverseasonal(id, context.oceanIce); return true;
        case FieldId::level_wind_east: data = levelfield(id, world, columns, rows, &atmospherelevel::eastwind); return true;
        case FieldId::level_wind_north: data = levelfield(id, world, columns, rows, &atmospherelevel::northwind); return true;
        case FieldId::level_geopotential_height: data = levelfield(id, world, columns, rows, &atmospherelevel::geopotentialheight); return true;
        case FieldId::level_vertical_velocity: data = levelfield(id, world, columns, rows, &atmospherelevel::verticalvelocity); return true;
        case FieldId::level_available: data = levelfield(id, world, columns, rows, &atmospherelevel::available); return true;
        case FieldId::column_water:
        case FieldId::moisture_transport_east:
        case FieldId::moisture_transport_north:
        case FieldId::moisture_convergence:
        case FieldId::precipitation_rate:
        case FieldId::ascent:
            return processfield(id, world, data);
        case FieldId::current_east:
        case FieldId::current_north:
        case FieldId::mixed_layer_temperature:
        case FieldId::storage_depth:
            return oceanfield(id, world, data);
        default: break;
        }

        if (!includeinputs)
            return false;

        const climateinputs& inputs = world.inputs_;
        auto inputgrid = [&](const auto* values)
        {
            using T = std::remove_cv_t<std::remove_pointer_t<decltype(values)>>;
            return encode::grid(id, columns, rows, inputview<T>{ values, inputs.rows });
        };

        switch (id)
        {
        case FieldId::seed: data = encode::scalar(id, inputs.seed); return true;
        case FieldId::width_index: data = encode::scalar(id, inputs.width); return true;
        case FieldId::height_index: data = encode::scalar(id, inputs.height); return true;
        case FieldId::prograde_rotation: data = encode::scalar(id, inputs.rotation); return true;
        case FieldId::axial_tilt: data = encode::scalar(id, inputs.tilt); return true;
        case FieldId::eccentricity: data = encode::scalar(id, inputs.eccentricity); return true;
        case FieldId::perihelion: data = encode::scalar(id, inputs.perihelion); return true;
        case FieldId::gravity: data = encode::scalar(id, inputs.gravity); return true;
        case FieldId::lunar_strength: data = encode::scalar(id, inputs.lunar); return true;
        case FieldId::lapse_rate: data = encode::scalar(id, inputs.tempdecrease); return true;
        case FieldId::north_polar_adjust: data = encode::scalar(id, inputs.northpolaradjust); return true;
        case FieldId::south_polar_adjust: data = encode::scalar(id, inputs.southpolaradjust); return true;
        case FieldId::average_temperature: data = encode::scalar(id, inputs.averagetemp); return true;
        case FieldId::water_pickup: data = encode::scalar(id, inputs.waterpickup); return true;
        case FieldId::max_elevation: data = encode::scalar(id, inputs.maxelevation); return true;
        case FieldId::sea_level: data = encode::scalar(id, inputs.sealevel); return true;
        case FieldId::land_cells: data = encode::scalar(id, inputs.landtotal); return true;
        case FieldId::sea_cells: data = encode::scalar(id, inputs.seatotal); return true;
        case FieldId::base_elevation: if (!inputs.nom) return false; data = inputgrid(inputs.nom); return true;
        case FieldId::mountain_height: if (!inputs.mountainheights) return false; data = inputgrid(inputs.mountainheights); return true;
        case FieldId::extra_elevation: if (!inputs.extraelev) return false; data = inputgrid(inputs.extraelev); return true;
        case FieldId::crater_rims: if (!inputs.craterrims) return false; data = inputgrid(inputs.craterrims); return true;
        case FieldId::crater_centres: if (!inputs.cratercentres) return false; data = inputgrid(inputs.cratercentres); return true;
        case FieldId::volcano: if (!inputs.volcano) return false; data = inputgrid(inputs.volcano); return true;
        case FieldId::ocean_ridge_height: if (!inputs.oceanridgeheight) return false; data = inputgrid(inputs.oceanridgeheight); return true;
        case FieldId::lake_surface: if (!inputs.lakesurface) return false; data = inputgrid(inputs.lakesurface); return true;
        case FieldId::rift_lake_surface: if (!inputs.riftlakesurface) return false; data = inputgrid(inputs.riftlakesurface); return true;
        case FieldId::special_features: if (!inputs.special) return false; data = inputgrid(inputs.special); return true;
        default: return false;
        }
    }

    // seasons x rows x columns on a solver's own grid; seasons without data are NaN.
    template<typename Value>
    static FieldData nativeseasonal(FieldId id, int columns, int rows, const Value& value)
    {
        worldmodel::encoding::checkelement<float>(id);
        const std::size_t cells = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
        FieldData data;
        data.id = id;
        data.shape = { static_cast<std::size_t>(CLIMATESEASONCOUNT), static_cast<std::size_t>(rows),
            static_cast<std::size_t>(columns) };
        data.bytes.resize(CLIMATESEASONCOUNT * cells * sizeof(float));
        auto* out = reinterpret_cast<float*>(data.bytes.data());
        std::fill(out, out + CLIMATESEASONCOUNT * cells, std::numeric_limits<float>::quiet_NaN());
        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
            value(season, out + season * cells);
        return data;
    }

    // Diagnostics of the moisture processes, as the circulation diagnostics export them.
    static bool processfield(FieldId id, const climateworld& world, FieldData& data)
    {
        const climatediagnostics* diagnostics = world.diagnostics();
        if (diagnostics == nullptr || !diagnostics->hasprocess)
            return false;
        int columns = 0, rows = 0;
        for (const auto& fields : diagnostics->process)
            if (fields.durationSeconds > 0.0 && fields.columns > 0)
            {
                columns = fields.columns;
                rows = fields.rows;
            }
        if (columns == 0)
            return false;
        const std::size_t cells = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
        data = nativeseasonal(id, columns, rows, [&](int season, float* out)
        {
            const auto& fields = diagnostics->process[season];
            if (fields.durationSeconds <= 0.0 || fields.columns != columns || fields.rows != rows)
                return;
            const auto copy = [&](const std::vector<float>& values, float sign)
            {
                if (values.size() == cells)
                    for (std::size_t cell = 0; cell < cells; cell++)
                        out[cell] = sign * values[cell];
            };
            if (id == FieldId::column_water)
                copy(fields.columnWaterMm, 1.0f);
            else if (id == FieldId::ascent)
                copy(fields.ascentHpaPerDay, 1.0f);
            else if (id == FieldId::precipitation_rate)
            {
                for (const auto& rain : fields.surfacePrecipitationMm)
                    if (rain.size() != cells)
                        return;
                for (std::size_t cell = 0; cell < cells; cell++)
                {
                    double total = 0.0;
                    for (const auto& rain : fields.surfacePrecipitationMm)
                        total += rain[cell] * 86400.0 / fields.durationSeconds;
                    out[cell] = static_cast<float>(total);
                }
            }
            else
            {
                const auto transport = climatehydrology::meanMoistureTransport(
                    fields, tuning::climate::atmosphere::referencePlanetRadiusMetres);
                if (id == FieldId::moisture_convergence)
                    copy(transport.convergenceMmPerDay, 1.0f);
                else
                {
                    // Lower plus upper layer; the solver's southward component becomes northward.
                    const auto& layers = id == FieldId::moisture_transport_east ?
                        transport.eastKgPerMetreSecond : transport.southKgPerMetreSecond;
                    const float sign = id == FieldId::moisture_transport_east ? 1.0f : -1.0f;
                    if (layers[0].size() == cells && layers[1].size() == cells)
                        for (std::size_t cell = 0; cell < cells; cell++)
                            out[cell] = sign * (layers[0][cell] + layers[1][cell]);
                }
            }
        });
        return true;
    }

    // Ocean state of the last solve on the ocean grid.
    static bool oceanfield(FieldId id, const climateworld& world, FieldData& data)
    {
        const climatediagnostics* diagnostics = world.diagnostics();
        if (diagnostics == nullptr || !diagnostics->hasocean || diagnostics->oceancolumns <= 0)
            return false;
        const int columns = diagnostics->oceancolumns, rows = diagnostics->oceanrows;
        const std::size_t cells = static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows);
        data = nativeseasonal(id, columns, rows, [&](int season, float* out)
        {
            const auto& state = diagnostics->ocean[season];
            const std::vector<float>* values = nullptr;
            float sign = 1.0f;
            if (id == FieldId::current_east)
                values = &state.eastCurrentMps;
            else if (id == FieldId::current_north)
            {
                values = &state.southCurrentMps;
                sign = -1.0f;
            }
            else if (id == FieldId::mixed_layer_temperature)
                values = state.meanSstC.size() == cells ? &state.meanSstC : &state.sstC;
            else
                values = &state.meanStorageDepthMetres;
            if (values->size() == cells)
                for (std::size_t cell = 0; cell < cells; cell++)
                    out[cell] = sign * (*values)[cell];
        });
        return true;
    }

    // levels x seasons x rows x columns; levels or seasons without data are zero.
    template<typename T>
    static FieldData levelfield(FieldId id, const climateworld& world, std::size_t columns, std::size_t rows,
        std::array<std::vector<T>, CLIMATESEASONCOUNT> atmospherelevel::* member)
    {
        worldmodel::encoding::checkelement<T>(id);
        const std::size_t cells = columns * rows;
        const std::size_t levels = world.levels_.size();
        FieldData data;
        data.id = id;
        data.shape = { levels, static_cast<std::size_t>(CLIMATESEASONCOUNT), rows, columns };
        data.bytes.assign(levels * CLIMATESEASONCOUNT * cells * sizeof(T), 0);
        auto* out = reinterpret_cast<T*>(data.bytes.data());
        for (std::size_t level = 0; level < levels; level++)
            for (std::size_t season = 0; season < CLIMATESEASONCOUNT; season++)
            {
                const std::vector<T>& values = (world.levels_[level].*member)[season];
                if (values.size() != cells)
                    continue;
                std::copy(values.begin(), values.end(), out + (level * CLIMATESEASONCOUNT + season) * cells);
            }
        return data;
    }

    static bool table(const climateworld& world, TableId id, std::string& text, std::string& extension)
    {
        if (id == TableId::atmosphere_levels)
        {
            std::ostringstream out;
            out << "index,name,pressure_hpa\n";
            for (std::size_t index = 0; index < world.levels_.size(); index++)
                out << index << ',' << world.levels_[index].name << ',' << world.levels_[index].pressurehpa << '\n';
            text = out.str();
            extension = ".csv";
            return true;
        }
        if (id != TableId::weather_state)
            return false;
        std::ostringstream out;
        out << std::hex << std::setfill('0');
        for (std::uint8_t byte : world.weatheranomalystatebytes)
            out << std::setw(2) << static_cast<int>(byte);
        text = out.str();
        extension = ".hex";
        return true;
    }
};

bool readclimateworldfield(const climateworld& world, FieldId id, FieldData& data, bool includeinputs)
{
    return climateworldfieldaccess::read(world, id, data, includeinputs);
}

bool loadclimateworldfield(climateworld& world, const FieldData& data)
{
    worldmodel::encoding::LoadOperation op;
    op.columns = static_cast<std::size_t>(world.width()) + 1;
    op.rows = static_cast<std::size_t>(world.height()) + 1;
    op.data = &data;
    return climateworldfieldaccess::state(world, data.id, op);
}

bool readclimateworldtable(const climateworld& world, TableId id, std::string& text, std::string& extension)
{
    return climateworldfieldaccess::table(world, id, text, extension);
}

bool isclimatediagnosticfield(FieldId id)
{
    switch (id)
    {
    case FieldId::column_water:
    case FieldId::moisture_transport_east:
    case FieldId::moisture_transport_north:
    case FieldId::moisture_convergence:
    case FieldId::precipitation_rate:
    case FieldId::ascent:
    case FieldId::current_east:
    case FieldId::current_north:
    case FieldId::mixed_layer_temperature:
    case FieldId::storage_depth:
        return true;
    default:
        return false;
    }
}

std::vector<FieldId> climateworldfields(bool includeinputs)
{
    std::vector<FieldId> fields =
    {
        FieldId::north_polar_temperature, FieldId::south_polar_temperature, FieldId::equatorial_temperature,
        FieldId::sea_ice, FieldId::temperature_january, FieldId::temperature_july, FieldId::precipitation_january,
        FieldId::precipitation_july, FieldId::mountain_precipitation_january, FieldId::mountain_precipitation_july,
        FieldId::mountain_precipitation_direction_january, FieldId::mountain_precipitation_direction_july,
        FieldId::prevailing_wind, FieldId::climate_class, FieldId::biome, FieldId::subtropical_high_latitudes,
        FieldId::temperature, FieldId::precipitation, FieldId::precipitation_exact, FieldId::pressure, FieldId::wind_east,
        FieldId::wind_south, FieldId::upper_height, FieldId::upper_wind_east, FieldId::upper_wind_south,
        FieldId::vertical_velocity, FieldId::evaporation, FieldId::maritime_influence, FieldId::maritime_thermal_anomaly,
        FieldId::maritime_fetch, FieldId::moisture, FieldId::convergence, FieldId::orographic_uplift,
        FieldId::orographic_subsidence, FieldId::ocean_current_east, FieldId::ocean_current_south,
        FieldId::sea_surface_temperature, FieldId::solver_monthly_temperature, FieldId::solver_monthly_precipitation,
        FieldId::solver_surface_wind_east, FieldId::solver_surface_wind_south, FieldId::solver_upper_wind_east,
        FieldId::solver_upper_wind_south, FieldId::solver_low_wind_east, FieldId::solver_low_wind_south,
        FieldId::solver_sea_surface_temperature, FieldId::solver_sea_ice, FieldId::level_wind_east, FieldId::level_wind_north,
        FieldId::level_geopotential_height, FieldId::level_vertical_velocity, FieldId::level_available,
        FieldId::column_water, FieldId::moisture_transport_east, FieldId::moisture_transport_north,
        FieldId::moisture_convergence, FieldId::precipitation_rate, FieldId::ascent, FieldId::current_east,
        FieldId::current_north, FieldId::mixed_layer_temperature, FieldId::storage_depth
    };
    if (includeinputs)
    {
        const FieldId inputs[] =
        {
            FieldId::seed, FieldId::width_index, FieldId::height_index, FieldId::prograde_rotation, FieldId::axial_tilt,
            FieldId::eccentricity, FieldId::perihelion, FieldId::gravity, FieldId::lunar_strength, FieldId::lapse_rate,
            FieldId::north_polar_adjust, FieldId::south_polar_adjust, FieldId::average_temperature, FieldId::water_pickup,
            FieldId::max_elevation, FieldId::sea_level, FieldId::land_cells, FieldId::sea_cells, FieldId::base_elevation,
            FieldId::mountain_height, FieldId::extra_elevation, FieldId::crater_rims, FieldId::crater_centres,
            FieldId::volcano, FieldId::ocean_ridge_height, FieldId::lake_surface, FieldId::rift_lake_surface,
            FieldId::special_features
        };
        fields.insert(fields.end(), std::begin(inputs), std::end(inputs));
    }
    std::sort(fields.begin(), fields.end());
    return fields;
}
