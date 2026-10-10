#pragma once

// Catalogue of world state. One X-macro list defines every field's id, owning
// component, element type, layout, units and meaning; FieldId, bundle names and
// declared process access all derive from it. Units describe the current
// storage: several climate fields keep a legacy integer encoding until they
// migrate to float32 physical units (docs/architecture/world-model.md).

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace worldmodel
{
enum class Component : std::uint8_t
{
    planetbody,
    geosphere,
    hydrosphere,
    cryosphere,
    atmosphere,
    biosphere,
    anthroposphere,
    presentation
};

enum class ElementType : std::uint8_t { b8, u8, i16, u16, i32, f32 };

// grid: W x H surface cells. seasonal/monthly: 4 or 12 grids. zonalprofile:
// six values per column. levelseasonal: levels x 4 x grid (atmosphere column).
// scalar: one value. Solver grids may differ from the
// surface grid; their shape is taken from the stored data.
enum class Layout : std::uint8_t { scalar, grid, seasonal, monthly, zonalprofile, levelseasonal };

// X(id, component, element, layout, units, description)
#define UW_WORLD_FIELD_LIST(X) \
    X(seed, planetbody, i32, scalar, "", "World seed") \
    X(size_class, planetbody, i32, scalar, "", "Planet size class: 1 tiny, 2 small, 3 Earthlike") \
    X(radius, planetbody, f32, scalar, "m", "Planet radius") \
    X(width_index, planetbody, i32, scalar, "cells", "Largest x index of the surface grid (width - 1)") \
    X(height_index, planetbody, i32, scalar, "cells", "Largest y index of the surface grid (height - 1)") \
    X(prograde_rotation, planetbody, b8, scalar, "", "True when the planet rotates like Earth") \
    X(axial_tilt, planetbody, f32, scalar, "degrees", "Axial tilt") \
    X(eccentricity, planetbody, f32, scalar, "", "Orbital eccentricity") \
    X(perihelion, planetbody, i32, scalar, "month", "Month of perihelion, 0 = January") \
    X(gravity, planetbody, f32, scalar, "g", "Surface gravity relative to Earth") \
    X(lunar_strength, planetbody, f32, scalar, "", "Lunar tidal strength relative to Earth") \
    X(lapse_rate, atmosphere, f32, scalar, "C/km", "Temperature decrease with elevation") \
    X(north_polar_adjust, atmosphere, i32, scalar, "C", "Adjustment to north polar temperature") \
    X(south_polar_adjust, atmosphere, i32, scalar, "C", "Adjustment to south polar temperature") \
    X(average_temperature, atmosphere, i32, scalar, "C", "Average world temperature") \
    X(north_polar_temperature, atmosphere, i32, scalar, "C", "Average north polar temperature") \
    X(south_polar_temperature, atmosphere, i32, scalar, "C", "Average south polar temperature") \
    X(equatorial_temperature, atmosphere, i32, scalar, "C", "Average equatorial temperature") \
    X(water_pickup, atmosphere, f32, scalar, "", "Ocean moisture pickup factor") \
    X(river_factor, hydrosphere, f32, scalar, "", "Divisor converting legacy river flow to m3/s") \
    X(river_land_reduce, hydrosphere, i32, scalar, "m", "Depth of river channels below the cell base") \
    X(estuary_limit, hydrosphere, i32, scalar, "legacy flow", "Minimum river flow for estuaries") \
    X(glacial_temperature, cryosphere, i32, scalar, "C", "Maximum mean temperature for glacially carved terrain") \
    X(glacier_temperature, cryosphere, i32, scalar, "C", "Maximum temperature for rivers to become glaciers") \
    X(mountain_reduce, geosphere, f32, scalar, "", "Global/regional mountain height matching factor") \
    X(climate_type_count, atmosphere, i32, scalar, "", "Number of climate types") \
    X(max_elevation, geosphere, i32, scalar, "m", "Largest elevation value of the datum") \
    X(sea_level, geosphere, i32, scalar, "m above datum", "Sea level on the elevation datum") \
    X(max_river_flow, hydrosphere, i32, scalar, "legacy flow", "Largest river flow") \
    X(land_cells, geosphere, i32, scalar, "cells", "Number of land cells") \
    X(sea_cells, geosphere, i32, scalar, "cells", "Number of sea cells") \
    X(crater_count, geosphere, i32, scalar, "", "Number of craters") \
    X(tectonic_time_origin_step, geosphere, i32, scalar, "updates", "Tectonic time origin") \
    X(tectonic_time, geosphere, f32, scalar, "Myr", "Tectonic time of the final state") \
    X(tectonic_time_step, geosphere, f32, scalar, "Myr", "Tectonic update length") \
    X(tectonic_plate_count, geosphere, i32, scalar, "", "Number of plates in the final state") \
    X(tectonic_sea_level, geosphere, i32, scalar, "m", "Engine sea level on its uint16 datum") \
    X(base_elevation, geosphere, f32, grid, "m above datum", "Terrain without mountains, volcanoes and craters (nom)") \
    X(mountain_height, geosphere, f32, grid, "m", "Mountain elevation added to the base") \
    X(mountain_ridges, geosphere, i16, grid, "direction bits", "Mountain ridge directions") \
    X(extra_elevation, geosphere, f32, grid, "m", "Extra terrain elevation (canyon uplift)") \
    X(crater_rims, geosphere, f32, grid, "m", "Crater rim elevation") \
    X(crater_centres, geosphere, f32, grid, "m", "Crater centre elevation") \
    X(volcano, geosphere, f32, grid, "m", "Isolated peaks; negative values are extinct volcanoes") \
    X(continental_shelf, geosphere, b8, grid, "", "Continental shelf: shallow sea margins of land") \
    X(stratovolcano, geosphere, b8, grid, "", "Stratovolcano flag") \
    X(island, geosphere, b8, grid, "", "Island flag") \
    X(mountain_island, geosphere, b8, grid, "", "Mountain island flag") \
    X(ocean_ridge, geosphere, i16, grid, "code", "Ocean ridge mask") \
    X(ocean_ridge_height, geosphere, f32, grid, "m", "Ocean ridge elevation") \
    X(ocean_rift, geosphere, i16, grid, "code", "Ocean rift mask") \
    X(ocean_ridge_offset, geosphere, i16, grid, "cells", "Ocean ridge offset") \
    X(ocean_ridge_angle, geosphere, i16, grid, "degrees", "Ocean ridge angle") \
    X(roughness, geosphere, f32, grid, "noise", "Seeded fractal roughness used by landforms and the regional map") \
    X(detail_noise, geosphere, i16, grid, "noise", "Terracing noise for the regional map") \
    X(geologic_regime, geosphere, u8, grid, "code", "Geologic regime from tectonic signals") \
    X(tectonic_convergence, geosphere, u8, grid, "0..100", "Boundary convergence score") \
    X(tectonic_divergence, geosphere, u8, grid, "0..100", "Boundary divergence score") \
    X(tectonic_shear, geosphere, u8, grid, "0..100", "Boundary shear score") \
    X(crust_age, geosphere, f32, grid, "Myr", "Crust age") \
    X(crust_thickness, geosphere, f32, grid, "engine units", "Crust thickness proxy") \
    X(crust_class, geosphere, u8, grid, "code", "Crust class: 0 none, 1 oceanic, 2 transitional, 3 continental") \
    X(uplift_tendency, geosphere, f32, grid, "engine units", "Native uplift tendency") \
    X(rock_uplift_rate, geosphere, f32, grid, "m/Myr", "Tectonic rock uplift rate, 10 Myr exponential mean of the engine's elevation change before surface processes, in the plate's frame (semi-Lagrangian: the elevation before each step and the mean are carried with the rock by the engine's own transport, exact for a translation where one column field varies, a bounded approximation at plate boundaries; the transport smears the mean by a fraction of a cell; converging rock mixes, divergent cells keep the rock that stayed, new crust has no history of its own); negative is subsidence. The mean starts at 0 with the run, so a run of t Myr shorter than the window underestimates it by the factor 1 - e^(-t/10)") \
    X(shortening, geosphere, f32, grid, "m", "Cumulative horizontal shortening of the rock over the tectonic run, in the plate's frame: minus the most compressive principal component of the engine's cumulative strain tensor (the sum of each update's strain sym(grad u) of the applied displacement, within a plate, the tensor rotated by the local spin and carried with the rock by the engine's own transport), times the native cell width (world_width_km / engine columns): strain x the native cell width, so it depends on the engine's resolution (divide by the native cell width for the strain) and is a small-strain measure (the strain reaches about 1.4 at 256 columns); 0 where no axis shortens. Jumps of velocity between plates are not strain; new ocean floor starts at 0. Nearest engine cell") \
    X(shortening_direction, geosphere, f32, grid, "rad", "Axis of the cumulative shortening, radians in [0, pi) from the grid's x axis (east) toward its y axis (south); an axis, not a vector (pi equals 0). Meaningful where geosphere.shortening is positive; 0 where the principal components are equal. Nearest engine cell, never interpolated") \
    X(subsidence_tendency, geosphere, f32, grid, "engine units", "Native subsidence tendency") \
    X(accumulated_strain, geosphere, f32, grid, "engine units", "Accumulated strain memory") \
    X(boundary_type, geosphere, u8, grid, "code", "Plate boundary type") \
    X(boundary_distance, geosphere, u16, grid, "cells", "Distance to the nearest boundary") \
    X(boundary_segment, geosphere, i32, grid, "id", "Boundary segment id") \
    X(nearest_boundary, geosphere, i32, grid, "id", "Nearest boundary segment id") \
    X(boundary_history, geosphere, f32, grid, "engine units", "Boundary lifecycle history") \
    X(deforming_region, geosphere, i32, grid, "id", "Deforming region id") \
    X(deforming_region_type, geosphere, u8, grid, "code", "Deforming region type") \
    X(deformation_rate, geosphere, f32, grid, "engine units", "Deformation rate") \
    X(deformation_velocity_x, geosphere, f32, grid, "world cells/update", "Deformation velocity, eastward") \
    X(deformation_velocity_y, geosphere, f32, grid, "world cells/update", "Deformation velocity, southward") \
    X(tectonic_height, geosphere, f32, grid, "engine units", "Engine height of the last tectonic update: 0.1 at the ocean floor datum, 1 at sea level on continental crust") \
    X(plate_id, geosphere, i32, grid, "id", "Plate index of the last tectonic update") \
    X(crust_age_steps, geosphere, i32, grid, "updates", "Crust age in tectonic updates") \
    X(continental_crust, geosphere, f32, grid, "m", "Continental crust of the isostatic column") \
    X(oceanic_crust, geosphere, f32, grid, "m", "Oceanic crust of the isostatic column") \
    X(sediment, geosphere, f32, grid, "m", "Sediment of the isostatic column") \
    X(rift_damage, geosphere, f32, grid, "0..1", "Advected rift damage") \
    X(plume_influence, geosphere, f32, grid, "0..1", "Proximity to an active mantle plume; 0 without plumes") \
    X(erosion_potential, geosphere, u8, grid, "0..100", "Relative erosion potential") \
    X(deposition_potential, geosphere, u8, grid, "0..100", "Relative deposition potential") \
    X(floodplain_fertility, geosphere, u8, grid, "0..100", "Floodplain and lowland fertility proxy") \
    X(erodibility, geosphere, f32, grid, "relative", "Erodibility relative to landscape.stream_power.erodibility") \
    X(erosion_rate, geosphere, f32, grid, "m/Myr", "Mean erosion rate over the erosion run") \
    X(deposition_rate, geosphere, f32, grid, "m/Myr", "Mean deposition rate over the erosion run") \
    X(denudation, geosphere, f32, grid, "m", "Cumulative denudation of the erosion run") \
    X(aggradation, geosphere, f32, grid, "m", "Cumulative aggradation of the erosion run") \
    X(glacial_denudation, geosphere, f32, grid, "m", "Glacial part of the denudation of the erosion run") \
    X(prograded, geosphere, u8, grid, "0/1", "Land built by marine deposition over the erosion run: sea at its start, land at its end") \
    X(relief_amplitude, geosphere, f32, grid, "m", "Grid-scale relief of a cell (terrain.subgrid_statistics): the steady-state relief E l^2 / (2 D) of the erosion run's rate E and grid-scale hillslope diffusivity D over half the cell size l, at most the cell's local relief; 0 at sea") \
    X(valley_spacing, geosphere, f32, grid, "km", "First-order valley spacing of Perron et al. 2009 (terrain.subgrid_statistics): c (D/K)^(1/(2 m_ref + 1)) m, D/K calibrated to a Table 1 site (terrain.subgrid.reference_*) and scaled per cell by K_A (erodibility x runoff^m_ref) relative to the reference, c terrain.subgrid.spacing_coefficient; a value equal to the cell size (sqrt of its area) means no sub-grid valleys resolved; 0 at sea and where no runoff or erodibility") \
    X(erosion_duration, geosphere, f32, scalar, "Myr", "Duration of the erosion run") \
    X(eroded_volume, geosphere, f32, scalar, "km3", "Volume eroded over the erosion run") \
    X(land_deposited_volume, geosphere, f32, scalar, "km3", "Volume deposited on land over the erosion run") \
    X(marine_deposited_volume, geosphere, f32, scalar, "km3", "Volume deposited in the sea over the erosion run") \
    X(metal_ore_potential, geosphere, u8, grid, "0..100", "Metal ore potential score") \
    X(placer_potential, geosphere, u8, grid, "0..100", "Placer potential score") \
    X(evaporite_potential, geosphere, u8, grid, "0..100", "Evaporite potential score") \
    X(volcanic_potential, geosphere, u8, grid, "0..100", "Volcanic resource potential score") \
    X(river_direction, hydrosphere, i16, grid, "direction code", "River flow direction") \
    X(river_flow_january, hydrosphere, i32, grid, "legacy flow", "January river flow") \
    X(river_flow_july, hydrosphere, i32, grid, "legacy flow", "July river flow") \
    X(lake_surface, hydrosphere, f32, grid, "m above datum", "Lake surface elevation") \
    X(rift_lake_surface, hydrosphere, f32, grid, "m above datum", "Rift lake surface elevation") \
    X(flow_direction, hydrosphere, u8, grid, "direction code", "Flow direction of the routing (1..8 as river_direction; 0 sea)") \
    X(drainage_area, hydrosphere, f32, grid, "km2", "Drainage area of the routing") \
    X(discharge, hydrosphere, f32, grid, "m3/s", "Annual mean discharge of the routing") \
    X(filled_elevation, hydrosphere, f32, grid, "m above datum", "Priority-Flood surface: equals the terrain outside depressions") \
    X(depression_id, hydrosphere, i32, grid, "id", "Depression of the routing (0 none)") \
    X(runoff, hydrosphere, f32, grid, "mm/yr", "Runoff depth the routing accumulates (land and sea cells)") \
    X(sediment_flux, hydrosphere, f32, grid, "m3/yr", "Sediment a land cell passes to its receiver (an outlet: into its sea cell), mean over the erosion run; 0 at sea") \
    X(imported_lake, hydrosphere, u8, grid, "0/1", "Lake water of an imported lake raster (Earth imports); 0 in procedural worlds") \
    X(rift_lake_bed, hydrosphere, f32, grid, "m above datum", "Rift lake bed elevation") \
    X(lake_start, hydrosphere, b8, grid, "", "Lake seed cell") \
    X(special_features, hydrosphere, i16, grid, "code", "100 salt lake, 110 salt pan, 120 dunes, 130/131/132 fresh/brackish/salt wetlands") \
    X(tide, hydrosphere, i16, grid, "legacy", "Tidal strength") \
    X(basin_class, hydrosphere, u8, grid, "code", "Drainage basin class") \
    X(delta_direction, hydrosphere, i16, grid, "direction code", "Delta channel direction") \
    X(delta_flow_january, hydrosphere, i32, grid, "legacy flow", "January delta flow") \
    X(delta_flow_july, hydrosphere, i32, grid, "legacy flow", "July delta flow") \
    X(sea_ice, cryosphere, i16, grid, "code", "Sea ice: 1 seasonal, 2 permanent") \
    X(equilibrium_line_altitude, cryosphere, f32, grid, "m above datum", "Glacier equilibrium-line altitude of the erosion run (its climate cooled by landscape.glacial.cooling_c)") \
    X(glaciation, cryosphere, f32, grid, "0..1", "Share of the erosion run's time a cell carried ice") \
    X(ice_thickness, cryosphere, f32, grid, "m", "Present-day ice on the final terrain (the ELA without the cooling)") \
    X(temperature_january, atmosphere, f32, grid, "C", "January temperature (legacy pair)") \
    X(temperature_july, atmosphere, f32, grid, "C", "July temperature (legacy pair)") \
    X(precipitation_january, atmosphere, f32, grid, "mm/month", "January precipitation (legacy pair)") \
    X(precipitation_july, atmosphere, f32, grid, "mm/month", "July precipitation (legacy pair)") \
    X(mountain_precipitation_january, atmosphere, f32, grid, "mm/month", "January orographic precipitation") \
    X(mountain_precipitation_july, atmosphere, f32, grid, "mm/month", "July orographic precipitation") \
    X(mountain_precipitation_direction_january, atmosphere, i16, grid, "direction code", "January orographic rain direction") \
    X(mountain_precipitation_direction_july, atmosphere, i16, grid, "direction code", "July orographic rain direction") \
    X(prevailing_wind, atmosphere, i32, grid, "legacy code", "Prevailing wind direction") \
    X(climate_class, atmosphere, i16, grid, "code", "Climate classification") \
    X(subtropical_high_latitudes, atmosphere, i32, zonalprofile, "rows", "Horse-latitude rows per column") \
    X(temperature, atmosphere, f32, seasonal, "C", "Seasonal surface temperature") \
    X(precipitation, atmosphere, f32, seasonal, "mm/month", "Seasonal precipitation: the coupled solve, with the January and July adjustment passes of a procedural world") \
    X(pressure, atmosphere, f32, seasonal, "hPa", "Seasonal surface pressure anomaly") \
    X(wind_east, atmosphere, f32, seasonal, "m/s", "Seasonal surface eastward wind") \
    X(wind_south, atmosphere, f32, seasonal, "m/s", "Seasonal surface southward wind") \
    X(upper_height, atmosphere, f32, seasonal, "m", "Seasonal 500 hPa height anomaly") \
    X(upper_wind_east, atmosphere, f32, seasonal, "m/s", "Seasonal 500 hPa eastward wind") \
    X(upper_wind_south, atmosphere, f32, seasonal, "m/s", "Seasonal 500 hPa southward wind") \
    X(vertical_velocity, atmosphere, f32, seasonal, "hPa/day", "Seasonal 500 hPa ascent") \
    X(evaporation, atmosphere, f32, seasonal, "mm/day", "Seasonal ocean evaporation") \
    X(potential_evaporation, atmosphere, f32, seasonal, "mm/month", "Seasonal potential evapotranspiration of the rainfall solver") \
    X(maritime_influence, atmosphere, f32, seasonal, "0..1", "Seasonal maritime influence") \
    X(maritime_thermal_anomaly, atmosphere, f32, seasonal, "C", "Seasonal maritime thermal anomaly") \
    X(maritime_fetch, atmosphere, i16, seasonal, "cells", "Seasonal maritime fetch distance") \
    X(moisture, atmosphere, f32, seasonal, "mm", "Seasonal column water vapour") \
    X(convergence, atmosphere, f32, seasonal, "mm/month", "Seasonal moisture convergence") \
    X(orographic_uplift, atmosphere, f32, seasonal, "m", "Seasonal orographic parcel uplift") \
    X(orographic_subsidence, atmosphere, f32, seasonal, "m", "Seasonal orographic parcel subsidence") \
    X(ocean_current_east, hydrosphere, f32, seasonal, "m/s", "Seasonal eastward ocean current") \
    X(ocean_current_south, hydrosphere, f32, seasonal, "m/s", "Seasonal southward ocean current") \
    X(sea_surface_temperature, hydrosphere, f32, seasonal, "C", "Seasonal sea surface temperature; 0 on land") \
    X(solver_monthly_temperature, atmosphere, f32, monthly, "C", "Climate solver monthly temperature") \
    X(solver_monthly_precipitation, atmosphere, f32, monthly, "mm/month", "Climate solver monthly precipitation") \
    X(solver_surface_wind_east, atmosphere, f32, seasonal, "m/s", "Climate solver surface eastward wind") \
    X(solver_surface_wind_south, atmosphere, f32, seasonal, "m/s", "Climate solver surface southward wind") \
    X(solver_upper_wind_east, atmosphere, f32, seasonal, "m/s", "Climate solver 500 hPa eastward wind") \
    X(solver_upper_wind_south, atmosphere, f32, seasonal, "m/s", "Climate solver 500 hPa southward wind") \
    X(solver_low_wind_east, atmosphere, f32, seasonal, "m/s", "Climate solver 850 hPa eastward wind") \
    X(solver_low_wind_south, atmosphere, f32, seasonal, "m/s", "Climate solver 850 hPa southward wind") \
    X(solver_sea_surface_temperature, hydrosphere, f32, seasonal, "C", "Climate solver seasonal SST") \
    X(solver_sea_ice, cryosphere, f32, seasonal, "fraction", "Climate solver seasonal sea ice") \
    X(level_wind_east, atmosphere, f32, levelseasonal, "m/s", "Eastward wind per atmosphere level (table atmosphere.levels)") \
    X(level_wind_north, atmosphere, f32, levelseasonal, "m/s", "Northward wind per atmosphere level") \
    X(level_geopotential_height, atmosphere, f32, levelseasonal, "m", "Geopotential height per level; 500 hPa holds the model height anomaly") \
    X(level_vertical_velocity, atmosphere, f32, levelseasonal, "hPa/day", "Ascent per level, positive upwards") \
    X(level_available, atmosphere, u8, levelseasonal, "", "1 where the level is above ground and solved") \
    X(column_water, atmosphere, f32, seasonal, "mm", "Column water vapour on the moisture-process grid; diagnostic of the last solve") \
    X(moisture_transport_east, atmosphere, f32, seasonal, "kg/(m s)", "Eastward column moisture transport, process grid; diagnostic") \
    X(moisture_transport_north, atmosphere, f32, seasonal, "kg/(m s)", "Northward column moisture transport, process grid; diagnostic") \
    X(moisture_convergence, atmosphere, f32, seasonal, "mm/day", "Convergence of the column moisture transport, process grid; diagnostic") \
    X(precipitation_rate, atmosphere, f32, seasonal, "mm/day", "Surface precipitation, stratiform + orographic + convective, process grid; diagnostic") \
    X(ascent, atmosphere, f32, seasonal, "hPa/day", "Mean ascent, process grid; diagnostic") \
    X(current_east, hydrosphere, f32, seasonal, "m/s", "Eastward mixed-layer current on the ocean grid; diagnostic") \
    X(current_north, hydrosphere, f32, seasonal, "m/s", "Northward mixed-layer current on the ocean grid; diagnostic") \
    X(mixed_layer_temperature, hydrosphere, f32, seasonal, "C", "Season-mean mixed-layer temperature on the ocean grid; diagnostic") \
    X(storage_depth, hydrosphere, f32, seasonal, "m", "Season-mean thermal storage depth on the ocean grid; diagnostic") \
    X(biome, biosphere, i16, grid, "code", "Holdridge biome slot") \
    X(fishery_potential, biosphere, u8, grid, "0..100", "Fishery potential score") \
    X(settlement_suitability, anthroposphere, u8, grid, "0..100", "Settlement suitability") \
    X(infrastructure, anthroposphere, u8, grid, "0..100", "Infrastructure diffusion") \
    X(agricultural_capacity, anthroposphere, u8, grid, "0..100", "Agricultural capacity") \
    X(route_traffic, anthroposphere, u8, grid, "0..100", "Route traffic") \
    X(river_access, anthroposphere, u8, grid, "0..100", "River access") \
    X(harbor_score, anthroposphere, u8, grid, "0..100", "Harbour score") \
    X(owner_settlement, anthroposphere, i32, grid, "id", "Owning settlement") \
    X(owner_polity, anthroposphere, i32, grid, "id", "Owning polity") \
    X(no_shade, presentation, b8, grid, "", "Cells excluded from relief shading") \
    X(debug_test, presentation, i32, grid, "", "Debug scratch map")

enum class FieldId : std::uint16_t
{
#define UW_WORLD_FIELD_ENUM(id, component, element, layout, units, description) id,
    UW_WORLD_FIELD_LIST(UW_WORLD_FIELD_ENUM)
#undef UW_WORLD_FIELD_ENUM
};

inline constexpr std::size_t kFieldCount = 0
#define UW_WORLD_FIELD_COUNT(id, component, element, layout, units, description) + 1
    UW_WORLD_FIELD_LIST(UW_WORLD_FIELD_COUNT)
#undef UW_WORLD_FIELD_COUNT
    ;

// Entity tables are lists of records rather than grids.
enum class TableId : std::uint8_t
{
    tectonic_boundary_segments,
    tectonic_deforming_regions,
    mineral_deposits,
    mineral_budget,
    craters,
    settlements,
    polities,
    route_edges,
    history_events,
    weather_state,
    atmosphere_levels,
    tectonic_snapshot
};
inline constexpr std::size_t kTableCount = 12;

struct FieldSpec
{
    FieldId id;
    Component component;
    ElementType element;
    Layout layout;
    const char* name;        // e.g. "base_elevation"
    const char* units;
    const char* description;
};

struct TableSpec
{
    TableId id;
    Component component;
    const char* name;
    const char* description;
};

const FieldSpec& fieldspec(FieldId id);
const TableSpec& tablespec(TableId id);
const char* componentname(Component component);
const char* elementname(ElementType element);   // NumPy dtype string, e.g. "<i2"
std::size_t elementsize(ElementType element);
const char* layoutname(Layout layout);
std::string qualifiedname(FieldId id);          // "geosphere.base_elevation"
std::string qualifiedname(TableId id);
bool findfield(std::string_view qualified, FieldId& id);
bool findtable(std::string_view qualified, TableId& id);
bool findcomponent(std::string_view name, Component& component);
}
