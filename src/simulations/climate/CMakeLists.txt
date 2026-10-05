# Climate solvers on their own value types, and the climate processes on a
# climateworld, one source per process id (climate_pressure.cpp is climate.pressure). The energy-map adapter still
# reads planet and tuning.
add_library(uw_climate STATIC
    climate_atmosphere.cpp
    climate_boundary_layer.cpp
    climate_classification.cpp
    climate_diagnostics.cpp
    climate_energy.cpp
    climate_fields.cpp
    climate_flow.cpp
    climate_hydrology.cpp
    climate_koppen.cpp
    climate_ocean_currents.cpp
    climate_ocean_dynamics.cpp
    climate_physics.cpp
    climate_pressure.cpp
    climate_rainfall_and_coupling.cpp
    climate_sea_surface_temperature.cpp
    climate_transport.cpp
    climate_weather.cpp
    climate_winds.cpp
    climate_world.cpp
    climate_world_fields.cpp
    legacy/ocean_feedback.cpp
    legacy/thermal_response.cpp
)
target_link_libraries(uw_climate PUBLIC uw_includes)
