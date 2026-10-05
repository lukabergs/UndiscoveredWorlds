# Climate solvers on their own value types, and the climate processes on a
# climateworld (coupling, fields, classification). The energy-map adapter still
# reads planet and tuning.
add_library(uw_climate STATIC
    climate_atmosphere.cpp
    climate_boundary_layer.cpp
    climate_classification.cpp
    climate_coupling.cpp
    climate_energy.cpp
    climate_fields.cpp
    climate_flow.cpp
    climate_hydrology.cpp
    climate_koppen.cpp
    climate_ocean_dynamics.cpp
    climate_physics.cpp
    climate_transport.cpp
    climate_weather.cpp
    climate_world.cpp
    climate_world_fields.cpp
    legacy/ocean_feedback.cpp
    legacy/thermal_response.cpp
)
target_link_libraries(uw_climate PUBLIC uw_includes)
