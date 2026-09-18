# Planet-level stages still share the legacy world/utilities boundary. Object
# libraries keep each domain's sources explicit while those callbacks remain
# supplied by the desktop. Pure numerical solvers live in uw_climate instead.
function(uw_world_module name)
    add_library(${name} OBJECT ${ARGN})
    target_include_directories(${name} PRIVATE ${UW_INCLUDE_DIRS})
    target_link_libraries(${name} PRIVATE sfml-graphics)
    set_target_properties(${name} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
endfunction()

uw_world_module(uw_geology src/simulations/geology/plate_tectonics_adapter.cpp)
target_link_libraries(uw_geology PRIVATE plate_tectonics)

uw_world_module(uw_terrain
    src/simulations/terrain/assetdata.cpp
    src/simulations/terrain/fastlem_mountains.cpp
    src/simulations/terrain/terrain_fractal.cpp
    src/simulations/terrain/terrain_shapes.cpp
    src/simulations/terrain/terrain_coasts.cpp
    src/simulations/terrain/terrain_ridges.cpp
    src/simulations/terrain/terrain_elevation.cpp
    src/simulations/terrain/terrain_shelves.cpp
    src/simulations/terrain/terrain_ocean_features.cpp
    src/simulations/terrain/terrain_mountain_import.cpp
    src/simulations/terrain/terrain_volcanoes.cpp
    src/simulations/terrain/terrain_craters.cpp
    src/simulations/terrain/climate_landforms.cpp
    src/simulations/terrain/land_distance.cpp
    src/simulations/terrain/surface_potentials.cpp
)

uw_world_module(uw_world_climate
    src/simulations/climate/climate_fields.cpp
    src/simulations/climate/climate_classification.cpp
    src/simulations/climate/climate_coupling.cpp
)
target_link_libraries(uw_world_climate PRIVATE uw_climate)

uw_world_module(uw_hydrology
    src/simulations/hydrology/drainage.cpp
    src/simulations/hydrology/drainage_basins.cpp
    src/simulations/hydrology/lakes.cpp
    src/simulations/hydrology/lake_effects.cpp
    src/simulations/hydrology/deltas.cpp
    src/simulations/hydrology/wetlands.cpp
    src/simulations/hydrology/tides.cpp
)

uw_world_module(uw_resources
    src/simulations/resources/mineral_resources.cpp
    src/simulations/resources/marine_resources.cpp
)

uw_world_module(uw_society
    src/simulations/society/social_sites.cpp
    src/simulations/society/social_routes.cpp
    src/simulations/society/social_polities.cpp
    src/simulations/society/social_trade.cpp
)
