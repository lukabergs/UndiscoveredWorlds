# Planet-level stages still share the world class and utilities in src/world.
# Object libraries keep each domain's sources explicit; none of them depends on
# SFML or ImGui. Pure numerical solvers live in uw_climate instead.
function(uw_world_module name)
    add_library(${name} OBJECT ${ARGN})
    target_include_directories(${name} PRIVATE ${UW_INCLUDE_DIRS})
endfunction()

uw_world_module(uw_geology src/simulations/geology/world_tectonics.cpp)
target_link_libraries(uw_geology PRIVATE plate_tectonics uw_mineral_tectonics uw_tectonics_snapshots)

uw_world_module(uw_terrain
    src/simulations/terrain/assetdata.cpp
    src/simulations/terrain/fastlem_mountains.cpp
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
    src/simulations/terrain/surface_potentials.cpp
)

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
    src/simulations/resources/gold_deposits.cpp
    src/simulations/resources/marine_resources.cpp
)
target_link_libraries(uw_resources PRIVATE uw_mineral_systems)

uw_world_module(uw_society
    src/simulations/society/social_sites.cpp
    src/simulations/society/social_routes.cpp
    src/simulations/society/social_polities.cpp
    src/simulations/society/social_trade.cpp
)
