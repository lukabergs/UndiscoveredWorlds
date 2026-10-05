add_library(uw_terrain STATIC
    climate_erosion.cpp
    climate_landforms.cpp
    coast_queries.cpp
    fastlem_mountains.cpp
    shape_templates.cpp
    surface_potentials.cpp
    terrain_coasts.cpp
    terrain_craters.cpp
    terrain_elevation.cpp
    terrain_mountain_import.cpp
    terrain_ocean_features.cpp
    terrain_ridges.cpp
    terrain_shapes.cpp
    terrain_shelves.cpp
    terrain_volcanoes.cpp
)
target_link_libraries(uw_terrain PUBLIC uw_includes)
