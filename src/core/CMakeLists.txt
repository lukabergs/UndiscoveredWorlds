# Grids, noise, random numbers, the field catalogue, environment paths, profiler.
add_library(uw_core STATIC
    app_environment.cpp
    classes.cpp
    climate_grid.cpp
    grid_coordinates.cpp
    grid_shapes.cpp
    parameters.cpp
    procedural_noise.cpp
    profiler.cpp
    text_values.cpp
    world/world_fields.cpp
    world/world_grid_check.cpp
)
target_link_libraries(uw_core PUBLIC uw_includes Threads::Threads)
