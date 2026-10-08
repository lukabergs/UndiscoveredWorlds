# Erosion and drainage kernels on a cell-centred sphere grid, planet-free: no
# world class, no tuning header (the processes that call them pass the values).
# world_tests links this library (tests/CMakeLists.txt).
add_library(uw_landscape STATIC
    analytical_relief.cpp
    channel_network.cpp
    depressions.cpp
    distributaries.cpp
    flow_graph.cpp
    glacial.cpp
    isostasy.cpp
    lake_balance.cpp
    ocean_mask.cpp
    regional_grid.cpp
    regional_pass.cpp
    regional_surface.cpp
    regional_water.cpp
    relief_roughness.cpp
    runoff.cpp
    sphere_metric.cpp
    spherical_noise.cpp
    stream_power.cpp
    upsampling.cpp
)
target_link_libraries(uw_landscape PUBLIC uw_includes PRIVATE uw_core)
