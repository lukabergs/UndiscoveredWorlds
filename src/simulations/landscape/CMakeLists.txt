# Erosion and drainage kernels on a cell-centred sphere grid, planet-free: no
# world class, no tuning header (the processes that call them pass the values;
# one exception: spherical_noise.cpp reads landscape.noise.* for a FractalNoise left
# at "configured", multi-scale plan S3.16). world_tests links this library
# (tests/CMakeLists.txt).
add_library(uw_landscape STATIC
    analytical_relief.cpp
    channel_network.cpp
    depressions.cpp
    distributaries.cpp
    flow_graph.cpp
    glacial.cpp
    isostasy.cpp
    lake_balance.cpp
    level_blocks.cpp
    ocean_mask.cpp
    opensimplex2s.cpp
    refinement_grid.cpp
    refinement_levels.cpp
    refinement_projection.cpp
    regional_cirques.cpp
    regional_coast.cpp
    regional_edges.cpp
    regional_folds.cpp
    regional_glacial.cpp
    regional_grid.cpp
    regional_ice.cpp
    regional_pass.cpp
    regional_solve.cpp
    regional_surface.cpp
    regional_tributaries.cpp
    regional_troughs.cpp
    regional_valleys.cpp
    regional_water.cpp
    relief_roughness.cpp
    runoff.cpp
    sphere_metric.cpp
    spherical_noise.cpp
    stream_power.cpp
    subgrid_statistics.cpp
    upsampling.cpp
)
target_link_libraries(uw_landscape PUBLIC uw_includes PRIVATE uw_core)
