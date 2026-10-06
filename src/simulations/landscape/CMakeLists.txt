# Erosion and drainage kernels on a cell-centred sphere grid, planet-free: no
# world class, no tuning header (the processes that call them pass the values).
# world_tests links this library (tests/CMakeLists.txt).
add_library(uw_landscape STATIC
    depressions.cpp
    flow_graph.cpp
    lake_balance.cpp
    sphere_metric.cpp
    stream_power.cpp
)
target_link_libraries(uw_landscape PUBLIC uw_includes PRIVATE uw_core)
