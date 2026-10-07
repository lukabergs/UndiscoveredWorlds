add_library(uw_hydrology STATIC
    drainage_basins.cpp
    flow_queries.cpp
    flow_routing.cpp
    lake_levels.cpp
    lakes.cpp
    rivers.cpp
    tides.cpp
    wetlands.cpp
)
target_link_libraries(uw_hydrology PUBLIC uw_includes)
