add_library(uw_hydrology STATIC
    deltas.cpp
    drainage.cpp
    drainage_basins.cpp
    lake_effects.cpp
    lakes.cpp
    tides.cpp
    wetlands.cpp
)
target_link_libraries(uw_hydrology PUBLIC uw_includes)
