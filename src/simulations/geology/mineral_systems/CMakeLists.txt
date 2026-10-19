# Mineral deposit history consumes engine-neutral tectonic frames; only the
# C API reader depends on the tectonics engine, so it is a library of its own.
add_library(uw_mineral_systems STATIC
    gold_formation.cpp
    mineral_deposit.cpp
    mineral_history.cpp
    placer_routing.cpp
    plate_markers.cpp
)
target_link_libraries(uw_mineral_systems PUBLIC uw_includes uw_core)
add_library(uw_mineral_tectonics STATIC tectonic_frame_reader.cpp)
target_link_libraries(uw_mineral_tectonics PUBLIC uw_mineral_systems uw_tectonics)
