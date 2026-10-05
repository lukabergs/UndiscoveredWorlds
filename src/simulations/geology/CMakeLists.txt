# The tectonics step of the world: runs the engine and writes its fields and deposits.
add_library(uw_geology STATIC world_tectonics.cpp)
target_link_libraries(uw_geology PUBLIC uw_includes uw_tectonics uw_mineral_tectonics)
