# The world class (planet components and facade), regions and shared utilities.
add_library(uw_world STATIC
    planet.cpp
    planet_bundle.cpp
    planet_fields.cpp
    region.cpp
    region_queries.cpp
    regionalmap.cpp
    world_setup.cpp
)
target_link_libraries(uw_world PUBLIC uw_includes)
