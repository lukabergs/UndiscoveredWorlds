# The world class (planet components and facade) and shared utilities.
add_library(uw_world STATIC
    planet.cpp
    planet_bundle.cpp
    planet_export.cpp
    planet_fields.cpp
    world_setup.cpp
)
target_link_libraries(uw_world PUBLIC uw_includes)
