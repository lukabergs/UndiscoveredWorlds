# The world class (planet components and facade), regions and shared utilities.
add_library(uw_world STATIC
    classes.cpp
    misc.cpp
    planet.cpp
    planet_bundle.cpp
    planet_fields.cpp
    region.cpp
    regionalmap.cpp
)
target_link_libraries(uw_world PUBLIC uw_includes)
