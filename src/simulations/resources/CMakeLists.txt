add_library(uw_resources STATIC
    gold_deposits.cpp
    marine_resources.cpp
    mineral_resources.cpp
)
target_link_libraries(uw_resources PUBLIC uw_includes uw_mineral_systems)
