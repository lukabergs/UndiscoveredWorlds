# The legacy regional map generator, moved from src/app in erosion plan E1: integer and
# factor 16 (consolidation plan section 11, "Later"). The seamless regional pass of
# Phase E is in uw_landscape.
add_library(uw_regional STATIC
    region.cpp
    region_queries.cpp
    regional_map.cpp
)
target_link_libraries(uw_regional PUBLIC uw_includes)
