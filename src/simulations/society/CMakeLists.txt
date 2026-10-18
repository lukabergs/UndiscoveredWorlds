add_library(uw_society STATIC
    social_geometry.cpp
    social_polities.cpp
    social_routes.cpp
    social_sites.cpp
    social_trade.cpp
)
target_link_libraries(uw_society PUBLIC uw_includes)
