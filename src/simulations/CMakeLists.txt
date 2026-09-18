# Numerical modules are shared by the desktop and headless test executables.
# Keep desktop headers and dependencies out of their public include paths.
add_library(uw_grid STATIC ../core/climate_grid.cpp)
find_package(Threads REQUIRED)
target_link_libraries(uw_grid PUBLIC Threads::Threads)
target_include_directories(uw_grid PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../core")

add_library(uw_climate STATIC
    climate/climate_atmosphere.cpp
    climate/climate_boundary_layer.cpp
    climate/climate_energy.cpp
    climate/climate_flow.cpp
    climate/climate_hydrology.cpp
    climate/climate_transport.cpp
    climate/climate_koppen.cpp
    climate/climate_ocean_dynamics.cpp
    climate/climate_physics.cpp
    climate/climate_weather.cpp
    climate/legacy/thermal_response.cpp
    climate/legacy/ocean_feedback.cpp
)
target_include_directories(uw_climate PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/climate")
# The energy-map adapter still reads planet and tuning; the other kernels use
# only their own value types and core grid/season definitions.
target_include_directories(uw_climate PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/../wip"
    "${CMAKE_CURRENT_SOURCE_DIR}/society"
)
target_link_libraries(uw_climate PUBLIC uw_grid)

add_library(uw_climate_io STATIC
    ../io/climate_reference.cpp
    ../io/climate_tiff.cpp
    ../io/climate_benchmark_outputs.cpp
)
target_include_directories(uw_climate_io PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../io")

add_library(uw_generation_stages STATIC ../pipeline/generation_stage_registry.cpp)
target_include_directories(uw_generation_stages PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../pipeline")

foreach(target uw_grid uw_climate uw_climate_io uw_generation_stages)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
endforeach()
