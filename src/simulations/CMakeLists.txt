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
add_library(uw_climate_erosion STATIC terrain/climate_erosion.cpp)
target_include_directories(uw_climate_erosion PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/terrain")
target_include_directories(uw_generation_stages PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../pipeline")

# Mineral deposit history consumes engine-neutral tectonic frames; only the
# C API reader depends on the tectonics engine.
add_library(uw_mineral_systems STATIC
    geology/mineral_systems/gold_formation.cpp
    geology/mineral_systems/mineral_deposit.cpp
    geology/mineral_systems/mineral_history.cpp
    geology/mineral_systems/placer_routing.cpp
    geology/mineral_systems/plate_markers.cpp
)
target_include_directories(uw_mineral_systems PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/geology/mineral_systems")
target_link_libraries(uw_mineral_systems PRIVATE uw_grid)
add_library(uw_mineral_tectonics STATIC geology/mineral_systems/tectonic_frame_reader.cpp)
target_link_libraries(uw_mineral_tectonics PUBLIC uw_mineral_systems plate_tectonics)

foreach(target uw_grid uw_climate uw_climate_io uw_generation_stages uw_climate_erosion uw_mineral_systems uw_mineral_tectonics)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
endforeach()
