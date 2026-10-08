include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
include(GoogleTest)

add_subdirectory(simulations/geology/plate_tectonics)
add_subdirectory(simulations/geology/mineral_systems)

if(TARGET uw_climate_world)
    add_executable(world_tests
        core/climate_grid_tests.cpp
        core/parallel_rows_tests.cpp
        io/climate_reference_tests.cpp
        io/climate_tiff_tests.cpp
        pipeline/climate_pipeline_tests.cpp
        pipeline/climate_world_tests.cpp
        simulations/climate/climate_atmosphere_tests.cpp
        simulations/climate/climate_boundary_layer_tests.cpp
        simulations/climate/climate_context_tests.cpp
        simulations/climate/climate_energy_tests.cpp
        simulations/climate/climate_flow_tests.cpp
        simulations/climate/climate_hydrology_tests.cpp
        simulations/climate/climate_koppen_tests.cpp
        simulations/climate/climate_ocean_dynamics_tests.cpp
        simulations/climate/climate_ocean_heat_tests.cpp
        simulations/climate/climate_physics_tests.cpp
        simulations/climate/climate_transport_tests.cpp
        simulations/climate/climate_weather_tests.cpp
        simulations/landscape/analytical_relief_tests.cpp
        simulations/landscape/channel_network_tests.cpp
        simulations/landscape/distributaries_tests.cpp
        simulations/landscape/glacial_tests.cpp
        simulations/landscape/isostasy_tests.cpp
        simulations/landscape/lake_balance_tests.cpp
        simulations/landscape/ocean_mask_tests.cpp
        simulations/landscape/regional_grid_tests.cpp
        simulations/landscape/regional_pass_tests.cpp
        simulations/landscape/regional_surface_tests.cpp
        simulations/landscape/regional_water_tests.cpp
        simulations/landscape/relief_roughness_tests.cpp
        simulations/landscape/routing_options_tests.cpp
        simulations/landscape/runoff_tests.cpp
        simulations/landscape/spherical_noise_tests.cpp
        simulations/landscape/landscape_tests.cpp
        simulations/landscape/stream_power_tests.cpp
        simulations/landscape/torus_metric_tests.cpp
        simulations/landscape/upsampling_tests.cpp
    )
    # Planet-free libraries only on its link line: climate, its stage runner, its
    # validation and the landscape kernels must not need the world class.
    target_link_libraries(world_tests PRIVATE uw_climate_world uw_landscape GTest::gtest GTest::gtest_main)
    target_include_directories(world_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    # A deadlock of the row executor must fail quickly.
    gtest_discover_tests(world_tests TEST_FILTER "ParallelRows.*" PROPERTIES TIMEOUT 15)
    gtest_discover_tests(world_tests TEST_FILTER "-ParallelRows.*")
else()
    message(STATUS "world_tests is not built: it needs UW_BUILD_WORLD or UW_BUILD_APP")
endif()
