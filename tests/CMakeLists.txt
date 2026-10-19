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
        core/mesh_tests.cpp
        core/mesh_window_tests.cpp
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
        simulations/landscape/cube_kernels_tests.cpp
        simulations/landscape/distributaries_tests.cpp
        simulations/landscape/glacial_tests.cpp
        simulations/landscape/isostasy_tests.cpp
        simulations/landscape/lake_balance_tests.cpp
        simulations/landscape/ocean_mask_tests.cpp
        simulations/landscape/refinement_grid_tests.cpp
        simulations/landscape/refinement_levels_tests.cpp
        simulations/landscape/refinement_projection_tests.cpp
        simulations/landscape/regional_cirques_tests.cpp
        simulations/landscape/regional_coast_tests.cpp
        simulations/landscape/regional_edges_tests.cpp
        simulations/landscape/regional_folds_tests.cpp
        simulations/landscape/regional_glacial_tests.cpp
        simulations/landscape/regional_grid_tests.cpp
        simulations/landscape/regional_ice_tests.cpp
        simulations/landscape/regional_pass_tests.cpp
        simulations/landscape/regional_solve_tests.cpp
        simulations/landscape/regional_surface_tests.cpp
        simulations/landscape/regional_tributaries_tests.cpp
        simulations/landscape/regional_troughs_tests.cpp
        simulations/landscape/regional_valleys_tests.cpp
        simulations/landscape/regional_water_tests.cpp
        simulations/landscape/relief_roughness_tests.cpp
        simulations/landscape/routing_options_tests.cpp
        simulations/landscape/runoff_tests.cpp
        simulations/landscape/spherical_noise_tests.cpp
        simulations/landscape/landscape_tests.cpp
        simulations/landscape/stream_power_tests.cpp
        simulations/landscape/subgrid_statistics_tests.cpp
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

if(TARGET world_pipeline)
    # The world class on its cell grid (docs/architecture/sphere-mesh.md section 4): needs planet and the pipeline.
    add_executable(planet_grid_tests
        world/planet_grid_tests.cpp
    )
    target_link_libraries(planet_grid_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(planet_grid_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(planet_grid_tests)

    # The coast, shelf, ocean-feature, arid-landform and import terrain functions on a cube
    # planet (sphere-mesh.md section 12, T8 W2-T2).
    add_executable(terrain_cube_coast_tests
        simulations/terrain/terrain_cube_coast_tests.cpp
    )
    target_link_libraries(terrain_cube_coast_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(terrain_cube_coast_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(terrain_cube_coast_tests)

    # Bundle v2 and the lat-lon export of a world (sphere-mesh.md sections 6 and 7).
    add_executable(bundle_v2_tests
        io/bundle_v2_tests.cpp
    )
    target_link_libraries(bundle_v2_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(bundle_v2_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(bundle_v2_tests)

    # The lat-lon export of a cube world for its consumers (sphere-mesh.md sections 6, 7 and 9; W2-IO).
    add_executable(export_tests
        io/export_tests.cpp
    )
    target_link_libraries(export_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(export_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(export_tests)
    # The terrain processes on a cube world (sphere-mesh.md sections 3, 11 and 12; T8 W2-T1).
    add_executable(terrain_cube_tests
        simulations/terrain/terrain_cube_tests.cpp
    )
    target_link_libraries(terrain_cube_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(terrain_cube_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(terrain_cube_tests)

    # Society on a cube planet (sphere-mesh.md section 12, W2-S).
    add_executable(society_cube_tests
        simulations/society/society_cube_tests.cpp
    )
    target_link_libraries(society_cube_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(society_cube_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(society_cube_tests)

    # The resource layers on a cube planet (sphere-mesh.md section 12, W2-R).
    add_executable(resources_cube_tests
        simulations/resources/resources_cube_tests.cpp
    )
    target_link_libraries(resources_cube_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(resources_cube_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(resources_cube_tests)
    # Hydrology on a cube planet (sphere-mesh.md section 12, W2-H).
    add_executable(hydrology_cube_tests
        simulations/hydrology/hydrology_cube_tests.cpp
    )
    target_link_libraries(hydrology_cube_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(hydrology_cube_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(hydrology_cube_tests)

    # The climate adapter on a cube world (sphere-mesh.md section 10).
    add_executable(climate_adapter_tests
        pipeline/climate_adapter_tests.cpp
    )
    target_link_libraries(climate_adapter_tests PRIVATE uw_headless GTest::gtest GTest::gtest_main)
    target_include_directories(climate_adapter_tests PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    gtest_discover_tests(climate_adapter_tests)
endif()
