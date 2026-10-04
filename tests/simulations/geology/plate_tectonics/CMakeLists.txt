include(FetchContent)
FetchContent_Declare(googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
    GIT_SHALLOW TRUE
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(PlateTectonicsTests
    test_world_grid.cpp
    test_bounds.cpp
    test_boundary_motion.cpp
    test_heightmap.cpp
    test_mass.cpp
    test_movement.cpp
    test_phase3_contract.cpp
    test_phase4_contract.cpp
    test_phase5_contract.cpp
    test_phase7_tooling.cpp
    test_plate.cpp
    test_placer_transport.cpp
    test_portability.cpp
    test_randomness.cpp
    test_rectangle.cpp
    test_tectonic_dynamics.cpp
    test_tectonic_columns.cpp
    test_tectonic_transport.cpp
    test_tectonic_events.cpp
    test_tectonic_deformation.cpp
    test_tectonic_ocean_relief.cpp
    test_tectonic_noise_pass.cpp
    test_tectonic_margins.cpp
    test_sqrdmd.cpp
    test_tectonic_junctions.cpp
    test_topography_codec.cpp
    test_topography_io.cpp
)
target_link_libraries(PlateTectonicsTests PRIVATE
    plate_tectonics uw_tectonics_heightmap_io uw_tectonics_snapshots
    GTest::gtest GTest::gtest_main)
target_compile_definitions(PlateTectonicsTests PRIVATE
    PLATE_TECTONICS_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
uw_tectonics_target_defaults(PlateTectonicsTests)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts")
add_test(NAME PlateTectonicsTests COMMAND PlateTectonicsTests)
set_tests_properties(PlateTectonicsTests PROPERTIES LABELS "tectonics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts")
