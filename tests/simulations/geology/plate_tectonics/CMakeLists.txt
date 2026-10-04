add_executable(tectonics_tests
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
target_link_libraries(tectonics_tests PRIVATE
    plate_tectonics uw_tectonics_heightmap_io uw_tectonics_snapshots
    GTest::gtest GTest::gtest_main)
uw_tectonics_target_defaults(tectonics_tests)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts")
gtest_discover_tests(tectonics_tests
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts"
    PROPERTIES LABELS "tectonics")
