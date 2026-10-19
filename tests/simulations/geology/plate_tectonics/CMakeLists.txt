add_executable(tectonics_tests
    test_tectonic_grid.cpp
    test_boundary_motion.cpp
    test_phase3_contract.cpp
    test_phase4_contract.cpp
    test_phase5_contract.cpp
    test_phase7_tooling.cpp
    test_placer_transport.cpp
    test_portability.cpp
    test_randomness.cpp
    test_tectonic_dynamics.cpp
    test_tectonic_columns.cpp
    test_tectonic_sphere_transport.cpp
    test_tectonic_events.cpp
    test_tectonic_deformation.cpp
    test_tectonic_ocean_relief.cpp
    test_tectonic_noise_pass.cpp
    test_tectonic_margins.cpp
    test_tectonic_options.cpp
    test_tectonic_junctions.cpp
    test_topography_codec.cpp
    test_topography_io.cpp
    # The scenario options of world_pipeline need only the engine.
    ${PROJECT_SOURCE_DIR}/src/cli/tectonic_options.cpp
)
target_link_libraries(tectonics_tests PRIVATE
    uw_tectonics uw_io uw_validation
    GTest::gtest GTest::gtest_main)
uw_tectonics_target_defaults(tectonics_tests)

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts")
gtest_discover_tests(tectonics_tests
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts"
    PROPERTIES LABELS "tectonics")
