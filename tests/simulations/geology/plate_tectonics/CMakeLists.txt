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
    test_acceptance.cpp
    test_bounds.cpp
    test_boundary_motion.cpp
    test_heightmap.cpp
    test_mass.cpp
    test_movement.cpp
    test_phase0_contract.cpp
    test_phase2_contract.cpp
    test_phase3_contract.cpp
    test_phase4_contract.cpp
    test_phase5_contract.cpp
    test_phase7_tooling.cpp
    test_plate.cpp
    test_placer_transport.cpp
    test_portability.cpp
    test_randomness.cpp
    test_rectangle.cpp
    test_regression.cpp
    test_surface_processes.cpp
    test_sqrdmd.cpp
    test_tectonic_junctions.cpp
    test_topography_codec.cpp
    test_topography_io.cpp
)
target_link_libraries(PlateTectonicsTests PRIVATE
    plate_tectonics PlateTectonicsHeightmapIO PlateTectonicsPipelineSupport
    GTest::gtest GTest::gtest_main)
target_compile_definitions(PlateTectonicsTests PRIVATE
    PLATE_TECTONICS_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/fixtures")
uw_tectonics_target_defaults(PlateTectonicsTests)

if(MSVC)
    # These historical tests retain their old baselines after GTEST_SKIP().
    # Debug MSVC reports the intentionally unreachable bodies as C4702.
    set_source_files_properties(test_phase0_contract.cpp test_regression.cpp
        PROPERTIES COMPILE_OPTIONS /wd4702)
endif()

set(tectonics_slow_filter "Regression.SimulationSeed12345_IsDeterministicAndWithinRange")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts")
add_test(NAME PlateTectonicsTests COMMAND PlateTectonicsTests
    "--gtest_filter=-${tectonics_slow_filter}")
set_tests_properties(PlateTectonicsTests PROPERTIES
    LABELS "tectonics"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts")

# Full-size deterministic replay and PNG artifacts are an explicit slow tier.
add_custom_target(check-tectonics-slow
    COMMAND $<TARGET_FILE:PlateTectonicsTests> "--gtest_filter=${tectonics_slow_filter}"
    DEPENDS PlateTectonicsTests
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/artifacts"
    USES_TERMINAL)
