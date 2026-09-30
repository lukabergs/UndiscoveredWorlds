add_executable(MineralSystemsTests
    test_plate_markers.cpp
    test_gold_formation.cpp
    test_mineral_history.cpp
    test_placer_routing.cpp
    test_lithosphere_markers.cpp
)
target_link_libraries(MineralSystemsTests PRIVATE uw_mineral_tectonics GTest::gtest GTest::gtest_main)
set_target_properties(MineralSystemsTests PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)
add_test(NAME MineralSystemsTests COMMAND MineralSystemsTests)
set_tests_properties(MineralSystemsTests PROPERTIES LABELS "geology")
