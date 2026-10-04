add_executable(minerals_tests
    test_plate_markers.cpp
    test_gold_formation.cpp
    test_mineral_history.cpp
    test_placer_routing.cpp
    test_lithosphere_markers.cpp
)
target_link_libraries(minerals_tests PRIVATE uw_mineral_tectonics GTest::gtest GTest::gtest_main)
gtest_discover_tests(minerals_tests PROPERTIES LABELS "geology")
