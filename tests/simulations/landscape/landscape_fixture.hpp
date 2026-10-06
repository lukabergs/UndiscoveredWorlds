#pragma once

// Suite Landscape (docs/architecture/erosion-hydrology-plan.md, section 3, Tests):
// one fixture shared by landscape_tests.cpp (SphereMetric, FlowGraph and
// Depressions; step A1) and stream_power_tests.cpp (StreamPower; step A2). A test
// is TEST_F(Landscape, <Kernel><What>), so that
// ctest -R "^Landscape\.<Kernel>" selects one kernel's tests. The fixture holds
// nothing: each test builds its own synthetic grid of a few hundred cells.

#include <gtest/gtest.h>

class Landscape : public ::testing::Test
{
};
