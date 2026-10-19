#pragma once

// Suite Landscape (docs/architecture/erosion-hydrology-plan.md, section 3, Tests):
// one fixture shared by landscape_tests.cpp (SphereMetric, FlowGraph and
// Depressions; step A1) and stream_power_tests.cpp (StreamPower; step A2). A test
// is TEST_F(Landscape, <Kernel><What>), so that
// ctest -R "^Landscape\.<Kernel>" selects one kernel's tests. The fixture holds
// nothing: each test builds its own synthetic grid of a few hundred cells.

#include <gtest/gtest.h>

#include <cstddef>

class Landscape : public ::testing::Test
{
};

// The lat-lon global column and row of a window's domain cell, from its level-0 ancestor
// (landscape::RegionalWindow::globalcell, landscape::LevelWindow::globalcell: J W + I), for the
// tests' raster checks.
template <typename Window>
int ancestorcolumn(const Window& window, int x)
{
    return static_cast<int>(window.globalcell(x, 0) % static_cast<std::size_t>(window.worldwidth));
}

template <typename Window>
int ancestorrow(const Window& window, int y)
{
    return static_cast<int>(window.globalcell(0, y) / static_cast<std::size_t>(window.worldwidth));
}
