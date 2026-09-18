#include "gtest/gtest.h"
#include "../../../../src/simulations/resources/placer_transport.hpp"

using physical_layers::detail::route_ore_potential;

TEST(PlacerTransport, DistantConnectedSourceReachesOutletButNotAdjacentCatchment) {
    const std::vector<int32_t> receiver{1, 2, 3, 4, 5, -1, 7, -1};
    const auto result = route_ore_potential(receiver, {80, 0, 0, 0, 0, 0, 20, 0});
    EXPECT_EQ(result, (std::vector<uint8_t>{80, 80, 80, 80, 80, 80, 20, 20}));
}

TEST(PlacerTransport, TributariesMergeWithoutSendingOreUpstream) {
    EXPECT_EQ(route_ore_potential({2, 2, 3, -1}, {20, 90, 0, 0}),
              (std::vector<uint8_t>{20, 90, 90, 90}));
}

TEST(PlacerTransport, CyclesAndInvalidOutletsTerminateDeterministically) {
    EXPECT_EQ(route_ore_potential({1, 2, 1, 99, 4}, {70, 0, 20, 5, 30}),
              (std::vector<uint8_t>{70, 70, 70, 5, 30}));
    EXPECT_TRUE(route_ore_potential({}, {}).empty());
    EXPECT_THROW(route_ore_potential({-1}, {}), std::invalid_argument);
}

TEST(PlacerTransport, CatchmentWithoutOreHasZeroPotential) {
    EXPECT_EQ(route_ore_potential({1, 2, -1}, {0, 0, 0}),
              (std::vector<uint8_t>{0, 0, 0}));
}
