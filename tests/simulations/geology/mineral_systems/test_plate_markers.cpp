#include "gtest/gtest.h"

#include "frame_fixture.hpp"
#include "plate_markers.hpp"

using namespace mineralsystems;

namespace
{
Marker markerat(double x, double y)
{
    Marker marker;
    marker.x = x;
    marker.y = y;
    return marker;
}

// Plate 0 owns x < split, plate 1 the rest.
void splitplates(FrameFixture& fixture, std::uint32_t split, std::uint64_t leftuid, std::uint64_t rightuid)
{
    fixture.frame.plates = { { leftuid }, { rightuid } };
    for (std::uint32_t y = 0; y < fixture.frame.height; ++y)
        for (std::uint32_t x = 0; x < fixture.frame.width; ++x)
            fixture.plate[fixture.cell(x, y)] = x < split ? 0u : 1u;
}
}

TEST(PlateMarkers, TranslationFollowsTheRecordedShiftAndWraps)
{
    FrameFixture fixture(16, 8);
    fixture.frame.plates = { { 11 } };
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(3.5, 4.5), markerat(15.5, 0.1) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    ASSERT_TRUE(tracker.attach(markers[1], fixture.view()));

    fixture.motion(1.5f, -0.25f);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.moved, 2u);
    EXPECT_NEAR(markers[0].x, 5.0, 1e-12);
    EXPECT_NEAR(markers[0].y, 4.25, 1e-12);
    EXPECT_NEAR(markers[1].x, 1.0, 1e-12);
    EXPECT_NEAR(markers[1].y, 7.85, 1e-12);
}

TEST(PlateMarkers, ContinentalMarkerTransfersWithCollidingCrust)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(3.5, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    EXPECT_TRUE(markers[0].continentalHost);

    fixture.motion(1.0f);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.rehomed, 1u);
    EXPECT_EQ(markers[0].hostUid, 2u);
    EXPECT_EQ(markers[0].state, MarkerState::attached);
}

TEST(PlateMarkers, ContinentalMarkerReturnsFromAnOceanicRasterSeam)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    for (std::uint32_t y = 0; y < 4; ++y)
        for (std::uint32_t x = 4; x < 8; ++x)
            fixture.crust[fixture.cell(x, y)] = crustcode::oceanic;
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(3.9, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));

    fixture.motion(0.3f);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.snapped, 1u);
    EXPECT_EQ(markers[0].hostUid, 1u);
    EXPECT_DOUBLE_EQ(markers[0].x, 3.75);
    EXPECT_DOUBLE_EQ(markers[0].y, 1.5);
}

TEST(PlateMarkers, OceanicMarkerIsSubductedUnderOlderCrustOfAnExistingPlate)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    fixture.frame.timeMyr = 100.0;
    for (std::size_t cell = 0; cell < fixture.plate.size(); ++cell)
    {
        fixture.crust[cell] = crustcode::oceanic;
        fixture.age[cell] = fixture.plate[cell] == 0 ? 10.0f : 50.0f;
    }
    MarkerTracker tracker(3.0);
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(3.5, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    EXPECT_FALSE(markers[0].continentalHost);
    EXPECT_DOUBLE_EQ(markers[0].hostCrustBirthMyr, 90.0);

    fixture.frame.timeMyr = 101.0;
    fixture.motion(1.0f);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.subducted, 1u);
    EXPECT_EQ(markers[0].state, MarkerState::subducted);
    EXPECT_DOUBLE_EQ(markers[0].stateTimeMyr, 101.0);
}

TEST(PlateMarkers, RelabelledOceanicCrustKeepsItsMarker)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    fixture.frame.timeMyr = 100.0;
    for (std::size_t cell = 0; cell < fixture.plate.size(); ++cell)
    {
        fixture.crust[cell] = crustcode::oceanic;
        fixture.age[cell] = 10.0f;
    }
    MarkerTracker tracker(3.0);
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(3.5, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));

    // The region changes owner with its crust; ages advance with time.
    fixture.frame.timeMyr = 101.0;
    for (float& age : fixture.age)
        age = 11.0f;
    fixture.plate[fixture.cell(3, 1)] = 1;
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.subducted, 0u);
    EXPECT_EQ(markers[0].hostUid, 2u);
}

TEST(PlateMarkers, NewbornPlateInheritsMarkersWithoutSubduction)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    for (std::size_t cell = 0; cell < fixture.plate.size(); ++cell)
    {
        fixture.crust[cell] = crustcode::oceanic;
        fixture.age[cell] = fixture.plate[cell] == 0 ? 10.0f : 50.0f;
    }
    MarkerTracker tracker(3.0);
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(1.5, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));

    fixture.frame.plates.push_back({ 3 });
    fixture.plate[fixture.cell(1, 1)] = 2;
    fixture.age[fixture.cell(1, 1)] = 70.0f;
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.subducted, 0u);
    EXPECT_EQ(markers[0].hostUid, 3u);
}

TEST(PlateMarkers, UncoveredCellsCannotBeAttached)
{
    FrameFixture fixture(4, 2);
    fixture.frame.plates = { { 1 } };
    fixture.plate[fixture.cell(1, 1)] = 0xFFFFFFFFu;
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    Marker marker = markerat(1.5, 1.5);
    EXPECT_FALSE(tracker.attach(marker, fixture.view()));
}
