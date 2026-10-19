#include "gtest/gtest.h"

#include "frame_fixture.hpp"
#include "simulations/geology/mineral_systems/plate_markers.hpp"

#include <cmath>

using namespace mineralsystems;

namespace
{
Marker markerat(const FrameFixture& fixture, double x, double y)
{
    Marker marker;
    marker.position = fixture.at(x, y);
    return marker;
}

// Plate 0 owns x < split, plate 1 the rest.
void splitplates(FrameFixture& fixture, std::uint32_t split, std::uint64_t leftuid, std::uint64_t rightuid)
{
    fixture.frame.plates = { { leftuid }, { rightuid } };
    for (std::uint32_t y = 0; y < fixture.height(); ++y)
        for (std::uint32_t x = 0; x < fixture.width(); ++x)
            fixture.plate[fixture.cell(x, y)] = x < split ? 0u : 1u;
}
}

TEST(PlateMarkers, TranslationFollowsTheRecordedShiftAndCrossesTheDateLine)
{
    FrameFixture fixture(16, 8);
    fixture.frame.plates = { { 11 } };
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    // Near the equator, so (east, north) cells move a marker by about as many columns and rows.
    std::vector<Marker> markers{ markerat(fixture, 3.5, 4.25), markerat(fixture, 15.5, 4.25) };
    const mesh::Vec3 start = markers[0].position;
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    ASSERT_TRUE(tracker.attach(markers[1], fixture.view()));

    fixture.motion(1.5, 0.25);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.moved, 2u);
    // The marker travels the shift's length along the sphere (to the projection of each
    // cell's tangent onto the marker's), and stays a unit vector.
    const double travelled = std::acos(std::min(1.0, mesh::dot(start, markers[0].position))) * FrameFixture::kRadiusKm;
    EXPECT_NEAR(travelled / (std::hypot(1.5, 0.25) * fixture.cellkm()), 1.0, 0.02);
    EXPECT_NEAR(mesh::length(markers[0].position), 1.0, 1e-12);
    const auto a = fixture.coordinates(markers[0].position);
    const auto b = fixture.coordinates(markers[1].position);
    EXPECT_NEAR(a[0], 5.0, 0.1);
    EXPECT_NEAR(a[1], 4.0, 0.1);
    EXPECT_NEAR(b[0], 1.0, 0.1);  // across the date line
    EXPECT_NEAR(b[1], 4.0, 0.1);
}

TEST(PlateMarkers, RigidRotationKeepsHighLatitudeMarkersOnTheirSmallCircle)
{
    // A plate turning about the pole by a quarter column per update carries a marker at 75 N
    // (and one at 80 S) along its parallel, with no drift toward the equator.
    FrameFixture fixture(64, 32);
    fixture.frame.plates = { { 5 } };
    const double step = 0.25 * 2.0 * FrameFixture::kPi / 64;
    fixture.rotation({ 0.0, 0.0, step });
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(fixture, 10.0, 32.0 * 15.0 / 180.0), markerat(fixture, 40.0, 32.0 * 170.0 / 180.0) };
    std::vector<mesh::Vec3> start{ markers[0].position, markers[1].position };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    ASSERT_TRUE(tracker.attach(markers[1], fixture.view()));
    for (int update = 0; update < 40; ++update)
        tracker.advance(markers, fixture.view());
    for (int m = 0; m < 2; ++m)
    {
        EXPECT_NEAR(mesh::latitudeof(markers[m].position), mesh::latitudeof(start[m]), 1e-12) << m;
        double turned = mesh::longitudeof(markers[m].position) - mesh::longitudeof(start[m]);
        if (turned < 0) turned += 2.0 * FrameFixture::kPi;
        EXPECT_NEAR(turned, 40 * step, 1e-12) << m;
    }
}

TEST(PlateMarkers, ContinentalMarkerTransfersWithCollidingCrust)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(fixture, 3.5, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    EXPECT_TRUE(markers[0].continentalHost);

    fixture.motion(1.0);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    EXPECT_EQ(summary.rehomed, 1u);
    EXPECT_EQ(markers[0].hostUid, 2u);
    EXPECT_EQ(markers[0].state, MarkerState::attached);
}

TEST(PlateMarkers, ContinentalMarkerReturnsFromOceanicCrustToItsHost)
{
    FrameFixture fixture(8, 4);
    splitplates(fixture, 4, 1, 2);
    for (std::uint32_t y = 0; y < 4; ++y)
        for (std::uint32_t x = 4; x < 8; ++x)
            fixture.crust[fixture.cell(x, y)] = crustcode::oceanic;
    MarkerTracker tracker;
    tracker.begin(fixture.view());
    std::vector<Marker> markers{ markerat(fixture, 3.9, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));

    fixture.motion(0.3);
    const MarkerStepSummary summary = tracker.advance(markers, fixture.view());

    // Back to the centre of the nearest host cell with continental crust: (3, 1).
    EXPECT_EQ(summary.snapped, 1u);
    EXPECT_EQ(markers[0].hostUid, 1u);
    const auto snapped = fixture.coordinates(markers[0].position);
    EXPECT_NEAR(snapped[0], 3.5, 1e-9);
    EXPECT_NEAR(snapped[1], 1.5, 1e-9);
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
    std::vector<Marker> markers{ markerat(fixture, 3.5, 1.5) };
    ASSERT_TRUE(tracker.attach(markers[0], fixture.view()));
    EXPECT_FALSE(markers[0].continentalHost);
    EXPECT_DOUBLE_EQ(markers[0].hostCrustBirthMyr, 90.0);

    fixture.frame.timeMyr = 101.0;
    fixture.motion(1.0);
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
    std::vector<Marker> markers{ markerat(fixture, 3.5, 1.5) };
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
    std::vector<Marker> markers{ markerat(fixture, 1.5, 1.5) };
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
    Marker marker = markerat(fixture, 1.5, 1.5);
    EXPECT_FALSE(tracker.attach(marker, fixture.view()));
}
