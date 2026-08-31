#include "climate_atmosphere.hpp"

#include <cmath>
#include <iostream>

namespace
{
int failures = 0;

void expect(bool condition, const char* message)
{
    if (condition)
        return;

    std::cerr << "FAIL: " << message << '\n';
    failures++;
}
}

int main()
{
    constexpr float earthRotation = 7.2921159e-5f;
    expect(
        std::abs(climateatmosphere::coriolisParameterPerSecond(0.0f, earthRotation)) < 1.0e-9f,
        "Coriolis acceleration must vanish at the equator");
    expect(
        climateatmosphere::coriolisParameterPerSecond(45.0f, earthRotation) > 0.0f &&
            climateatmosphere::coriolisParameterPerSecond(-45.0f, earthRotation) < 0.0f,
        "Coriolis sign must reverse across the equator");

    const float heightResponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        100000.0f, 50000.0f);
    expect(
        std::abs(heightResponse - 20.27f) < 0.05f,
        "1000-to-500 hPa thickness must change by about 20.27 metres per kelvin");

    const float earthHadleyEdge = climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
        60.0f, 10000.0f, 288.0f, 9.80665f, earthRotation, 6371000.0f);
    expect(
        earthHadleyEdge > 20.0f && earthHadleyEdge < 26.0f,
        "Held-Hou scaling must place the Earth-like Hadley edge in the subtropics");
    expect(
        climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
            60.0f, 10000.0f, 288.0f, 9.80665f, 0.0f, 6371000.0f) == 90.0f,
        "a non-rotating atmosphere must permit pole-to-pole overturning");

    const float warmSurfacePressure = climateatmosphere::thermalSurfacePressureAnomalyHpa(
        10.0f, 1000.0f, 288.0f, 0.12f);
    const float coldSurfacePressure = climateatmosphere::thermalSurfacePressureAnomalyHpa(
        -10.0f, 1000.0f, 288.0f, 0.12f);
    expect(
        warmSurfacePressure < 0.0f && coldSurfacePressure > 0.0f &&
            std::abs(warmSurfacePressure + coldSurfacePressure) < 0.0001f,
        "warm columns must form thermal lows and cold columns thermal highs");

    const float equatorialPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(0.0f, 0.0f, 25.0f, 10.0f);
    const float subtropicalPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(25.0f, 0.0f, 25.0f, 10.0f);
    const float subpolarPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(57.5f, 0.0f, 25.0f, 10.0f);
    const float polarPressure =
        climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(90.0f, 0.0f, 25.0f, 10.0f);
    expect(
        equatorialPressure < 0.0f && subtropicalPressure > 0.0f &&
            subpolarPressure < 0.0f && polarPressure > 0.0f,
        "axisymmetric overturning must produce equatorial and subpolar lows with subtropical and polar highs");

    constexpr int responseColumns = 72;
    constexpr int responseRows = 37;
    const auto responseIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * responseColumns + x);
    };
    std::vector<float> thermalImpulse(
        static_cast<size_t>(responseColumns * responseRows),
        0.0f);
    const int impulseColumn = responseColumns / 2;
    const int equatorRow = responseRows / 2;
    thermalImpulse[responseIndex(impulseColumn, equatorRow)] = 1.0f;
    const auto progradeResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        1.0f);
    const auto retrogradeResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        -1.0f);
    expect(
        progradeResponse[responseIndex(impulseColumn - 1, equatorRow)] >
            progradeResponse[responseIndex(impulseColumn + 1, equatorRow)] &&
        retrogradeResponse[responseIndex(impulseColumn + 1, equatorRow)] >
            retrogradeResponse[responseIndex(impulseColumn - 1, equatorRow)],
        "the tropical response must extend westward relative to planetary rotation");
    const int extratropicalRow = 10;
    thermalImpulse.assign(static_cast<size_t>(responseColumns * responseRows), 0.0f);
    thermalImpulse[responseIndex(impulseColumn, extratropicalRow)] = 1.0f;
    const auto extratropicalResponse = climateatmosphere::nonlocalThermalResponse(
        responseColumns,
        responseRows,
        thermalImpulse,
        20.0f,
        20.0f,
        8.0f,
        1.0f);
    expect(
        std::abs(
            extratropicalResponse[responseIndex(impulseColumn - 1, extratropicalRow)] -
            extratropicalResponse[responseIndex(impulseColumn + 1, extratropicalRow)]) < 1.0e-6f,
        "the extratropical thermal response must be zonally symmetric");
    double equatorialResponseSum = 0.0;
    for (int x = 0; x < responseColumns; x++)
        equatorialResponseSum += progradeResponse[responseIndex(x, equatorRow)];
    expect(
        equatorialResponseSum > 0.0 && equatorialResponseSum < 1.0,
        "the nonlocal response must redistribute an impulse without creating a new extremum");

    constexpr int topographicColumns = 21;
    constexpr int topographicRows = 9;
    const auto topographicIndex = [=](int x, int y)
    {
        return static_cast<size_t>(y * topographicColumns + x);
    };
    std::vector<float> ridgeTerrain(
        static_cast<size_t>(topographicColumns * topographicRows),
        0.0f);
    std::vector<float> ridgeEastWind(ridgeTerrain.size(), 10.0f);
    std::vector<float> ridgeSouthWind(ridgeTerrain.size(), 0.0f);
    for (int y = 0; y < topographicRows; y++)
        ridgeTerrain[topographicIndex(10, y)] = 2000.0f;
    const auto ridgeForcing = climateatmosphere::mechanicalTopographicPressureForcingHpa(
        topographicColumns,
        topographicRows,
        ridgeTerrain,
        ridgeEastWind,
        ridgeSouthWind,
        3.0f,
        2000.0f,
        2.0f,
        1.0f,
        8.0f,
        10.0f,
        30.0f);
    expect(
        ridgeForcing[topographicIndex(7, 2)] > 0.0f &&
            ridgeForcing[topographicIndex(13, 2)] < 0.0f,
        "westerly flow across a ridge must produce a windward high and lee trough");
    ridgeEastWind.assign(ridgeTerrain.size(), -10.0f);
    const auto reversedRidgeForcing =
        climateatmosphere::mechanicalTopographicPressureForcingHpa(
            topographicColumns,
            topographicRows,
            ridgeTerrain,
            ridgeEastWind,
            ridgeSouthWind,
            3.0f,
            2000.0f,
            2.0f,
            1.0f,
            8.0f,
            10.0f,
            30.0f);
    expect(
        reversedRidgeForcing[topographicIndex(13, 2)] > 0.0f &&
            reversedRidgeForcing[topographicIndex(7, 2)] < 0.0f,
        "topographic pressure forcing must reverse when the background wind reverses");

    const auto spacing = climateatmosphere::cellSpacingMetres(0.0f, 2048, 1025, 6371000.0f);
    expect(
        spacing.zonalMetres > 19000.0f && spacing.zonalMetres < 20000.0f &&
            spacing.meridionalMetres > 19000.0f && spacing.meridionalMetres < 20000.0f,
        "Earth benchmark cells must be about 19.5 kilometres at the equator");

    const auto unrotated = climateatmosphere::steadyRayleighCoriolisWind(
        0.001f, 0.002f, 0.0f, 1000.0f, earthRotation);
    expect(
        std::abs(unrotated.eastMetresPerSecond - 1.0f) < 0.001f &&
            std::abs(unrotated.southMetresPerSecond + 2.0f) < 0.001f,
        "Rayleigh flow at the equator must follow the pressure-gradient acceleration");

    const auto northern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, 45.0f, 86400.0f, earthRotation);
    const auto southern = climateatmosphere::steadyRayleighCoriolisWind(
        0.0f, 0.001f, -45.0f, 86400.0f, earthRotation);
    expect(
        northern.eastMetresPerSecond > 0.0f && southern.eastMetresPerSecond < 0.0f,
        "the same meridional height force must produce opposite zonal flow across the equator");

    return failures == 0 ? 0 : 1;
}
