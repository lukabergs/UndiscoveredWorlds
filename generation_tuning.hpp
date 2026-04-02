#pragma once

#include <array>

namespace tuning
{
namespace worlddefaults
{
inline constexpr int size = 2;
inline constexpr int type = 2;
inline constexpr int width = 2047;
inline constexpr int height = 1024;
inline constexpr bool rotation = true;
inline constexpr float tilt = 22.5f;
inline constexpr float eccentricity = 0.0167f;
inline constexpr short perihelion = 0;
inline constexpr float gravity = 1.0f;
inline constexpr float lunar = 1.0f;
inline constexpr float tempDecrease = 6.5f;
inline constexpr int northPolarAdjust = 0;
inline constexpr int southPolarAdjust = 0;
inline constexpr int averageTemp = 14;
inline constexpr float waterPickup = 1.0f;
inline constexpr float riverFactor = 15.0f;
inline constexpr int riverLandReduce = 20;
inline constexpr int estuaryLimit = 20;
inline constexpr int glacialTemp = 4;
inline constexpr int glacierTemp = -1;
inline constexpr float mountainReduce = 0.75f;
inline constexpr int climateCount = 31;
inline constexpr int maxElevation = 48000;
inline constexpr int seaLevel = 24000;
inline constexpr int craterCount = 0;
}

namespace climate
{
inline constexpr std::array<int, 10> windZoneBorders = { 28, 32, 58, 64, 86, 94, 116, 122, 148, 152 };
inline constexpr int windVariationMin = 10;
inline constexpr int windVariationMax = 30;
inline constexpr int windNodeStep = 20;
inline constexpr int windLineColour = 1000000;

namespace oceancurrents
{
inline constexpr float equatorialBand = 18.0f;
inline constexpr float midLatitudeBand = 55.0f;
inline constexpr float polarBand = 72.0f;
inline constexpr float counterCurrentBand = 4.0f;
inline constexpr float equatorialSpeed = 72.0f;
inline constexpr float counterCurrentSpeed = 18.0f;
inline constexpr float midLatitudeSpeed = 56.0f;
inline constexpr float polarSpeed = 24.0f;
inline constexpr float westernBoundarySpeed = 44.0f;
inline constexpr float easternBoundarySpeed = 34.0f;
inline constexpr float subpolarBoundarySpeed = 24.0f;
inline constexpr float seasonalShiftFactor = 0.35f;
inline constexpr float retainedBaseStrength = 0.55f;
inline constexpr float smoothingBlend = 0.35f;
inline constexpr float blockedComponentFactor = 0.25f;
inline constexpr int coastalSearchDistance = 8;
inline constexpr int smoothingIterations = 6;
}

namespace sst
{
inline constexpr float advectionBlend = 0.65f;
inline constexpr float westernBoundaryWarming = 2.4f;
inline constexpr float easternBoundaryCooling = 3.0f;
inline constexpr float evaporationScale = 5.0f;
inline constexpr float evaporationCurrentBoost = 0.04f;
inline constexpr float minimumSst = -8.0f;
inline constexpr float maximumSst = 38.0f;
inline constexpr float advectionSampleDistance = 24.0f;
inline constexpr int smoothingIterations = 2;
}

namespace pressure
{
inline constexpr float equatorialLow = 32.0f;
inline constexpr float equatorialWidth = 15.5f;
inline constexpr float subtropicalLatitude = 27.0f;
inline constexpr float subtropicalHigh = 20.0f;
inline constexpr float subtropicalWidth = 9.0f;
inline constexpr float subpolarLatitude = 57.0f;
inline constexpr float subpolarLow = 20.0f;
inline constexpr float subpolarWidth = 10.0f;
inline constexpr float polarLatitude = 82.0f;
inline constexpr float polarHigh = 16.0f;
inline constexpr float polarWidth = 8.0f;
inline constexpr float landThermalResponse = 1.08f;
inline constexpr float oceanThermalResponse = 0.7f;
inline constexpr float elevationResponse = 0.0011f;
inline constexpr float seasonalShiftFactor = 0.35f;
inline constexpr int smoothingIterations = 2;
}

namespace atmosphere
{
inline constexpr float directFlowFactor = 0.75f;
inline constexpr float geostrophicFactor = 2.3f;
inline constexpr float coriolisLatitude = 35.0f;
inline constexpr float maxVectorWind = 42.0f;
inline constexpr float blockedComponentFactor = 0.35f;
inline constexpr float scalarWindDivisor = 4.0f;
inline constexpr float minimumScalarZonalWind = 0.75f;
inline constexpr int smoothingIterations = 2;
inline constexpr int topographySmoothingIterations = 6;
inline constexpr int topographyIterations = 2;
inline constexpr float topographyMinimumWindSpeed = 1.5f;
inline constexpr float topographyMinimumRelief = 250.0f;
inline constexpr float topographyGradientScale = 1800.0f;
inline constexpr float topographyLookaheadRiseScale = 5000.0f;
inline constexpr int topographyLookaheadDistance = 6;
inline constexpr float topographySideSampleDistance = 3.0f;
inline constexpr float topographyDeflectionFactor = 0.85f;
inline constexpr float topographyChannelFactor = 0.45f;
inline constexpr float topographyDownslopeAcceleration = 0.15f;
inline constexpr float topographySpeedReduction = 0.18f;
inline constexpr float topographyVerticalMotionWindScale = 18.0f;
inline constexpr float topographyVerticalMotionStorageScale = 100.0f;
}

namespace moistureadvection
{
inline constexpr int iterations = 30;
inline constexpr float sourceScale = 0.16f;
inline constexpr float landEvaporationFactor = 0.07f;
inline constexpr float carryRetention = 0.992f;
inline constexpr float advectionDistanceScale = 1.32f;
inline constexpr int maxAdvectionDistance = 38;
inline constexpr float transportMaxFraction = 0.94f;
inline constexpr float oceanCondensation = 0.014f;
inline constexpr float landCondensation = 0.024f;
inline constexpr float calmOceanBoost = 0.004f;
inline constexpr float convergenceFactor = 0.053f;
inline constexpr float lowPressureFactor = 0.028f;
inline constexpr float highPressureFactor = 0.0010f;
inline constexpr float upliftScale = 1800.0f;
inline constexpr float upliftFactor = 0.24f;
inline constexpr float descentFactor = 0.022f;
inline constexpr float overflowFactor = 0.82f;
inline constexpr float moistureCapacityBase = 24.0f;
inline constexpr float moistureCapacityTemperatureFactor = 1.6f;
inline constexpr float rainfallScale = 2.05f;
inline constexpr float seaIceFactor = 0.3f;
inline constexpr float saturationThreshold = 0.55f;
inline constexpr float humidityCondensationFactor = 0.20f;
inline constexpr float minimumForcedHumidity = 0.22f;
inline constexpr float convergenceStorageScale = 100.0f;
}

namespace maritime
{
inline constexpr float influenceStorageScale = 1000.0f;
inline constexpr float thermalAnomalyStorageScale = 10.0f;
inline constexpr int maxSearchDistance = 48;
inline constexpr int oceanSampleCount = 8;
inline constexpr float minimumWindStrength = 1.5f;
inline constexpr float windScale = 16.0f;
inline constexpr float evaporationScale = 220.0f;
inline constexpr float maxThermalAnomaly = 20.0f;
inline constexpr float minimumRainfallModeration = 0.26f;
inline constexpr float rainfallModerationFactor = 0.90f;
inline constexpr float continentalFetchScale = 80.0f;
inline constexpr float continentalMaritimeReduction = 0.94f;
}

namespace coastalclimate
{
inline constexpr float coldSeasonThermalFactor = 0.40f;
inline constexpr float warmSeasonThermalFactor = 0.28f;
inline constexpr float rangeModerationFactor = 0.64f;
inline constexpr float fetchModerationFactor = 0.45f;
inline constexpr float minimumInfluence = 0.02f;
inline constexpr float maximumSeasonalShift = 14.0f;
}

namespace lapse
{
inline constexpr float moistLapseRate = 5.0f;
inline constexpr float dryLapseRate = 7.4f;
inline constexpr float rainfallHumidityScale = 150.0f;
inline constexpr float moistureHumidityScale = 115.0f;
inline constexpr int reliefRadius = 2;
inline constexpr float reliefScale = 1400.0f;
inline constexpr float reliefCoolingFactor = 0.8f;
inline constexpr float upliftCoolingFactor = 0.45f;
inline constexpr float subsidenceWarmingFactor = 0.65f;
inline constexpr float maximumAdditionalCooling = 2.5f;
inline constexpr float maximumLeeWarming = 3.5f;
}

namespace oceanrain
{
inline constexpr int landMultiplier = 5;
inline constexpr int fractalRange = 60;
inline constexpr int maxRain = 1500;
inline constexpr int landShadowFactor = 4;
}

namespace prevailingrain
{
inline constexpr int seaMultiplier = 60;
inline constexpr float oceanTemperatureFactor = 80.0f;
inline constexpr float minimumTemperatureMultiplier = 0.15f;
inline constexpr int dumpRate = 80;
inline constexpr int landPickupRate = 40;
inline constexpr int swerveChance = 3;
inline constexpr int spreadChance = 2;
inline constexpr float newSeedProportion = 0.95f;
inline constexpr float horseSeedProportion = 0.4f;
inline constexpr int splashSize = 1;
inline constexpr float elevationFactor = 0.002f;
inline constexpr int slopeMinimum = 300;
inline constexpr float seasonalVariationTiltDivisor = 3750.0f;
inline constexpr float tropicalSeasonalVariationTiltDivisor = 2250.0f;
inline constexpr int maxSeasonalDistance = 50;
inline constexpr float heatPickupRate = 0.008f;
inline constexpr float heatDepositRate = 0.15f;
inline constexpr float summerFactorTiltDivisor = 75.0f;
inline constexpr float maxWinterHeatFactor = 3.5f;
inline constexpr float winterIceFactor = 0.05f;
inline constexpr float tropicalRainReduction = 0.4f;
inline constexpr float slopeBase = 160.0f;
}

namespace lakerain
{
inline constexpr int lakeMultiplier = 15;
inline constexpr float temperatureFactor = 20.0f;
inline constexpr float minimumTemperatureMultiplier = 0.15f;
inline constexpr int dumpRate = 5;
inline constexpr int pickupRate = 60;
inline constexpr int swerveChance = 3;
inline constexpr int splashSize = 1;
inline constexpr float slopeFactor = 5.0f;
inline constexpr float elevationFactor = 0.002f;
inline constexpr int slopeMinimum = 200;
inline constexpr float seasonalVariation = 0.02f;
inline constexpr int maxRain = 800;
inline constexpr float capFactor = 0.1f;
}

namespace riftlakerain
{
inline constexpr int lakeMultiplier = 15;
inline constexpr float temperatureFactor = 20.0f;
inline constexpr float minimumTemperatureMultiplier = 0.15f;
inline constexpr int dumpRate = 5;
inline constexpr int pickupRate = 300;
inline constexpr int swerveChance = 3;
inline constexpr int splashSize = 1;
inline constexpr int slopeFactor = 50;
inline constexpr float elevationFactor = 0.002f;
inline constexpr int slopeMinimum = 200;
inline constexpr float seasonalVariation = 0.02f;
inline constexpr int maxRain = 800;
inline constexpr float capFactor = 0.1f;
}

namespace desertworld
{
inline constexpr float slopeFactor = 130.0f;
inline constexpr float idealTemperature = 20.0f;
inline constexpr float maxTemperatureDifference = 40.0f;
inline constexpr int maxRainMin = 20;
inline constexpr int maxRainMax = 200;
inline constexpr int fractalGrain = 8;
inline constexpr float fractalValueMod = 0.2f;
inline constexpr int fractalValueMod2Min = 1;
inline constexpr int fractalValueMod2Max = 4;
inline constexpr int warpFactor = 60;
inline constexpr int blurDistance = 1;
}

namespace monsoon
{
inline constexpr float strengthCenterTilt = 32.5f;
inline constexpr float strengthTiltDivisor = 10.0f;
inline constexpr int minimumTemperatureDifference = 2;
inline constexpr int minimumAverageTemperature = 15;
inline constexpr float temperatureDifferenceFactor = 1.0f;
inline constexpr float temperatureFactor = 400.0f;
inline constexpr float inlandTemperatureFactor = 25.0f;
inline constexpr float minimumTideFactor = 0.2f;
inline constexpr float initialIncrease = 1.8f;
inline constexpr float increaseDecreasePerTick = 0.015f;
inline constexpr float minimumIncrease = 0.99f;
inline constexpr int dumpRate = 15;
inline constexpr float elevationFactor = 0.01f;
inline constexpr int slopeMinimum = 200;
inline constexpr int swerveChance = 2;
inline constexpr float maxSummerRain = 410.0f;
inline constexpr int equatorDistance = 30;
inline constexpr float equatorFactor = 0.25f;
inline constexpr int minimumEquatorDistance = 10;
inline constexpr int maxTicks = 1500;
inline constexpr int slopeBase = 77;
inline constexpr float elevationRainReductionFactor = 8.0f;
}

namespace mediterranean
{
inline constexpr float strengthCenterTilt = 31.5f;
inline constexpr float strengthTiltDivisor = 10.0f;
inline constexpr int targetMaxTemperature = 19;
inline constexpr float maxTemperatureDifferenceFactor = 0.001f;
inline constexpr float minimumColdTemperature = 2.0f;
inline constexpr float minimumColdTemperatureFactor = 0.001f;
inline constexpr int maxRain = 500;
inline constexpr float maxRainDifferenceFactor = 0.001f;
inline constexpr int maxInlandDistance = 40;
inline constexpr float maxInlandDifferenceFactor = 0.01f;
inline constexpr float horseLatitudeDifferenceFactor = 0.02f;
inline constexpr int smoothDistance = 3;
}

namespace equatorialrain
{
inline constexpr float strengthCenterTilt = 32.5f;
inline constexpr float strengthTiltDivisor = 10.0f;
inline constexpr float winterAdditionFactor = 0.3f;
inline constexpr float summerAdditionFactor = 0.4f;
}

namespace temperaturerainfall
{
inline constexpr float optimumTilt = 22.5f;
inline constexpr float tiltStrengthDivisor = 10.0f;
inline constexpr float optimumAverageTemperature = 14.0f;
inline constexpr float averageTemperatureStrengthDivisor = 20.0f;
inline constexpr float winterRainWarmth = 0.10f;
inline constexpr float maxWinterRainWarmth = 10.0f;
inline constexpr float summerRainCooling = 0.0025f;
inline constexpr float noRainSummerHeatMultiplier = 1.3f;
inline constexpr float noRainWinterColdMultiplier = 1.08f;
inline constexpr float maxWinterVariation = 15.0f;
inline constexpr float maxSummerVariation = 10.0f;
inline constexpr float maxAffectedTemperature = 40.0f;
inline constexpr float minAffectedTemperature = -20.0f;
inline constexpr float inlandOffset = 20.0f;
inline constexpr float inlandFactorNumerator = 5.0f;
inline constexpr float changeMultiplier = 100.0f;
}

namespace continentality
{
inline constexpr float maxAffectedTemperature = 30.0f;
inline constexpr float minAffectedTemperature = -10.0f;
inline constexpr float winterRemovalFactor = 0.6f;
inline constexpr float summerAdditionFactor = 0.05f;
inline constexpr float strengthTiltBaseline = 22.5f;
inline constexpr float maxStrength = 1.5f;
inline constexpr float maxWinterRemoval = 20.0f;
inline constexpr float maxSummerAddition = 10.0f;
inline constexpr float tundraTweakTiltBaseline = 22.5f;
inline constexpr float tundraTweakTiltDivisor = 30.0f;
inline constexpr int tundraMaxTemperature = 12;
inline constexpr float tundraWarmFactor = 10.0f;
inline constexpr float tundraElevationDivisor = 400.0f;
inline constexpr int tundraClampedMaxTemperature = 13;
}

namespace saltlakes
{
inline constexpr int maxRain = 30;
inline constexpr int minimumTemperature = 15;
inline constexpr int minimumRiverFlow = 200;
inline constexpr int lakeChance = 100;
inline constexpr int depressionChance = 100;
inline constexpr int edgeMargin = 20;
inline constexpr int largeShapeChance = 8;
inline constexpr int largeShapeMin = 5;
inline constexpr int largeShapeMax = 11;
inline constexpr int smallShapeMin = 2;
inline constexpr int smallShapeMax = 5;
inline constexpr int depthMin = 5;
inline constexpr int depthMax = 30;
}

namespace rivers
{
inline constexpr int mountainHeightLimit = 600;
inline constexpr int minimumFlow = 40;
inline constexpr int maxRepeatDirection = 2;
}
}

namespace terrain
{
namespace shared
{
inline constexpr int commonFractalValueMod2Min = 3;
inline constexpr int commonFractalValueMod2Max = 6;
}

namespace oceanridges
{
inline constexpr int boundaryMaxSourceDifference = 400;
inline constexpr int gridSize = 16;
inline constexpr int pointShiftFractalGrain = 8;
inline constexpr float pointShiftFractalValueMod = 0.2f;
inline constexpr int pointShiftFractalValueMod2Min = 3;
inline constexpr int pointShiftFractalValueMod2Max = 6;
inline constexpr int maxShift = 40;
inline constexpr int maxAdditionalShift = 6;
inline constexpr int minAdditionalShift = 1;
inline constexpr int firstPassMaxRadius = 70;
inline constexpr int firstPassHeightMultiplier = 6;
inline constexpr int maxVolcanoRadius = 6;
inline constexpr int secondPassMaxRadius = 50;
inline constexpr int secondPassHeightMultiplier = 6;
inline constexpr short ridgeAngleSearchDistance = 3;
inline constexpr int regionalDisplacement = 16;
inline constexpr int regionalDisplacementFractalGrain = 128;
inline constexpr float regionalDisplacementValueMod = 4.0f;
inline constexpr float regionalDisplacementValueMod2 = 8.0f;
inline constexpr int faultPasses = 4;
inline constexpr int faultStep = 8;
inline constexpr int faultVariation = 3;
inline constexpr int faultLookDistance = 4;
}

namespace platetectonics
{
inline constexpr float seaLevelBias = 0.0f;
inline constexpr uint32_t erosionPeriod = 0;
inline constexpr float foldingRatio = 0.02f;
inline constexpr uint32_t aggregationOverlapAbsolute = 1000000;
inline constexpr float aggregationOverlapRelative = 0.33f;
inline constexpr uint32_t cycleCount = 4;
inline constexpr uint32_t plateCount = 4;
inline constexpr uint32_t maximumSimulationSteps = 20000;
inline constexpr int minimumOceanDepth = 1;
inline constexpr int coastalOceanOffset = 10;
inline constexpr int landStartOffset = 50;
inline constexpr float oceanExponent = 1.15f;
inline constexpr float landExponent = 0.85f;
inline constexpr float outputBlend = 1.0f;
inline constexpr float landRetentionSeaBias = 0.08f;
inline constexpr float convergentBoundaryThreshold = 0.18f;
inline constexpr int collisionUplift = 1800;
inline constexpr int collisionMinimumPeak = 4500;
inline constexpr int collisionMaximumPeak = 18000;
}

namespace edgeseams
{
inline constexpr int bandSize = 36;
inline constexpr int seaDepthBase = 500;
inline constexpr int seaDepthVariation = 900;
inline constexpr float edgeExponent = 1.8f;
inline constexpr float bandNoiseMinScale = 0.6f;
inline constexpr float bandNoiseRange = 0.8f;
inline constexpr float distanceJitter = 0.9f;
}

namespace fastlem
{
inline constexpr int cellSize = 6;
inline constexpr int baseLatticeStep = 5;
inline constexpr int outletMinimumDistanceBlocks = 18;
inline constexpr float minimumCandidateSpacingBlocks = 3.0f;
inline constexpr float siteWeightMin = 0.80f;
inline constexpr float siteWeightMax = 1.40f;
inline constexpr float clusterNoiseThreshold = 0.62f;
inline constexpr int minimumLandTilesPerSite = 2;
inline constexpr int fallbackNeighbourRadius = 3;
inline constexpr int minimumConnections = 3;
inline constexpr int iterations = 24;
inline constexpr float mExponent = 0.5f;
inline constexpr float baseUplift = 0.60f;
inline constexpr float inlandUplift = 6.20f;
inline constexpr float reliefUplift = 1.80f;
inline constexpr float noiseUplift = 1.60f;
inline constexpr float baseErodibility = 0.32f;
inline constexpr float coastalErodibility = 0.22f;
inline constexpr float noiseErodibility = 0.10f;
inline constexpr float minimumErodibility = 0.08f;
inline constexpr float maximumErodibility = 2.5f;
inline constexpr float maxSlopeRadians = 0.38f;
inline constexpr int minimumRidgeCoastDistance = 2;
inline constexpr int minimumPeakCoastDistance = 4;
inline constexpr float minimumRidgeElevationNormalised = 0.14f;
inline constexpr float minimumPeakElevationNormalised = 0.32f;
inline constexpr float minimumRidgeScore = 0.22f;
inline constexpr float elevationScoreWeight = 0.60f;
inline constexpr float coastScoreWeight = 0.40f;
inline constexpr float peakHeightExponent = 0.95f;
inline constexpr int minimumPeakHeight = 1400;
inline constexpr int maximumPeakHeight = 10500;
inline constexpr int minimumCandidateSites = 18;
inline constexpr int minimumBoundaryCandidates = 8;
inline constexpr int postRiverMaximumDistance = 18;
inline constexpr int postRiverMountainDistance = 16;
inline constexpr int postRiverMaximumUplift = 650;
inline constexpr int postRiverMinimumUplift = 24;
inline constexpr int postRiverMinimumMountainHeight = 500;
inline constexpr float postRiverDistanceExponent = 1.25f;
inline constexpr float postRiverMountainExponent = 0.90f;
inline constexpr int postRiverElevationRange = 5000;
}
}

namespace regional
{
inline constexpr float submarineWarpFactor = 100.0f;
inline constexpr int submarineExtraMargin = 20;
inline constexpr int ridgeRadiationExtraMargin = 20;
inline constexpr int saltPanMaxRain = 30;
inline constexpr float mountainPrecipitationSlopeFactor = 400.0f;
}
}
