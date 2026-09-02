#pragma once

#include <array>

namespace tuning
{
namespace worlddefaults
{
inline constexpr int size = 2;
inline constexpr int width = 599;
inline constexpr int height = 399;
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

namespace climateresolution
{
inline constexpr float gridScale(int width, int height)
{
    const float widthscale = static_cast<float>(width + 1) / static_cast<float>(worlddefaults::width + 1);
    const float heightscale = static_cast<float>(height + 1) / static_cast<float>(worlddefaults::height + 1);
    return (widthscale + heightscale) / 2.0f;
}

inline constexpr int scaleDistance(int distance, int width, int height)
{
    const int scaled = static_cast<int>(static_cast<float>(distance) * gridScale(width, height) + 0.5f);
    return scaled > 0 ? scaled : 1;
}

inline constexpr float scaleDistance(float distance, int width, int height)
{
    return distance * gridScale(width, height);
}

static_assert(scaleDistance(30, worlddefaults::width, worlddefaults::height) == 30);
static_assert(scaleDistance(30, 2047, 1024) == 90);
}

namespace climate
{
inline constexpr std::array<int, 10> windZoneBorders = { 28, 32, 58, 64, 86, 94, 116, 122, 148, 152 };
inline constexpr int windVariationMin = 10;
inline constexpr int windVariationMax = 30;
inline constexpr int windNodeStep = 20;
inline constexpr int windLineColour = 1000000;

namespace energybalance
{
inline constexpr double earthLongwaveInterceptWm2 = 210.0;
inline constexpr double longwaveSlopeWm2K = 2.0;
inline constexpr double meridionalTransportWm2K = 3.0;
inline constexpr double zonalLandOceanExchangeWm2K = 12.0;
inline constexpr double landHeatCapacityJm2K = 7.5e6;
inline constexpr double landDeepLayerHeatCapacityJm2K = 3.0e7;
inline constexpr double landDeepLayerCouplingWm2K = 1.5;
inline constexpr double oceanMixedLayerHeatCapacityJm2K = 1.25e8;
inline constexpr double landAlbedo = 0.20;
inline constexpr double oceanAlbedo = 0.27;
inline constexpr double snowAlbedo = 0.55;
inline constexpr double permanentIceAlbedo = 0.80;
inline constexpr double permanentIceWarmTransitionC = -4.0;
inline constexpr double permanentIceColdTransitionC = -14.0;
inline constexpr double iceSheetSeedMinimumElevationKm = 0.75;
inline constexpr double iceSheetSeedFullElevationKm = 1.75;
inline constexpr int localIceCouplingMaximumIterations = 80;
inline constexpr double localIceCouplingRelaxation = 0.25;
inline constexpr double localIceCouplingTolerance = 0.005;
inline constexpr int localIceCouplingStableIterations = 5;
inline constexpr double temperatureCalibrationRelaxation = 0.60;
inline constexpr double temperatureCalibrationToleranceC = 0.002;
inline constexpr double seaIceAlbedo = 0.58;
inline constexpr double fractalTemperatureRangeC = 4.0;
}

namespace oceancurrents
{
inline constexpr int internalHorizontalCells = 64;
inline constexpr float currentStorageCentimetresPerSecond = 100.0f;
inline constexpr float waterDensityKgM3 = 1025.0f;
inline constexpr float mixedLayerDepthMetres = 60.0f;
inline constexpr float barotropicDragPerSecond = 1.5e-6f;
inline constexpr float heatDiffusivityM2S = 750.0f;
inline constexpr float surfaceHeatExchangeWm2K = 18.0f;
inline constexpr int streamfunctionIterations = 800;
inline constexpr int heatStepsPerCouplingIteration = 30;
inline constexpr int couplingIterations = 40;
inline constexpr float couplingUnderRelaxation = 0.35f;
inline constexpr float couplingTolerance = 1.0e-3f;
inline constexpr float sstWindFeedbackMpsPerK = 0.08f;
inline constexpr float maximumCurrentMps = 2.5f;
inline constexpr bool oneWayDiagnosticsOnly = false;
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

namespace atmosphere
{
inline constexpr float referencePlanetRadiusMetres = 6371000.0f;
inline constexpr float rotationRatePerSecond = 7.2921159e-5f;
inline constexpr float gravityMetresPerSecondSquared = 9.80665f;
inline constexpr float surfaceAirDensityKgM3 = 1.225f;
inline constexpr float pressurePascalsPerHectopascal = 100.0f;
inline constexpr float surfaceBoundaryLayerMomentumDepthMetres = 300.0f;
inline constexpr float oceanMomentumDragCoefficient = 0.0013f;
inline constexpr float landMomentumDragCoefficient = 0.0030f;
inline constexpr float highReliefMomentumDragCoefficient = 0.0060f;
inline constexpr float maxVectorWind = 42.0f;
inline constexpr float blockedComponentFactor = 0.35f;
inline constexpr float scalarWindDivisor = 4.0f;
inline constexpr float minimumScalarZonalWind = 0.75f;
inline constexpr int smoothingIterations = 2;
inline constexpr int topographySmoothingIterations = 6;
inline constexpr int topographySmoothingReferenceHorizontalCells = 512;
inline constexpr int maximumTopographySmoothingIterations = 48;
inline constexpr int landmaskSmoothingIterations = 8;
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
inline constexpr float maximumOrographicParcelDisplacementMetres = 7000.0f;
inline constexpr float topographyVerticalMotionStorageScale = 4.0f;
}

namespace circulation
{
inline constexpr int iterations = 24;
inline constexpr float surfaceReferencePressurePa = 100000.0f;
inline constexpr float upperReferencePressurePa = 50000.0f;
// Thermal anomalies taper to zero at the top of the trade-cumulus layer.
// This thermodynamic depth is separate from the surface momentum depth.
inline constexpr float surfaceThermalTopPressurePa = 70000.0f;
inline constexpr float referenceTemperatureK = 288.0f;
inline constexpr float troposphereHeightMetres = 16000.0f; // Held-Hou overturning depth
inline constexpr float overturningMassRedistributionEfficiency = 0.060f;
inline constexpr float mechanicalTopographicPressureAmplitudeHpa = 6.0f;
inline constexpr float mechanicalTopographicTerrainScaleMetres = 2000.0f;
inline constexpr float mechanicalTopographicSampleDistanceAtReferenceGrid = 4.0f;
inline constexpr float mechanicalTopographicMinimumWindMps = 1.5f;
inline constexpr float mechanicalTopographicFullStrengthWindMps = 8.0f;
inline constexpr float mechanicalTopographicMinimumLatitudeDegrees = 10.0f;
inline constexpr float mechanicalTopographicFullStrengthLatitudeDegrees = 30.0f;
inline constexpr int mechanicalTopographicSmoothingAtReferenceGrid = 4;
inline constexpr float minimumOverturningPressureAmplitudeHpa = 4.0f;
inline constexpr float maximumOverturningPressureAmplitudeHpa = 18.0f;
inline constexpr float minimumHadleyHalfWidthDegrees = 10.0f;
inline constexpr float maximumHadleyHalfWidthDegrees = 60.0f;
inline constexpr float maximumThermalEquatorShiftHadleyFraction = 0.65f;
inline constexpr int zonalTemperatureSmoothingIterations = 4;
inline constexpr float upperLayerDragTimeSeconds = 86400.0f;
inline constexpr int thermalHeightSmoothingIterations = 10;
inline constexpr int divergenceSmoothingIterations = 6;
inline constexpr float secondsPerDay = 86400.0f;
inline constexpr float timeStepDays = 0.125f;
inline constexpr float layerPressureDepthHpa = 500.0f;
inline constexpr float surfacePressureReferenceHpa = 1000.0f;
inline constexpr float surfacePressureRestoringTimeDays = 6.0f;
inline constexpr float divergentSurfacePressureTendencyResponse = 0.0f;
inline constexpr float maximumSurfacePressureTendencyHpaPerDay = 12.0f;
inline constexpr float maximumSurfacePressureAnomalyHpa = 35.0f;
inline constexpr int stationaryWaveLongitudeCells = 64;
inline constexpr float stationaryWaveLinearizationWindMps = 5.0f;
inline constexpr float surfaceModeBruntVaisalaPerSecond = 0.012f;
inline constexpr float surfaceModeDepthMetres = 10000.0f;
inline constexpr float upperModeBruntVaisalaPerSecond = 0.018f;
inline constexpr float upperModeDepthMetres = 7000.0f;
inline constexpr float stationaryWaveNondimensionalDamping = 0.10f;
inline constexpr bool stationaryWavePreserveZonalMean = true;
inline constexpr int stationaryWaveMaximumIterations = 2000;
inline constexpr float stationaryWaveRelativeTolerance = 1.0e-4f;
inline constexpr int maximumResolvedStationaryZonalWavenumber = 8;
inline constexpr int maximumResolvedStationaryMeridionalWavenumber = 24;
inline constexpr float diagnosedForcingScaleMetres = 1200000.0f;
// Legacy scalar projection/response controls are retained for experiments;
// production uses column exchanges and mass/cp/hypsometric mode conversion.
inline constexpr float diabaticVerticalProjection = 0.35f;
inline constexpr float sensibleHeatingCoefficientWm2K = 5.0f;
inline constexpr float stationaryLatentHeatProjectionFraction = 1.0f; // active latent-heating isolation control
inline constexpr float surfaceHeatingPressureResponseHpaPerWm2 = -0.04f;
inline constexpr float upperHeatingHeightResponseMetresPerWm2 = 1.5f;
// Active grey-column radiation and climate-coupling controls.
inline constexpr float lowerLongwaveOpticalDepth = 1.0f;
inline constexpr float upperLongwaveOpticalDepth = 0.5f;
inline constexpr float lowerShortwaveOpticalDepth = 0.12f;
inline constexpr float upperShortwaveOpticalDepth = 0.10f;
inline constexpr int maximumCouplingIterations = 24;
inline constexpr int minimumCouplingIterations = 2;
inline constexpr float couplingRelativeTolerance = 0.02f;
inline constexpr float couplingHeatingRelaxation = 0.5f;
inline constexpr float couplingMinimumHeatingRelaxation = 0.1f;
inline constexpr float surfaceToUpperModeCoupling = 0.10f;
inline constexpr float upperOrographicProjection = 0.20f;
inline constexpr bool enableZonalMode = true;
inline constexpr bool enableStationaryMode = true;
inline constexpr bool enableSurfaceMode = true;
inline constexpr bool enableUpperMode = true;
inline constexpr bool enableLaggedDiabaticCoupling = true;
inline constexpr bool enablePrognosticClimateAtmosphere = true;
inline constexpr int prognosticClimateSpinupSteps = 12;
inline constexpr int prognosticClimateSamples = 30;
inline constexpr float prognosticClimateStepSeconds = 21600.0f;
inline constexpr float prognosticClimateBlend = 0.30f;
inline constexpr float unresolvedCirculationGradientRetention = 0.25f;
inline constexpr float upperHeightRelaxation = 0.22f;
inline constexpr float upperHeightDiffusion = 0.16f;
inline constexpr int windSmoothingIterations = 2;
inline constexpr float verticalRelaxation = 0.24f;
inline constexpr float verticalDivergenceResponse = 0.32f;
inline constexpr float surfacePressureDiffusion = 0.30f;
inline constexpr float maximumVerticalVelocity = 100.0f;
inline constexpr float verticalVelocityStorageScale = 100.0f;
}

namespace moistureadvection
{
inline constexpr int iterations = 30;
inline constexpr int internalClimateHorizontalCells = 128;
inline constexpr int weatherPhaseCount = 3;
inline constexpr bool enableEvolvingWeatherAnomalies = true;
inline constexpr float maximumClimateStepDays = 1.0f;
inline constexpr bool enablePrescribedSynopticPerturbations = false;
inline constexpr bool enablePrescribedCoastalDiurnalCycle = true;
inline constexpr int weatherAnomalyHorizontalCells = 32;
inline constexpr float weatherAnomalyEquivalentDepthMetres = 100.0f;
inline constexpr float weatherAnomalyDragTimeDays = 3.0f;
inline constexpr float weatherAnomalyRelaxationTimeDays = 10.0f;
inline constexpr float weatherAnomalyStochasticHeightMps = 2.0e-5f;
inline constexpr float maximumWeatherAnomalyWindMps = 6.0f;
inline constexpr int minimumSpinupCycles = 2;
inline constexpr int maximumSpinupCycles = 8;
inline constexpr float spinupRelativeStorageTolerance = 0.012f;
inline constexpr float landSoilMoistureCapacity = 80.0f;
inline constexpr float initialSoilMoistureFraction = 0.5f;
inline constexpr float landInfiltrationFraction = 1.0f;
inline constexpr float advectionTimeStepSeconds = 86400.0f;
inline constexpr float surfaceExchangeCoefficient = 0.0013f;
inline constexpr float minimumSurfaceWind = 1.0f;
inline constexpr float landSurfaceResistance = 1.0f;
inline constexpr float snowSublimationResistance = 0.50f;
inline constexpr float exchangeReferenceHeightMetres = 10.0f;
inline constexpr float oceanRoughnessLengthMetres = 0.0002f;
inline constexpr float landRoughnessLengthMetres = 0.05f;
inline constexpr float iceRoughnessLengthMetres = 0.001f;
inline constexpr float minimumStabilityExchangeMultiplier = 0.25f;
inline constexpr float maximumStabilityExchangeMultiplier = 1.8f;
inline constexpr float soilMoistureCriticalFraction = 0.5f;
inline constexpr float soilMoistureStressExponent = 0.5f;
inline constexpr float upperWindTransportFraction = 0.0f;
inline constexpr float boundaryLayerCapacityFraction = 0.45f;
inline constexpr float freeTroposphereCapacityFraction = 0.55f;
inline constexpr int freeTroposphereTransportCadence = 1;
inline constexpr float backgroundVerticalExchangeTimeDays = 6.0f;
inline constexpr float ascentVerticalExchangeFraction = 0.35f;
inline constexpr float subsidenceVerticalExchangeFraction = 0.20f;
inline constexpr float convectiveVerticalExchangeFraction = 0.25f;
inline constexpr float maximumVerticalExchangeFraction = 0.65f;
inline constexpr float transientEddyMixingLengthMetres = 25000.0f;
inline constexpr float backgroundMoistureDiffusivityM2S = 10000.0f;
inline constexpr float maximumTransientEddyDiffusivityM2S = 250000.0f;
inline constexpr float transientEddyMinimumLatitudeDegrees = 0.0f;
inline constexpr float transientEddyFullStrengthLatitudeDegrees = 0.0f;
inline constexpr float referencePlanetRadiusMetres = 6371000.0f;
inline constexpr float transportMaximumCourant = 0.80f;
inline constexpr int transportCorrectivePasses = 1;
inline constexpr float freeTroposphereUpperWindDepartureFraction = 0.65f;
inline constexpr float landBoundaryLayerTransportWindFraction = 0.85f;
inline constexpr float polarUpperWindBlendTaperStartDegrees = 60.0f;
inline constexpr float polarUpperWindBlendTaperEndDegrees = 75.0f;
inline constexpr float polarMeridionalTransportTaperStartDegrees = 60.0f;
inline constexpr float polarMeridionalTransportTaperEndDegrees = 72.0f;
inline constexpr float polarConvectionTaperStartDegrees = 55.0f;
inline constexpr float polarConvectionTaperEndDegrees = 75.0f;
inline constexpr float convectiveConvergenceMixingFraction = 0.50f;
inline constexpr int convectiveConvergenceSmoothingPasses = 2;
inline constexpr float vapourHighWavenumberDamping = 0.08f;
inline constexpr float convergenceHighWavenumberDamping = 0.12f;
inline constexpr int highWavenumberDampingPasses = 1;
inline constexpr float cloudOnsetRelativeHumidity = 0.85f;
inline constexpr float stratiformCriticalRelativeHumidity = 0.98f;
inline constexpr float condensationConversionTimeDays = 2.0f;
inline constexpr float convectiveConversionEfficiency = 0.75f;
inline constexpr float kuoCriticalRelativeHumidity = 0.0f;
inline constexpr float kuoHumidityExponent = 3.0f;
inline constexpr float elevatedMoistureAccessionFraction = 0.50f;
inline constexpr float freeTroposphereEnvironmentalLapseC = 12.0f;
inline constexpr float convectiveActivationBuoyancyC = 0.0f;
inline constexpr float convectiveFullStrengthBuoyancyC = 6.0f;
inline constexpr float shallowConvectionMixingTimeDays = 2.0f;
inline constexpr float maximumShallowExchangeFraction = 0.30f;
inline constexpr float dryConvectionMixingTimeDays = 1.0f;
inline constexpr float maximumDryExchangeFraction = 0.25f;
inline constexpr float dryConvectionActivationBuoyancyC = 10.0f;
inline constexpr float dryConvectionFullStrengthBuoyancyC = 18.0f;
inline constexpr float shallowConvectionHumidityOnset = 0.55f;
inline constexpr float shallowConvectionFullHumidity = 0.85f;
inline constexpr float shallowConvectionFullShearMps = 20.0f;
inline constexpr float subgridOrographicExtremeRetention = 0.35f;
inline constexpr float synopticWindRotationDegrees = 7.0f;
inline constexpr float synopticWindPerturbationMps = 1.5f;
inline constexpr float daytimeOnshoreCoastalWindMps = 2.5f;
inline constexpr float nighttimeOffshoreCoastalWindMps = 1.8f;
inline constexpr float daytimeLandTemperatureAnomalyC = 2.0f;
inline constexpr float nighttimeLandTemperatureAnomalyC = -2.0f;
inline constexpr float daytimeSeaTemperatureAnomalyC = 0.35f;
inline constexpr float nighttimeSeaTemperatureAnomalyC = -0.35f;
inline constexpr float persistentHeatingRetentionPerMillimetreC = 0.05f;
inline constexpr float persistentHeatingDecayTimeDays = 3.0f;
inline constexpr float maximumPersistentHeatingC = 2.5f;
inline constexpr float cloudMemoryTimeDays = 2.0f;
inline constexpr float maximumCloudRadiativeCoolingC = 1.25f;
inline constexpr float heatingVerticalVelocityResponse = 8.0f;
inline constexpr int moistAdjustmentIterations = 2;
inline constexpr float latentHeatingCPerMillimetre = 0.35f;
inline constexpr float saturationCapacityTemperatureSensitivityPerC = 0.065f;
inline constexpr float maximumFallingPrecipitationReevaporationFraction = 0.40f;
inline constexpr float allSnowTemperatureC = -1.0f;
inline constexpr float allRainTemperatureC = 2.0f;
inline constexpr float degreeDaySnowMeltMmPerDegreeC = 3.0f;
inline constexpr float maximumSnowStorageMm = 5000.0f;
inline constexpr float snowCoverConvergenceThresholdMm = 1.0f;
inline constexpr float fullSnowCoverWaterEquivalentMm = 10.0f;
inline constexpr float dynamicVerticalCooling = 0.025f;
inline constexpr float topographicUpliftCoolingCPerMetre = 0.0040f;
inline constexpr float dynamicSubsidenceWarming = 0.018f;
inline constexpr float topographicDescentWarmingCPerMetre = 0.0049f;
inline constexpr float maximumParcelTemperatureAdjustment = 8.0f;
inline constexpr float rainfallScale = 1.0f;
inline constexpr float seaIceFactor = 0.08f;
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
inline constexpr uint32_t erosionPeriod = 60;
inline constexpr float foldingRatio = 0.08f;
inline constexpr uint32_t aggregationOverlapAbsolute = 64;
inline constexpr float aggregationOverlapRelative = 0.20f;
inline constexpr uint32_t cycleCount = 2;
inline constexpr uint32_t plateCount = 10;
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
