#ifndef NOMINMAX
#define NOMINMAX
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include "heightmap_io.hpp"
#include "material.hpp"
#include "platecapi.hpp"
#include "sqrdmd.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <system_error>
#include <vector>

namespace fs = std::filesystem;

#ifndef UW_TECTONICS_OUTPUT_ROOT
#define UW_TECTONICS_OUTPUT_ROOT "."
#endif

namespace {

constexpr float kSeaLevel = 0.65f;
constexpr uint32_t kErosionPeriod = 60;
constexpr float kDefaultErosionStrength = 1.0f;
constexpr float kFoldingRatio = 0.08f;
constexpr uint32_t kMinAggregationOverlapAbs = 64;
constexpr uint32_t kAggregationOverlapAreaDivisor = 1000;
constexpr float kAggrOverlapRel = 0.20f;
constexpr uint32_t kDefaultCycleCount = 2;
constexpr uint32_t kDefaultCycleStepLimit = 600;
constexpr uint32_t kDefaultPlateCount = 10;
constexpr float kDefaultRotationStrength = 1.0f;
constexpr float kDefaultLandmassRotation = 0.20f;
constexpr float kDefaultSubductionStrength = 1.0f;
constexpr float kDefaultDivergentCarve = 0.015f;
constexpr uint16_t kAnimationDelayCs = 8;
constexpr float kDefaultMineralSurfaceTempC = 15.0f;
constexpr float kDefaultMineralGeothermalGradientCPerKm = 28.0f;
constexpr float kDefaultMineralPressureGradientKbarPerKm = 0.30f;
constexpr float kDefaultMineralHydrothermalActivity = 1.0f;
constexpr float kDefaultMineralFluidFlux = 1.0f;
constexpr float kDefaultMineralMetamorphismGain = 1.0f;
constexpr float kDefaultMineralNoise = 0.35f;

const fs::path kOutputRoot(UW_TECTONICS_OUTPUT_ROOT);
const fs::path kDefaultInputDir = kOutputRoot / "img" / "in";
const fs::path kDefaultOutputDir = kOutputRoot / "img" / "out";
const fs::path kDefaultAnimationDir = kOutputRoot / "img" / "anim";
const fs::path kDefaultFramesDir = kOutputRoot / "img" / "frames";

struct Params {
    uint32_t seed;
    uint32_t width;
    uint32_t height;
    bool colors;
    bool show_boundaries;
    uint32_t cycles;
    uint32_t cycle_steps;
    uint32_t plates;
    bool plates_explicit;
    uint32_t animated_step;
    uint32_t aggregation_overlap_abs;
    bool aggregation_overlap_abs_explicit;
    float aggregation_overlap_rel;
    uint32_t erosion_period;
    float folding_ratio;
    float erosion_strength;
    float rotation_strength;
    float landmass_rotation;
    float subduction_strength;
    float divergent_carve;
    double delta_time_myr;
    double cycle_duration_myr;
    float gravity_mps2;
    float glacial_erosion_strength;
    uint32_t glacial_erosion_period;
    float collision_uplift_ratio;
    float continental_compression_gain;
    float continental_boundary_fluidity;
    float movement_energy;
    int32_t crust_type_boundary_m;
    uint16_t max_initial_height_m;
    uint32_t hf_noise_period;
    float hf_noise_strength;
    uint32_t lf_noise_period;
    float lf_noise_strength;
    float mineral_surface_temp_c;
    float mineral_geothermal_gradient_c_per_km;
    float mineral_pressure_gradient_kbar_per_km;
    float mineral_hydrothermal_activity;
    float mineral_fluid_flux;
    float mineral_metamorphism_gain;
    float mineral_noise;
    uint32_t step;
    bool has_sea_level_m;
    uint16_t sea_level_m;
};

struct BoundaryOverlayData {
    std::vector<uint8_t> boundary_type;

    bool valid() const
    {
        return !boundary_type.empty();
    }
};

bool capture_boundary_overlay(void* simulation, int width, int height, BoundaryOverlayData& data);
void overlay_boundaries(const BoundaryOverlayData& data, int width, int height,
                        std::vector<png_byte>& rgb);

struct CheckpointRecord {
    std::string label;
    uint32_t step = 0;
    uint32_t iteration_count = 0;
    uint32_t cycle_count = 0;
    double time_myr = 0.0;
    std::string frame_file;
};

enum class MapArtifact : uint32_t {
    Height = 0U,
    Material = 1U,
    PlateType = 2U,
    PlateId = 3U,
    Boundary = 4U,
    CrustAge = 5U,
    Mineral = 6U,
    Tension = 7U,
};

constexpr size_t kMapArtifactCount = 8U;

const char* map_suffix(MapArtifact artifact)
{
    switch (artifact) {
    case MapArtifact::Height:
        return "HT";
    case MapArtifact::Material:
        return "MT";
    case MapArtifact::PlateType:
        return "PT";
    case MapArtifact::PlateId:
        return "ID";
    case MapArtifact::Boundary:
        return "B";
    case MapArtifact::CrustAge:
        return "CA";
    case MapArtifact::Mineral:
        return "MN";
    case MapArtifact::Tension:
        return "TN";
    default:
        return "B";
    }
}

[[noreturn]] void fail(const std::string& message)
{
    fprintf(stderr, "error: %s\n", message.c_str());
    exit(1);
}

std::string display_path(fs::path path)
{
    path.make_preferred();
    return path.string();
}

void ensure_console_output()
{
#ifdef _WIN32
    const HANDLE stdout_handle = GetStdHandle(STD_OUTPUT_HANDLE);
    const HANDLE stderr_handle = GetStdHandle(STD_ERROR_HANDLE);
    const bool stdout_ready =
        stdout_handle != nullptr && stdout_handle != INVALID_HANDLE_VALUE;
    const bool stderr_ready =
        stderr_handle != nullptr && stderr_handle != INVALID_HANDLE_VALUE;

    if (!stdout_ready || !stderr_ready) {
        if (!AttachConsole(ATTACH_PARENT_PROCESS)) {
            AllocConsole();
        }

        FILE* stream = nullptr;
        if (freopen_s(&stream, "CONOUT$", "w", stdout) == 0 && stream != nullptr) {
            setvbuf(stdout, nullptr, _IONBF, 0);
        }
        stream = nullptr;
        if (freopen_s(&stream, "CONOUT$", "w", stderr) == 0 && stream != nullptr) {
            setvbuf(stderr, nullptr, _IONBF, 0);
        }
    }
#endif
}

uint32_t parse_u32(const char* flag, const char* value)
{
    char* end = nullptr;
    const unsigned long parsed = strtoul(value, &end, 10);
    if (end == value || *end != '\0' || parsed == 0 ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint32_t>(parsed);
}

uint32_t parse_nonnegative_u32(const char* flag, const char* value)
{
    char* end = nullptr;
    const unsigned long parsed = strtoul(value, &end, 10);
    if (end == value || *end != '\0' ||
        parsed > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint32_t>(parsed);
}

uint16_t parse_u16(const char* flag, const char* value)
{
    const uint32_t parsed = parse_nonnegative_u32(flag, value);
    if (parsed > TopographyCodec::kMaxHeightMeters) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return static_cast<uint16_t>(parsed);
}

float parse_nonnegative_float(const char* flag, const char* value)
{
    char* end = nullptr;
    const float parsed = strtof(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed) || parsed < 0.0f) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

double parse_double(const char* flag, const char* value)
{
    char* end = nullptr;
    const double parsed = std::strtod(value, &end);
    if (end == value || *end != '\0' || !std::isfinite(parsed)) {
        fail(std::string("invalid ") + flag + ": " + value);
    }
    return parsed;
}

float parse_unit_float(const char* flag, const char* value)
{
    const float parsed = parse_nonnegative_float(flag, value);
    if (parsed > 1.0f) {
        fail(std::string("invalid value for ") + flag + ": " + value);
    }
    return parsed;
}

uint32_t default_aggregation_overlap_abs(uint32_t width, uint32_t height)
{
    const uint64_t area = static_cast<uint64_t>(width) * static_cast<uint64_t>(height);
    const uint32_t scaled = static_cast<uint32_t>(area / kAggregationOverlapAreaDivisor);
    return std::max(kMinAggregationOverlapAbs, scaled);
}

std::string json_escape(const std::string& value)
{
    std::string escaped;
    escaped.reserve(value.size());
    for (char ch : value) {
        switch (ch) {
        case '\\':
            escaped += "\\\\";
            break;
        case '\"':
            escaped += "\\\"";
            break;
        case '\n':
            escaped += "\\n";
            break;
        case '\r':
            escaped += "\\r";
            break;
        case '\t':
            escaped += "\\t";
            break;
        default:
            escaped += ch;
            break;
        }
    }
    return escaped;
}

float hash01(uint32_t seed, uint32_t x, uint32_t y, uint32_t salt)
{
    uint32_t value = seed ^ (x * 0x9E3779B9u) ^ (y * 0x85EBCA6Bu) ^ (salt * 0xC2B2AE35u);
    value ^= value >> 16;
    value *= 0x7FEB352Du;
    value ^= value >> 15;
    value *= 0x846CA68Bu;
    value ^= value >> 16;
    return static_cast<float>(value & 0x00FFFFFFu) / static_cast<float>(0x01000000u);
}

bool is_continental_material(uint8_t material_index)
{
    return platec::material::from_index(material_index) != platec::material::Type::Basalt;
}

void black_out_boundaries(const uint8_t* boundary_type_map, int width, int height,
                          std::vector<png_byte>& rgb)
{
    if (boundary_type_map == nullptr) {
        return;
    }

    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    for (size_t i = 0; i < pixel_count; ++i) {
        png_byte* ptr = &rgb[i * 3U];
        switch (static_cast<platec::contract::BoundaryType>(boundary_type_map[i])) {
        case platec::contract::BoundaryType::Convergent:
            ptr[0] = 255;
            ptr[1] = 64;
            ptr[2] = 64;
            break;
        case platec::contract::BoundaryType::Divergent:
            ptr[0] = 64;
            ptr[1] = 192;
            ptr[2] = 255;
            break;
        case platec::contract::BoundaryType::Transform:
            ptr[0] = 255;
            ptr[1] = 215;
            ptr[2] = 0;
            break;
        case platec::contract::BoundaryType::PassiveMargin:
            ptr[0] = 196;
            ptr[1] = 160;
            ptr[2] = 96;
            break;
        case platec::contract::BoundaryType::None:
        default:
            break;
        }
    }
}

void print_help()
{
    printf(" -h --help           : show this message\n");
    printf(" -s --seed SEED      : use the given SEED\n");
    printf(" -d --dim WIDTH HEIGHT: use the given width and height when no input image is used\n");
    printf(" -sl --sea-level-m N : metric coastline threshold in [0, 65535]\n");
    printf(" -co --color         : colorize preview frames and animated PNG output\n");
    printf(" -c --cycles N       : number of simulation cycles to run (default %u)\n",
           kDefaultCycleCount);
    printf(" --cycle-steps N     : max updates per cycle before restart; 0 disables (default %u)\n",
           kDefaultCycleStepLimit);
    printf(" -p --plates N       : total number of tectonic plates (default %u)\n",
           kDefaultPlateCount);
    printf(" -aa --aggregation-overlap-abs N: overlap pixels needed to aggregate continents\n");
    printf("                       default max(%u, area/%u); %u at 600x400\n",
           kMinAggregationOverlapAbs, kAggregationOverlapAreaDivisor,
           default_aggregation_overlap_abs(600, 400));
    printf(" -ar --aggregation-overlap-rel X: overlap ratio needed to aggregate continents in [0, 1]\n");
    printf("                       default %.2f\n", kAggrOverlapRel);
    printf(" -fr --folding-ratio X: fraction of overlapping continental crust turned into uplift in [0, 1]\n");
    printf("                       default %.2f\n", kFoldingRatio);
    printf(" -wep --wind-erosion-period N: number of simulation updates between wind erosion passes (default %u)\n",
           kErosionPeriod);
    printf(" -wes --wind-erosion-strength X: scale wind erosion amount (0 disables, default %.2f)\n",
           kDefaultErosionStrength);
    printf(" -rs --rotation-strength X: scale angular plate motion (default %.2f)\n",
           kDefaultRotationStrength);
    printf(" -lr --landmass-rotation X: scale visible crust rotation (0 disables, default %.2f)\n",
           kDefaultLandmassRotation);
    printf(" -ss --subduction-strength X: scale oceanic crust removal during subduction in [0, 1]\n");
    printf("                       default %.2f\n", kDefaultSubductionStrength);
    printf(" -dc --divergent-carve X: downward carve per divergent regeneration step (default %.3f)\n",
           kDefaultDivergentCarve);
    printf(" -dt --delta-time-myr X: geological time per simulation update (default %.2f)\n",
           platec::scenario::kDefaultDeltaTimeMyr);
    printf(" -cd --cycle-duration-myr X: cap one cycle to X Myr by deriving the step limit (0 disables)\n");
    printf(" -gr --gravity X     : planetary gravity in m/s^2 (default %.5f)\n",
           platec::scenario::kDefaultGravityMps2);
    printf(" -ges --glacial-erosion-strength X: glacial buzzsaw erosion strength (default %.2f)\n",
           platec::scenario::kDefaultGlacialErosionStrength);
    printf(" -gep --glacial-erosion-period N: number of simulation updates between glacial erosion passes (default %u)\n",
           platec::scenario::kDefaultGlacialErosionPeriod);
    printf(" -cq --collision-uplift-ratio X: convergence uplift ratio q in [0, 1] (default %.2f)\n",
           kFoldingRatio);
    printf(" -cg --cc-compression-gain X: exponential continental compression uplift gain (default %.2f)\n",
           platec::scenario::kDefaultContinentalCompressionGain);
    printf(" -cf --cc-boundary-fluidity X: continental boundary fluidity in [0, 1] (default %.2f)\n",
           platec::scenario::kDefaultContinentalBoundaryFluidity);
    printf(" -me --movement-energy X: scales retained plate motion after deformation (default %.2f)\n",
           platec::scenario::kDefaultMovementEnergy);
    printf(" -ct --crust-type-boundary N: initial continental/oceanic threshold in metres [0, 65535]\n");
    printf(" -max --max-initial-height N: initial terrain amplitude mapped from uint16 domain [0, 65535]\n");
    printf("                       default %u\n", TopographyCodec::kMaxHeightMeters);
    printf(" -hp --hf-noise-period N: apply high-frequency deterministic noise every N updates (0 disables)\n");
    printf(" -hs --hf-noise-strength X: high-frequency deterministic noise strength (default %.2f)\n",
           platec::scenario::kDefaultPeriodicNoiseStrength);
    printf(" -lp --lf-noise-period N: apply low-frequency deterministic noise every N updates (0 disables)\n");
    printf(" -ls --lf-noise-strength X: low-frequency deterministic noise strength (default %.2f)\n",
           platec::scenario::kDefaultPeriodicNoiseStrength);
    printf(" -mt --mineral-surface-temp-c X: mineral surface temperature baseline (default %.2f)\n",
           kDefaultMineralSurfaceTempC);
    printf(" -mg --mineral-geothermal-gradient X: mineral geothermal gradient C/km (default %.2f)\n",
           kDefaultMineralGeothermalGradientCPerKm);
    printf(" -mp --mineral-pressure-gradient X: mineral pressure gradient kbar/km (default %.2f)\n",
           kDefaultMineralPressureGradientKbarPerKm);
    printf(" -mh --mineral-hydrothermal X: mineral hydrothermal activity scale (default %.2f)\n",
           kDefaultMineralHydrothermalActivity);
    printf(" -mf --mineral-fluid-flux X: mineral subduction fluid flux scale (default %.2f)\n",
           kDefaultMineralFluidFlux);
    printf(" -mm --mineral-metamorphism X: mineral metamorphism gain (default %.2f)\n",
           kDefaultMineralMetamorphismGain);
    printf(" -mn --mineral-noise X: mineral deterministic noise strength (default %.2f)\n",
           kDefaultMineralNoise);
    printf(" -st --step X        : save intermediate maps every X steps\n");
    printf(" -a --animated X     : export animated PNGs sampling every X steps\n");
    printf(" -b --boundaries     : export transparent boundary PNGs and overlay preview/animation frames\n");
    printf(" output files        : img/in/<YYMMDDHHMM>_<seed>_HT.tiff\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_HT.tiff\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_MT.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_MT.png\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_PT.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_PT.png\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_ID.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_ID.png\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_B.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_B.png\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_CA.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_CA.png\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_MN.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_MN.png\n");
    printf("                       img/in/<YYMMDDHHMM>_<seed>_TN.png\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>_TN.png\n");
    printf("                       TIFF is grayscale uint16 metric elevation in metres\n");
    printf("                       img/out/<YYMMDDHHMM>_<seed>.run.json\n");
    printf("                       img/anim/<YYMMDDHHMM>_<seed>_<MAP>.png when --animated is used\n");
    printf("                       img/frames/<YYMMDDHHMM>_<seed>_<frame>_<MAP>.png for preview frames\n");
}

Params fill_params(int argc, char* argv[])
{
    srand(static_cast<unsigned int>(time(nullptr)));

    Params params = {
        static_cast<uint32_t>(rand()),
        600,
        400,
        false,
        false,
        kDefaultCycleCount,
        kDefaultCycleStepLimit,
        kDefaultPlateCount,
        false,
        0,
        default_aggregation_overlap_abs(600, 400),
        false,
        kAggrOverlapRel,
        kErosionPeriod,
        kFoldingRatio,
        kDefaultErosionStrength,
        kDefaultRotationStrength,
        kDefaultLandmassRotation,
        kDefaultSubductionStrength,
        kDefaultDivergentCarve,
        platec::scenario::kDefaultDeltaTimeMyr,
        0.0,
        platec::scenario::kDefaultGravityMps2,
        platec::scenario::kDefaultGlacialErosionStrength,
        platec::scenario::kDefaultGlacialErosionPeriod,
        kFoldingRatio,
        platec::scenario::kDefaultContinentalCompressionGain,
        platec::scenario::kDefaultContinentalBoundaryFluidity,
        platec::scenario::kDefaultMovementEnergy,
        platec::scenario::kAutoCrustTypeBoundaryMeters,
        TopographyCodec::kMaxHeightMeters,
        platec::scenario::kDefaultPeriodicNoisePeriod,
        platec::scenario::kDefaultPeriodicNoiseStrength,
        platec::scenario::kDefaultPeriodicNoisePeriod,
        platec::scenario::kDefaultPeriodicNoiseStrength,
        kDefaultMineralSurfaceTempC,
        kDefaultMineralGeothermalGradientCPerKm,
        kDefaultMineralPressureGradientKbarPerKm,
        kDefaultMineralHydrothermalActivity,
        kDefaultMineralFluidFlux,
        kDefaultMineralMetamorphismGain,
        kDefaultMineralNoise,
        0,
        false,
        0
    };

    int p = 1;
    while (p < argc) {
        if (strcmp(argv[p], "--help") == 0 || strcmp(argv[p], "-h") == 0) {
            print_help();
            exit(0);
        } else if (strcmp(argv[p], "-s") == 0 || strcmp(argv[p], "--seed") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --seed");
            }
            params.seed = parse_u32("--seed", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--dim") == 0 || strcmp(argv[p], "-d") == 0) {
            if (p + 2 >= argc) {
                fail("two parameters should follow --dim");
            }
            params.width = parse_u32("--dim", argv[p + 1]);
            params.height = parse_u32("--dim", argv[p + 2]);
            if (params.width < 5 || params.height < 5) {
                fail("dimensions have to be >= 5");
            }
            p += 3;
        } else if (strcmp(argv[p], "--sea-level-m") == 0 || strcmp(argv[p], "-sl") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --sea-level-m");
            }
            params.has_sea_level_m = true;
            params.sea_level_m = parse_u16("--sea-level-m", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--color") == 0 || strcmp(argv[p], "--colors") == 0 ||
                   strcmp(argv[p], "-co") == 0) {
            params.colors = true;
            p += 1;
        } else if (strcmp(argv[p], "--cycles") == 0 || strcmp(argv[p], "-c") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --cycles");
            }
            params.cycles = parse_u32("--cycles", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--cycle-steps") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --cycle-steps");
            }
            params.cycle_steps = parse_nonnegative_u32("--cycle-steps", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--plates") == 0 || strcmp(argv[p], "-p") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --plates");
            }
            params.plates = parse_u32("--plates", argv[p + 1]);
            params.plates_explicit = true;
            p += 2;
        } else if (strcmp(argv[p], "--aggregation-overlap-abs") == 0 ||
                   strcmp(argv[p], "-aa") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --aggregation-overlap-abs");
            }
            params.aggregation_overlap_abs =
                parse_nonnegative_u32("--aggregation-overlap-abs", argv[p + 1]);
            params.aggregation_overlap_abs_explicit = true;
            p += 2;
        } else if (strcmp(argv[p], "--aggregation-overlap-rel") == 0 ||
                   strcmp(argv[p], "-ar") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --aggregation-overlap-rel");
            }
            params.aggregation_overlap_rel =
                parse_unit_float("--aggregation-overlap-rel", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--folding-ratio") == 0 || strcmp(argv[p], "-fr") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --folding-ratio");
            }
            params.folding_ratio = parse_unit_float("--folding-ratio", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--wind-erosion-period") == 0 ||
                   strcmp(argv[p], "-wep") == 0 || strcmp(argv[p], "--erosion-period") == 0 ||
                   strcmp(argv[p], "-ep") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --wind-erosion-period");
            }
            params.erosion_period = parse_u32("--wind-erosion-period", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--wind-erosion-strength") == 0 ||
                   strcmp(argv[p], "-wes") == 0 || strcmp(argv[p], "--erosion-strength") == 0 ||
                   strcmp(argv[p], "-es") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --wind-erosion-strength");
            }
            params.erosion_strength =
                parse_nonnegative_float("--wind-erosion-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--rotation-strength") == 0 || strcmp(argv[p], "-rs") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --rotation-strength");
            }
            params.rotation_strength = parse_nonnegative_float("--rotation-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--landmass-rotation") == 0 ||
                   strcmp(argv[p], "-lr") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --landmass-rotation");
            }
            params.landmass_rotation =
                parse_nonnegative_float("--landmass-rotation", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--subduction-strength") == 0 ||
                   strcmp(argv[p], "-ss") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --subduction-strength");
            }
            params.subduction_strength =
                parse_unit_float("--subduction-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--divergent-carve") == 0 || strcmp(argv[p], "-dc") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --divergent-carve");
            }
            params.divergent_carve =
                parse_nonnegative_float("--divergent-carve", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--delta-time-myr") == 0 || strcmp(argv[p], "-dt") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --delta-time-myr");
            }
            params.delta_time_myr = parse_double("--delta-time-myr", argv[p + 1]);
            if (params.delta_time_myr <= 0.0) {
                fail("--delta-time-myr must be > 0");
            }
            p += 2;
        } else if (strcmp(argv[p], "--cycle-duration-myr") == 0 || strcmp(argv[p], "-cd") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --cycle-duration-myr");
            }
            params.cycle_duration_myr = parse_double("--cycle-duration-myr", argv[p + 1]);
            if (params.cycle_duration_myr < 0.0) {
                fail("--cycle-duration-myr must be >= 0");
            }
            p += 2;
        } else if (strcmp(argv[p], "--gravity") == 0 || strcmp(argv[p], "-gr") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --gravity");
            }
            params.gravity_mps2 = parse_nonnegative_float("--gravity", argv[p + 1]);
            if (params.gravity_mps2 <= 0.0f) {
                fail("--gravity must be > 0");
            }
            p += 2;
        } else if (strcmp(argv[p], "--glacial-erosion-strength") == 0 ||
                   strcmp(argv[p], "-ges") == 0 || strcmp(argv[p], "--glacial-erosion") == 0 ||
                   strcmp(argv[p], "-ge") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --glacial-erosion-strength");
            }
            params.glacial_erosion_strength =
                parse_nonnegative_float("--glacial-erosion-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--glacial-erosion-period") == 0 ||
                   strcmp(argv[p], "-gep") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --glacial-erosion-period");
            }
            params.glacial_erosion_period =
                parse_u32("--glacial-erosion-period", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--collision-uplift-ratio") == 0 ||
                   strcmp(argv[p], "-cq") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --collision-uplift-ratio");
            }
            params.collision_uplift_ratio =
                parse_unit_float("--collision-uplift-ratio", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--cc-compression-gain") == 0 ||
                   strcmp(argv[p], "-cg") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --cc-compression-gain");
            }
            params.continental_compression_gain =
                parse_nonnegative_float("--cc-compression-gain", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--cc-boundary-fluidity") == 0 ||
                   strcmp(argv[p], "-cf") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --cc-boundary-fluidity");
            }
            params.continental_boundary_fluidity =
                parse_unit_float("--cc-boundary-fluidity", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--movement-energy") == 0 ||
                   strcmp(argv[p], "-me") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --movement-energy");
            }
            params.movement_energy =
                parse_nonnegative_float("--movement-energy", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--crust-type-boundary") == 0 ||
                   strcmp(argv[p], "-ct") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --crust-type-boundary");
            }
            params.crust_type_boundary_m =
                static_cast<int32_t>(parse_u16("--crust-type-boundary", argv[p + 1]));
            p += 2;
        } else if (strcmp(argv[p], "--max-initial-height") == 0 ||
                   strcmp(argv[p], "-max") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --max-initial-height");
            }
            params.max_initial_height_m =
                parse_u16("--max-initial-height", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--hf-noise-period") == 0 || strcmp(argv[p], "-hp") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --hf-noise-period");
            }
            params.hf_noise_period =
                parse_nonnegative_u32("--hf-noise-period", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--hf-noise-strength") == 0 || strcmp(argv[p], "-hs") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --hf-noise-strength");
            }
            params.hf_noise_strength =
                parse_nonnegative_float("--hf-noise-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--lf-noise-period") == 0 || strcmp(argv[p], "-lp") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --lf-noise-period");
            }
            params.lf_noise_period =
                parse_nonnegative_u32("--lf-noise-period", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--lf-noise-strength") == 0 || strcmp(argv[p], "-ls") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --lf-noise-strength");
            }
            params.lf_noise_strength =
                parse_nonnegative_float("--lf-noise-strength", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-surface-temp-c") == 0 ||
                   strcmp(argv[p], "-mt") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-surface-temp-c");
            }
            params.mineral_surface_temp_c =
                parse_nonnegative_float("--mineral-surface-temp-c", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-geothermal-gradient") == 0 ||
                   strcmp(argv[p], "-mg") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-geothermal-gradient");
            }
            params.mineral_geothermal_gradient_c_per_km =
                parse_nonnegative_float("--mineral-geothermal-gradient", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-pressure-gradient") == 0 ||
                   strcmp(argv[p], "-mp") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-pressure-gradient");
            }
            params.mineral_pressure_gradient_kbar_per_km =
                parse_nonnegative_float("--mineral-pressure-gradient", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-hydrothermal") == 0 ||
                   strcmp(argv[p], "-mh") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-hydrothermal");
            }
            params.mineral_hydrothermal_activity =
                parse_nonnegative_float("--mineral-hydrothermal", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-fluid-flux") == 0 ||
                   strcmp(argv[p], "-mf") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-fluid-flux");
            }
            params.mineral_fluid_flux =
                parse_nonnegative_float("--mineral-fluid-flux", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-metamorphism") == 0 ||
                   strcmp(argv[p], "-mm") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-metamorphism");
            }
            params.mineral_metamorphism_gain =
                parse_nonnegative_float("--mineral-metamorphism", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--mineral-noise") == 0 || strcmp(argv[p], "-mn") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --mineral-noise");
            }
            params.mineral_noise =
                parse_nonnegative_float("--mineral-noise", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--step") == 0 || strcmp(argv[p], "-st") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --step");
            }
            params.step = parse_u32("--step", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--animated") == 0 || strcmp(argv[p], "-a") == 0) {
            if (p + 1 >= argc) {
                fail("a parameter should follow --animated");
            }
            params.animated_step = parse_u32("--animated", argv[p + 1]);
            p += 2;
        } else if (strcmp(argv[p], "--boundaries") == 0 || strcmp(argv[p], "-b") == 0) {
            params.show_boundaries = true;
            p += 1;
        } else {
            fail(std::string("unexpected param '") + argv[p] + "', use -h to display a list of params");
        }
    }

    if (params.cycle_duration_myr > 0.0) {
        const double derived_steps = std::ceil(params.cycle_duration_myr / params.delta_time_myr);
        if (!std::isfinite(derived_steps) ||
            derived_steps > static_cast<double>(std::numeric_limits<uint32_t>::max())) {
            fail("derived cycle step limit is out of range");
        }
        params.cycle_steps = std::max<uint32_t>(1U, static_cast<uint32_t>(derived_steps));
    }

    return params;
}

void ensure_output_directories()
{
    fs::create_directories(kDefaultInputDir);
    fs::create_directories(kDefaultOutputDir);
    fs::create_directories(kDefaultAnimationDir);
    fs::create_directories(kDefaultFramesDir);
}

std::string sanitize_file_component(const std::string& value)
{
    std::string sanitized;
    sanitized.reserve(value.size());

    for (char ch : value) {
        const unsigned char byte = static_cast<unsigned char>(ch);
        if (std::isalnum(byte) != 0 || ch == '-' || ch == '_') {
            sanitized += ch;
        } else if (sanitized.empty() || sanitized.back() != '_') {
            sanitized += '_';
        }
    }

    while (!sanitized.empty() && sanitized.back() == '_') {
        sanitized.pop_back();
    }

    return sanitized.empty() ? "run" : sanitized;
}

std::string make_run_timestamp(std::time_t now)
{
    std::tm local_tm = {};
#ifdef _WIN32
    localtime_s(&local_tm, &now);
#else
    localtime_r(&now, &local_tm);
#endif

    char buffer[32];
    if (std::strftime(buffer, sizeof(buffer), "%y%m%d%H%M%S", &local_tm) == 0) {
        fail("failed to format run timestamp");
    }
    return buffer;
}

std::string make_run_id(const Params& params, std::time_t started_at)
{
    return make_run_timestamp(started_at) + "_" +
           sanitize_file_component(std::to_string(params.seed));
}

std::string frame_file_name(const std::string& run_id, uint32_t frame_index,
                            MapArtifact artifact)
{
    char buffer[128];
    snprintf(buffer, sizeof(buffer), "%s_%05u_%s.png", run_id.c_str(), frame_index,
             map_suffix(artifact));
    return buffer;
}

std::string animation_file_name(const std::string& run_id, MapArtifact artifact)
{
    return run_id + "_" + map_suffix(artifact) + ".png";
}

std::vector<uint16_t> to_metric_heightmap(const float* heightmap, int width, int height,
                                          uint16_t sea_level_m)
{
    const size_t sample_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<uint16_t> metric_heightmap(sample_count);
    for (size_t i = 0; i < sample_count; ++i) {
        metric_heightmap[i] = TopographyCodec::internal_to_meters(heightmap[i], sea_level_m);
    }
    return metric_heightmap;
}

void export_metric_heightmap(const float* heightmap, int width, int height, uint16_t sea_level_m,
                             const fs::path& image_path)
{
    const std::vector<uint16_t> metric_heightmap =
        to_metric_heightmap(heightmap, width, height, sea_level_m);

    if (writeImageGrayTiff16(image_path.string().c_str(), width, height,
                             metric_heightmap.data()) != 0) {
        fail("failed to write metric TIFF16: " + display_path(image_path));
    }

}

void export_material_map(const uint8_t* material_map, int width, int height,
                         const fs::path& image_path)
{
    if (material_map == nullptr) {
        fail("material map is not available");
    }

    std::vector<png_byte> rgb;
    renderMaterialMapRgb(width, height, material_map, rgb);
    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Material Map") != 0) {
        fail("failed to write material map: " + display_path(image_path));
    }
}

void export_plate_id_map(const uint32_t* plate_map, int width, int height,
                         const fs::path& image_path)
{
    if (plate_map == nullptr) {
        fail("plate map is not available");
    }

    std::vector<png_byte> rgb;
    renderPlateIdMapRgb(width, height, plate_map, rgb);
    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Plate IDs") != 0) {
        fail("failed to write plate ID map: " + display_path(image_path));
    }
}

void export_crust_type_map(void* simulation, int width, int height, const fs::path& image_path,
                           bool overlay_boundaries)
{
    const uint8_t* surface_crust_type_map = platec_api_get_surface_crust_type_map(simulation);
    if (surface_crust_type_map == nullptr) {
        fail("surface crust type map is not available");
    }

    std::vector<png_byte> rgb;
    renderCrustClassMapRgb(width, height, surface_crust_type_map, rgb);
    if (overlay_boundaries) {
        black_out_boundaries(platec_api_get_boundary_type_map(simulation), width, height, rgb);
    }
    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Crust Types") != 0) {
        fail("failed to write crust type map: " + display_path(image_path));
    }
}

void export_crust_age_map(const float* crust_age_myr_map, int width, int height,
                          const fs::path& image_path)
{
    if (crust_age_myr_map == nullptr) {
        fail("crust age map is not available");
    }

    std::vector<png_byte> rgb;
    renderCrustAgeMapRgb(width, height, crust_age_myr_map, rgb);
    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Crust Age") != 0) {
        fail("failed to write crust age map: " + display_path(image_path));
    }
}

void export_tension_map(const float* tension_map, int width, int height,
                        const fs::path& image_path)
{
    if (tension_map == nullptr) {
        fail("tension map is not available");
    }

    const size_t pixel_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    std::vector<float> normalized(pixel_count, 0.0f);
    float max_tension = 0.0f;
    for (size_t index = 0; index < pixel_count; ++index) {
        max_tension = std::max(max_tension, tension_map[index]);
    }

    const float denom = std::log1pf(std::max(1.0f, max_tension));
    for (size_t index = 0; index < pixel_count; ++index) {
        normalized[index] =
            denom > 0.0f ? std::log1pf(std::max(0.0f, tension_map[index])) / denom : 0.0f;
    }

    std::vector<png_byte> rgb;
    renderImageGrayRgb(width, height, normalized.data(), rgb);
    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Tension") != 0) {
        fail("failed to write tension map: " + display_path(image_path));
    }
}

void export_mineral_map(void* simulation, const Params& params, int width, int height,
                        const fs::path& image_path)
{
    const float* heightmap = platec_api_get_heightmap(simulation);
    const uint8_t* material_map = platec_api_get_material_map(simulation);
    const uint8_t* crust_class_map = platec_api_get_crust_class_map(simulation);
    const float* crust_age_myr_map = platec_api_get_crust_age_myr_map(simulation);
    const float* crust_thickness_map = platec_api_get_crust_thickness_map(simulation);
    const float* accumulated_strain_map = platec_api_get_accumulated_strain_map(simulation);
    const uint8_t* boundary_type_map = platec_api_get_boundary_type_map(simulation);
    const uint8_t* geologic_regime_map = platec_api_get_geologic_regime_map(simulation);
    if (heightmap == nullptr || material_map == nullptr || crust_class_map == nullptr ||
        crust_age_myr_map == nullptr || crust_thickness_map == nullptr ||
        accumulated_strain_map == nullptr || boundary_type_map == nullptr ||
        geologic_regime_map == nullptr) {
        fail("mineral map inputs are not available");
    }

    struct MineralColor {
        png_byte r;
        png_byte g;
        png_byte b;
    };
    static const MineralColor palette[] = {
        {0, 0, 0},
        {0xD9, 0xF0, 0xFF},
        {0xF2, 0xB7, 0xA1},
        {0x8B, 0x6F, 0x9C},
        {0x3F, 0x6F, 0x67},
        {0x6A, 0x8B, 0x52},
        {0xA7, 0xC9, 0x57},
        {0xFF, 0xF2, 0xCC},
        {0x8F, 0x1D, 0x21},
    };

    std::vector<png_byte> rgb(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U, 0U);
    const uint16_t sea_level_m = platec_api_get_sea_level_m(simulation);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t index =
                static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x);
            png_byte* ptr = &rgb[index * 3U];
            if (crust_class_map[index] == static_cast<uint8_t>(platec::contract::CrustClass::None)) {
                continue;
            }

            const float elevation_m =
                static_cast<float>(TopographyCodec::internal_to_meters(heightmap[index], sea_level_m)) -
                static_cast<float>(sea_level_m);
            const bool continental = is_continental_material(material_map[index]);
            const float age_myr = std::max(0.0f, crust_age_myr_map[index]);
            const float strain = std::max(0.0f, accumulated_strain_map[index]);
            const float thickness = std::max(0.0f, crust_thickness_map[index]);
            const auto boundary_type =
                static_cast<platec::contract::BoundaryType>(boundary_type_map[index]);
            const auto regime =
                static_cast<platec::contract::GeologicRegime>(geologic_regime_map[index]);
            const float noise =
                (hash01(params.seed, static_cast<uint32_t>(x), static_cast<uint32_t>(y), 17U) * 2.0f -
                 1.0f) *
                params.mineral_noise;

            float burial_km =
                std::max(0.0f, thickness * 6.0f + strain * 10.0f + age_myr * 0.020f -
                                   elevation_m / 2000.0f);
            if (boundary_type == platec::contract::BoundaryType::Convergent) {
                burial_km += 9.0f * params.mineral_metamorphism_gain;
            } else if (boundary_type == platec::contract::BoundaryType::Transform) {
                burial_km += 2.0f * params.mineral_metamorphism_gain;
            } else if (boundary_type == platec::contract::BoundaryType::Divergent) {
                burial_km += 1.5f * params.mineral_hydrothermal_activity;
            }

            float temperature_c =
                params.mineral_surface_temp_c +
                params.mineral_geothermal_gradient_c_per_km * burial_km;
            if (boundary_type == platec::contract::BoundaryType::Divergent ||
                regime == platec::contract::GeologicRegime::MidOceanRidge) {
                temperature_c += 220.0f * params.mineral_hydrothermal_activity;
            }
            if (boundary_type == platec::contract::BoundaryType::Convergent) {
                temperature_c += 90.0f * params.mineral_metamorphism_gain;
            }

            float pressure_kbar =
                params.mineral_pressure_gradient_kbar_per_km * burial_km;
            if (regime == platec::contract::GeologicRegime::TrenchAdjacent) {
                pressure_kbar += 3.0f * params.mineral_fluid_flux;
            }

            temperature_c *= 1.0f + 0.15f * noise;
            pressure_kbar *= 1.0f + 0.10f * noise;

            size_t palette_index = 1U;
            if (!continental && temperature_c >= 950.0f) {
                palette_index = 6U;
            } else if ((!continental && temperature_c >= 700.0f) || elevation_m < -2000.0f) {
                palette_index = 5U;
            } else if ((boundary_type == platec::contract::BoundaryType::Convergent ||
                        regime == platec::contract::GeologicRegime::TrenchAdjacent) &&
                       temperature_c >= 500.0f) {
                palette_index = pressure_kbar >= 5.0f ? 8U : 4U;
            } else if (platec::material::from_index(material_map[index]) ==
                           platec::material::Type::Sedimentary &&
                       temperature_c <= 320.0f && pressure_kbar <= 4.0f) {
                palette_index = 7U;
            } else if (temperature_c >= 450.0f && continental) {
                palette_index = 2U;
            } else if ((pressure_kbar >= 1.0f || strain >= 0.3f) && temperature_c >= 250.0f) {
                palette_index = 3U;
            }

            ptr[0] = palette[palette_index].r;
            ptr[1] = palette[palette_index].g;
            ptr[2] = palette[palette_index].b;
        }
    }

    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Mineral Map") != 0) {
        fail("failed to write mineral map: " + display_path(image_path));
    }
}

void render_boundary_map_rgb(const BoundaryOverlayData& data, int width, int height,
                             std::vector<png_byte>& rgb)
{
    rgb.assign(static_cast<size_t>(width) * static_cast<size_t>(height) * 3U, 0U);
    if (data.valid()) {
        overlay_boundaries(data, width, height, rgb);
    }
}

void export_boundary_map(void* simulation, int width, int height, const fs::path& image_path,
                         const BoundaryOverlayData* fallback_boundaries = nullptr)
{
    BoundaryOverlayData boundaries;
    if (!capture_boundary_overlay(simulation, width, height, boundaries) &&
        fallback_boundaries != nullptr && fallback_boundaries->valid()) {
        boundaries = *fallback_boundaries;
    }

    std::vector<png_byte> rgb;
    render_boundary_map_rgb(boundaries, width, height, rgb);
    if (writeImageRgb(image_path.string().c_str(), width, height, rgb.data(), "Plate Boundaries") != 0) {
        fail("failed to write boundary map: " + display_path(image_path));
    }
}

void write_run_manifest(const fs::path& path, const Params& params, const std::string& run_id,
                        uint16_t sea_level_m,
                        uint32_t time_origin_step, uint32_t iteration_count,
                        uint32_t cycle_count, double time_myr,
                        const std::vector<CheckpointRecord>& checkpoints,
                        const fs::path& initial_image_output_path,
                        const fs::path& initial_material_output_path,
                        const fs::path& initial_plate_output_path,
                        const fs::path& initial_plate_class_output_path,
                        const fs::path* initial_boundary_output_path,
                        const fs::path& final_image_output_path,
                        const fs::path& final_material_output_path,
                        const fs::path& final_plate_output_path,
                        const fs::path& final_plate_class_output_path,
                        const fs::path* final_boundary_output_path,
                        const fs::path* animation_output_path)
{
    std::ofstream output(path);
    if (!output) {
        fail("failed to write run manifest: " + display_path(path));
    }

    output << std::setprecision(17);
    output << "{\n";
    output << "  \"schema\": \"plate-tectonics-run/v1\",\n";
    output << "  \"run_id\": \"" << json_escape(run_id) << "\",\n";
    output << "  \"seed\": " << params.seed << ",\n";
    output << "  \"width\": " << params.width << ",\n";
    output << "  \"height\": " << params.height << ",\n";
    output << "  \"sea_level_m\": " << sea_level_m << ",\n";
    output << "  \"time_origin_step\": " << time_origin_step << ",\n";
    output << "  \"iteration_count\": " << iteration_count << ",\n";
    output << "  \"cycle_count\": " << cycle_count << ",\n";
    output << "  \"time_myr\": " << time_myr << ",\n";
    output << "  \"delta_time_myr\": " << params.delta_time_myr << ",\n";
    output << "  \"input\": \"\",\n";
    output << "  \"scenario\": {\n";
    output << "    \"cycles\": " << params.cycles << ",\n";
    output << "    \"cycle_steps\": " << params.cycle_steps << ",\n";
    output << "    \"plates\": " << params.plates << ",\n";
    output << "    \"aggregation_overlap_abs\": " << params.aggregation_overlap_abs << ",\n";
    output << "    \"aggregation_overlap_rel\": " << params.aggregation_overlap_rel << ",\n";
    output << "    \"wind_erosion_period\": " << params.erosion_period << ",\n";
    output << "    \"folding_ratio\": " << params.folding_ratio << ",\n";
    output << "    \"wind_erosion_strength\": " << params.erosion_strength << ",\n";
    output << "    \"rotation_strength\": " << params.rotation_strength << ",\n";
    output << "    \"landmass_rotation\": " << params.landmass_rotation << ",\n";
    output << "    \"subduction_strength\": " << params.subduction_strength << ",\n";
    output << "    \"divergent_carve\": " << params.divergent_carve << ",\n";
    output << "    \"delta_time_myr\": " << params.delta_time_myr << ",\n";
    output << "    \"cycle_duration_myr\": " << params.cycle_duration_myr << ",\n";
    output << "    \"gravity_mps2\": " << params.gravity_mps2 << ",\n";
    output << "    \"glacial_erosion_strength\": " << params.glacial_erosion_strength << ",\n";
    output << "    \"glacial_erosion_period\": " << params.glacial_erosion_period << ",\n";
    output << "    \"collision_uplift_ratio\": " << params.collision_uplift_ratio << ",\n";
    output << "    \"continental_compression_gain\": " << params.continental_compression_gain
           << ",\n";
    output << "    \"continental_boundary_fluidity\": "
           << params.continental_boundary_fluidity << ",\n";
    output << "    \"movement_energy\": " << params.movement_energy << ",\n";
    output << "    \"crust_type_boundary_m\": " << params.crust_type_boundary_m << ",\n";
    output << "    \"max_initial_height_m\": " << params.max_initial_height_m << ",\n";
    output << "    \"hf_noise_period\": " << params.hf_noise_period << ",\n";
    output << "    \"hf_noise_strength\": " << params.hf_noise_strength << ",\n";
    output << "    \"lf_noise_period\": " << params.lf_noise_period << ",\n";
    output << "    \"lf_noise_strength\": " << params.lf_noise_strength << ",\n";
    output << "    \"mineral_surface_temp_c\": " << params.mineral_surface_temp_c << ",\n";
    output << "    \"mineral_geothermal_gradient_c_per_km\": "
           << params.mineral_geothermal_gradient_c_per_km << ",\n";
    output << "    \"mineral_pressure_gradient_kbar_per_km\": "
           << params.mineral_pressure_gradient_kbar_per_km << ",\n";
    output << "    \"mineral_hydrothermal_activity\": "
           << params.mineral_hydrothermal_activity << ",\n";
    output << "    \"mineral_fluid_flux\": " << params.mineral_fluid_flux << ",\n";
    output << "    \"mineral_metamorphism_gain\": " << params.mineral_metamorphism_gain
           << ",\n";
    output << "    \"mineral_noise\": " << params.mineral_noise << ",\n";
    output << "    \"use_material_map_for_height_limit\": false\n";
    output << "  },\n";
    output << "  \"outputs\": {\n";
    output << "    \"initial_heightmap_tiff\": \""
           << json_escape(initial_image_output_path.filename().string()) << "\",\n";
    output << "    \"initial_material_png\": \""
           << json_escape(initial_material_output_path.filename().string()) << "\",\n";
    output << "    \"initial_id_png\": \""
           << json_escape(initial_plate_output_path.filename().string()) << "\",\n";
    output << "    \"initial_crust_type_png\": \""
           << json_escape(initial_plate_class_output_path.filename().string()) << "\",\n";
    output << "    \"initial_crust_age_png\": \"" << json_escape(run_id + "_CA.png") << "\",\n";
    output << "    \"initial_mineral_png\": \"" << json_escape(run_id + "_MN.png") << "\",\n";
    output << "    \"initial_tension_png\": \"" << json_escape(run_id + "_TN.png") << "\",\n";
    if (initial_boundary_output_path != nullptr) {
        output << "    \"initial_boundary_png\": \""
               << json_escape(initial_boundary_output_path->filename().string()) << "\",\n";
    }
    output << "    \"final_heightmap_tiff\": \""
           << json_escape(final_image_output_path.filename().string()) << "\",\n";
    output << "    \"final_material_png\": \""
           << json_escape(final_material_output_path.filename().string()) << "\",\n";
    output << "    \"final_id_png\": \""
           << json_escape(final_plate_output_path.filename().string()) << "\",\n";
    output << "    \"final_crust_type_png\": \""
           << json_escape(final_plate_class_output_path.filename().string()) << "\"";
    output << ",\n    \"final_crust_age_png\": \"" << json_escape(run_id + "_CA.png") << "\"";
    output << ",\n    \"final_mineral_png\": \"" << json_escape(run_id + "_MN.png") << "\"";
    output << ",\n    \"final_tension_png\": \"" << json_escape(run_id + "_TN.png") << "\"";
    if (final_boundary_output_path != nullptr) {
        output << ",\n    \"final_boundary_png\": \""
               << json_escape(final_boundary_output_path->filename().string()) << "\"";
    }
    if (animation_output_path != nullptr) {
        output << ",\n    \"height_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::Height)) << "\"";
        output << ",\n    \"material_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::Material)) << "\"";
        output << ",\n    \"crust_type_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::PlateType)) << "\"";
        output << ",\n    \"id_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::PlateId)) << "\"";
        output << ",\n    \"boundary_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::Boundary)) << "\"";
        output << ",\n    \"crust_age_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::CrustAge)) << "\"";
        output << ",\n    \"mineral_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::Mineral)) << "\"";
        output << ",\n    \"tension_animation_png\": \""
               << json_escape(animation_file_name(run_id, MapArtifact::Tension)) << "\"";
    }
    output << "\n  },\n";
    output << "  \"checkpoints\": [\n";
    for (size_t i = 0; i < checkpoints.size(); ++i) {
        const CheckpointRecord& checkpoint = checkpoints[i];
        output << "    {\n";
        output << "      \"label\": \"" << json_escape(checkpoint.label) << "\",\n";
        output << "      \"step\": " << checkpoint.step << ",\n";
        output << "      \"iteration_count\": " << checkpoint.iteration_count << ",\n";
        output << "      \"cycle_count\": " << checkpoint.cycle_count << ",\n";
        output << "      \"time_myr\": " << checkpoint.time_myr;
        if (!checkpoint.frame_file.empty()) {
            output << ",\n      \"frame\": \"" << json_escape(checkpoint.frame_file) << "\"\n";
        } else {
            output << "\n";
        }
        output << "    }" << (i + 1 < checkpoints.size() ? "," : "") << "\n";
    }
    output << "  ]\n";
    output << "}\n";
    if (!output) {
        fail("failed to write run manifest: " + display_path(path));
    }
}

void save_image(const float* heightmap, const fs::path& filename, int width, int height, bool colors)
{
    std::vector<float> copy(heightmap,
                            heightmap + static_cast<size_t>(width) * static_cast<size_t>(height));

    const int result = colors
        ? writeImageColors(filename.string().c_str(), width, height, copy.data(), "Plate Tectonics")
        : writeImageGray(filename.string().c_str(), width, height, copy.data(), "Plate Tectonics");

    if (result != 0) {
        fail("failed to write image: " + display_path(filename));
    }
}

enum class BoundaryType {
    None,
    Convergent,
    Divergent,
    Transform,
    PassiveMargin
};

void apply_boundary_color(std::vector<png_byte>& rgb, size_t pixel_index, BoundaryType boundary_type)
{
    png_byte* ptr = &rgb[pixel_index * 3U];
    switch (boundary_type) {
    case BoundaryType::Convergent:
        ptr[0] = 255;
        ptr[1] = 64;
        ptr[2] = 64;
        break;
    case BoundaryType::Divergent:
        ptr[0] = 64;
        ptr[1] = 192;
        ptr[2] = 255;
        break;
    case BoundaryType::Transform:
        ptr[0] = 255;
        ptr[1] = 215;
        ptr[2] = 0;
        break;
    case BoundaryType::PassiveMargin:
        ptr[0] = 196;
        ptr[1] = 160;
        ptr[2] = 96;
        break;
    case BoundaryType::None:
    default:
        break;
    }
}

bool capture_boundary_overlay(void* simulation, int width, int height, BoundaryOverlayData& data)
{
    const size_t cell_count = static_cast<size_t>(width) * static_cast<size_t>(height);
    const uint8_t* boundary_type_map = platec_api_get_boundary_type_map(simulation);
    if (boundary_type_map == nullptr) {
        return false;
    }

    data.boundary_type.assign(boundary_type_map, boundary_type_map + cell_count);
    return true;
}

void overlay_boundaries(const BoundaryOverlayData& data, int width, int height, std::vector<png_byte>& rgb)
{
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index = static_cast<size_t>(y) * static_cast<size_t>(width) +
                                       static_cast<size_t>(x);
            BoundaryType boundary_type = BoundaryType::None;
            switch (static_cast<platec::contract::BoundaryType>(data.boundary_type[pixel_index])) {
            case platec::contract::BoundaryType::Convergent:
                boundary_type = BoundaryType::Convergent;
                break;
            case platec::contract::BoundaryType::Divergent:
                boundary_type = BoundaryType::Divergent;
                break;
            case platec::contract::BoundaryType::Transform:
                boundary_type = BoundaryType::Transform;
                break;
            case platec::contract::BoundaryType::PassiveMargin:
                boundary_type = BoundaryType::PassiveMargin;
                break;
            case platec::contract::BoundaryType::None:
            default:
                boundary_type = BoundaryType::None;
                break;
            }

            apply_boundary_color(rgb, pixel_index, boundary_type);
        }
    }
}

void render_boundary_overlay_rgba(const BoundaryOverlayData& data, int width, int height,
                                  std::vector<png_byte>& rgba)
{
    rgba.assign(static_cast<size_t>(width) * static_cast<size_t>(height) * 4U, 0);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            const size_t pixel_index = static_cast<size_t>(y) * static_cast<size_t>(width) +
                                       static_cast<size_t>(x);
            BoundaryType boundary_type = BoundaryType::None;
            switch (static_cast<platec::contract::BoundaryType>(data.boundary_type[pixel_index])) {
            case platec::contract::BoundaryType::Convergent:
                boundary_type = BoundaryType::Convergent;
                break;
            case platec::contract::BoundaryType::Divergent:
                boundary_type = BoundaryType::Divergent;
                break;
            case platec::contract::BoundaryType::Transform:
                boundary_type = BoundaryType::Transform;
                break;
            case platec::contract::BoundaryType::PassiveMargin:
                boundary_type = BoundaryType::PassiveMargin;
                break;
            case platec::contract::BoundaryType::None:
            default:
                boundary_type = BoundaryType::None;
                break;
            }

            if (boundary_type == BoundaryType::None) {
                continue;
            }

            const size_t base = pixel_index * 4U;
            png_byte* ptr = &rgba[base];
            ptr[3] = 255;
            switch (boundary_type) {
            case BoundaryType::Convergent:
                ptr[0] = 255;
                ptr[1] = 64;
                ptr[2] = 64;
                break;
            case BoundaryType::Divergent:
                ptr[0] = 64;
                ptr[1] = 192;
                ptr[2] = 255;
                break;
            case BoundaryType::Transform:
                ptr[0] = 255;
                ptr[1] = 215;
                ptr[2] = 0;
                break;
            case BoundaryType::PassiveMargin:
                ptr[0] = 196;
                ptr[1] = 160;
                ptr[2] = 96;
                break;
            case BoundaryType::None:
            default:
                break;
            }
        }
    }
}

void save_boundary_overlay(const fs::path& filename, int width, int height,
                           const BoundaryOverlayData& boundaries)
{
    std::vector<png_byte> rgba;
    render_boundary_overlay_rgba(boundaries, width, height, rgba);
    if (writeImageRgba(filename.string().c_str(), width, height, rgba.data(),
                       "Plate Boundaries") != 0) {
        fail("failed to write boundary overlay: " + display_path(filename));
    }
}

void save_image(void* simulation, const fs::path& filename, int width, int height, bool colors,
                bool show_boundaries,
                const BoundaryOverlayData* fallback_boundaries = nullptr)
{
    if (!show_boundaries) {
        save_image(platec_api_get_heightmap(simulation), filename, width, height, colors);
        return;
    }

    std::vector<float> copy(platec_api_get_heightmap(simulation),
                            platec_api_get_heightmap(simulation) +
                                static_cast<size_t>(width) * static_cast<size_t>(height));

    std::vector<png_byte> rgb;
    if (colors) {
        renderImageColorsRgb(width, height, copy.data(), rgb);
    } else {
        renderImageGrayRgb(width, height, copy.data(), rgb);
    }

    BoundaryOverlayData current_boundaries;
    if (capture_boundary_overlay(simulation, width, height, current_boundaries)) {
        overlay_boundaries(current_boundaries, width, height, rgb);
    } else if (fallback_boundaries != nullptr && fallback_boundaries->valid()) {
        overlay_boundaries(*fallback_boundaries, width, height, rgb);
    }

    if (writeImageRgb(filename.string().c_str(), width, height, rgb.data(), "Plate Tectonics") != 0) {
        fail("failed to write image: " + display_path(filename));
    }
}

void save_map_png(void* simulation, const Params& params, MapArtifact artifact,
                  const fs::path& filename, int width, int height, bool colors, bool show_boundaries,
                  const BoundaryOverlayData* fallback_boundaries = nullptr)
{
    switch (artifact) {
    case MapArtifact::Height:
        save_image(simulation, filename, width, height, colors, show_boundaries,
                   fallback_boundaries);
        return;
    case MapArtifact::Material:
        export_material_map(platec_api_get_material_map(simulation), width, height, filename);
        return;
    case MapArtifact::PlateType:
        export_crust_type_map(simulation, width, height, filename, show_boundaries);
        return;
    case MapArtifact::PlateId:
        export_plate_id_map(platec_api_get_platesmap(simulation), width, height, filename);
        return;
    case MapArtifact::Boundary:
        export_boundary_map(simulation, width, height, filename, fallback_boundaries);
        return;
    case MapArtifact::CrustAge:
        export_crust_age_map(platec_api_get_crust_age_myr_map(simulation), width, height, filename);
        return;
    case MapArtifact::Mineral:
        export_mineral_map(simulation, params, width, height, filename);
        return;
    case MapArtifact::Tension:
        export_tension_map(platec_api_get_tension_map(simulation), width, height, filename);
        return;
    default:
        fail("unsupported map artifact");
    }
}

struct RgbFrame {
    uint32_t width = 0;
    uint32_t height = 0;
    std::vector<png_byte> rgb;
};

void append_u16_be(std::vector<png_byte>& bytes, uint16_t value)
{
    bytes.push_back(static_cast<png_byte>((value >> 8U) & 0xFFU));
    bytes.push_back(static_cast<png_byte>(value & 0xFFU));
}

void append_u32_be(std::vector<png_byte>& bytes, uint32_t value)
{
    bytes.push_back(static_cast<png_byte>((value >> 24U) & 0xFFU));
    bytes.push_back(static_cast<png_byte>((value >> 16U) & 0xFFU));
    bytes.push_back(static_cast<png_byte>((value >> 8U) & 0xFFU));
    bytes.push_back(static_cast<png_byte>(value & 0xFFU));
}

void write_bytes(FILE* fp, const void* data, size_t byte_count, const fs::path& path)
{
    if (byte_count == 0U) {
        return;
    }
    if (fwrite(data, 1U, byte_count, fp) != byte_count) {
        fail("failed to write animated PNG: " + display_path(path));
    }
}

void write_u32_be(FILE* fp, uint32_t value, const fs::path& path)
{
    const png_byte bytes[] = {
        static_cast<png_byte>((value >> 24U) & 0xFFU),
        static_cast<png_byte>((value >> 16U) & 0xFFU),
        static_cast<png_byte>((value >> 8U) & 0xFFU),
        static_cast<png_byte>(value & 0xFFU),
    };
    write_bytes(fp, bytes, sizeof(bytes), path);
}

void write_png_chunk(FILE* fp, const char* type, const std::vector<png_byte>& payload,
                     const fs::path& path)
{
    if (payload.size() > static_cast<size_t>(std::numeric_limits<uInt>::max())) {
        fail("animated PNG chunk is too large: " + display_path(path));
    }

    write_u32_be(fp, static_cast<uint32_t>(payload.size()), path);
    write_bytes(fp, type, 4U, path);
    write_bytes(fp, payload.data(), payload.size(), path);

    uLong crc = crc32(0L, Z_NULL, 0U);
    crc = crc32(crc, reinterpret_cast<const Bytef*>(type), 4U);
    if (!payload.empty()) {
        crc = crc32(crc, payload.data(), static_cast<uInt>(payload.size()));
    }
    write_u32_be(fp, static_cast<uint32_t>(crc), path);
}

RgbFrame read_png_rgb_frame(const fs::path& path)
{
    RgbFrame frame;
    FILE* fp = nullptr;
#ifdef _WIN32
    FILE* fp_temp = nullptr;
    const errno_t err = fopen_s(&fp_temp, path.string().c_str(), "rb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(path.string().c_str(), "rb");
    if (fp == nullptr) {
#endif
        fail("failed to open PNG frame: " + display_path(path));
    }

    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (png_ptr == nullptr) {
        fclose(fp);
        fail("failed to allocate PNG reader");
    }

    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == nullptr) {
        png_destroy_read_struct(&png_ptr, nullptr, nullptr);
        fclose(fp);
        fail("failed to allocate PNG info");
    }

    png_bytep raw = nullptr;
    png_bytepp rows = nullptr;
    const auto cleanup = [&]() {
        if (raw != nullptr) {
            free(raw);
        }
        if (rows != nullptr) {
            free(rows);
        }
        if (png_ptr != nullptr) {
            png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
        }
        if (fp != nullptr) {
            fclose(fp);
        }
    };

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4611)
#endif
    if (setjmp(png_jmpbuf(png_ptr))) {
        cleanup();
        fail("failed to read PNG frame: " + display_path(path));
    }
#ifdef _MSC_VER
#pragma warning(pop)
#endif

    png_init_io(png_ptr, fp);
    png_read_info(png_ptr, info_ptr);

    png_uint_32 png_width = 0;
    png_uint_32 png_height = 0;
    int bit_depth = 0;
    int color_type = 0;
    int interlace_method = 0;
    png_get_IHDR(png_ptr, info_ptr, &png_width, &png_height, &bit_depth, &color_type,
                 &interlace_method, nullptr, nullptr);
    if (png_width == 0U || png_height == 0U) {
        cleanup();
        fail("PNG frame has invalid dimensions: " + display_path(path));
    }

    if (bit_depth == 16) {
        png_set_strip_16(png_ptr);
    }
    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png_ptr);
    }
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(png_ptr);
    }
    if (png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS) != 0) {
        png_set_tRNS_to_alpha(png_ptr);
    }
    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(png_ptr);
    }
    png_set_strip_alpha(png_ptr);
    png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);

    png_get_IHDR(png_ptr, info_ptr, &png_width, &png_height, &bit_depth, &color_type,
                 &interlace_method, nullptr, nullptr);
    if (bit_depth != 8 || color_type != PNG_COLOR_TYPE_RGB) {
        cleanup();
        fail("PNG frame could not be converted to RGB8: " + display_path(path));
    }

    const png_size_t row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    const size_t expected_row_bytes = static_cast<size_t>(png_width) * 3U;
    if (row_bytes != expected_row_bytes) {
        cleanup();
        fail("PNG frame row stride is unsupported: " + display_path(path));
    }
    if (expected_row_bytes > std::numeric_limits<size_t>::max() / static_cast<size_t>(png_height)) {
        cleanup();
        fail("PNG frame is too large: " + display_path(path));
    }

    const size_t raw_size = expected_row_bytes * static_cast<size_t>(png_height);
    raw = static_cast<png_bytep>(malloc(raw_size));
    rows = static_cast<png_bytepp>(malloc(static_cast<size_t>(png_height) * sizeof(png_bytep)));
    if (raw == nullptr || rows == nullptr) {
        cleanup();
        fail("failed to allocate PNG frame buffer");
    }

    for (png_uint_32 y = 0; y < png_height; ++y) {
        rows[y] = raw + static_cast<size_t>(y) * expected_row_bytes;
    }

    png_read_image(png_ptr, rows);
    png_read_end(png_ptr, nullptr);

    frame.width = static_cast<uint32_t>(png_width);
    frame.height = static_cast<uint32_t>(png_height);
    frame.rgb.assign(raw, raw + raw_size);
    cleanup();
    return frame;
}

std::vector<png_byte> compress_png_rgb_frame(const RgbFrame& frame, const fs::path& path)
{
    const size_t row_bytes = static_cast<size_t>(frame.width) * 3U;
    const size_t expected_size = row_bytes * static_cast<size_t>(frame.height);
    if (frame.rgb.size() != expected_size) {
        fail("animated PNG frame has invalid RGB data: " + display_path(path));
    }
    if (row_bytes + 1U < row_bytes ||
        row_bytes + 1U > std::numeric_limits<size_t>::max() / static_cast<size_t>(frame.height)) {
        fail("animated PNG frame is too large: " + display_path(path));
    }

    const size_t filtered_row_bytes = row_bytes + 1U;
    std::vector<png_byte> filtered(filtered_row_bytes * static_cast<size_t>(frame.height));
    for (uint32_t y = 0; y < frame.height; ++y) {
        png_byte* filtered_row = filtered.data() + static_cast<size_t>(y) * filtered_row_bytes;
        filtered_row[0] = 0U;
        memcpy(filtered_row + 1U, frame.rgb.data() + static_cast<size_t>(y) * row_bytes,
               row_bytes);
    }

    if (filtered.size() > static_cast<size_t>(std::numeric_limits<uLong>::max())) {
        fail("animated PNG frame is too large to compress: " + display_path(path));
    }

    const uLong source_size = static_cast<uLong>(filtered.size());
    uLongf compressed_size = compressBound(source_size);
    std::vector<png_byte> compressed(static_cast<size_t>(compressed_size));
    const int rc = compress2(compressed.data(), &compressed_size, filtered.data(), source_size,
                             Z_BEST_COMPRESSION);
    if (rc != Z_OK) {
        fail("failed to compress animated PNG frame: " + display_path(path));
    }
    compressed.resize(static_cast<size_t>(compressed_size));
    return compressed;
}

void write_frame_control_chunk(FILE* fp, const fs::path& output_path, uint32_t sequence_number,
                               uint32_t width, uint32_t height)
{
    std::vector<png_byte> frame_control;
    frame_control.reserve(26U);
    append_u32_be(frame_control, sequence_number);
    append_u32_be(frame_control, width);
    append_u32_be(frame_control, height);
    append_u32_be(frame_control, 0U);
    append_u32_be(frame_control, 0U);
    append_u16_be(frame_control, kAnimationDelayCs);
    append_u16_be(frame_control, 100U);
    frame_control.push_back(0U);
    frame_control.push_back(0U);
    write_png_chunk(fp, "fcTL", frame_control, output_path);
}

void create_animated_png(const fs::path& output_path, const std::vector<fs::path>& frame_paths)
{
    if (frame_paths.empty()) {
        fail("cannot create animated PNG without frames");
    }
    if (frame_paths.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
        fail("animated PNG has too many frames: " + display_path(output_path));
    }

    FILE* fp = nullptr;
#ifdef _WIN32
    FILE* fp_temp = nullptr;
    const errno_t err = fopen_s(&fp_temp, output_path.string().c_str(), "wb");
    fp = fp_temp;
    if (err != 0 || fp == nullptr) {
#else
    fp = fopen(output_path.string().c_str(), "wb");
    if (fp == nullptr) {
#endif
        fail("failed to open animated PNG for writing: " + display_path(output_path));
    }

    const png_byte signature[] = {0x89U, 0x50U, 0x4EU, 0x47U, 0x0DU, 0x0AU, 0x1AU, 0x0AU};
    write_bytes(fp, signature, sizeof(signature), output_path);

    const RgbFrame first_frame = read_png_rgb_frame(frame_paths[0]);

    std::vector<png_byte> ihdr;
    ihdr.reserve(13U);
    append_u32_be(ihdr, first_frame.width);
    append_u32_be(ihdr, first_frame.height);
    ihdr.push_back(8U);
    ihdr.push_back(PNG_COLOR_TYPE_RGB);
    ihdr.push_back(PNG_COMPRESSION_TYPE_BASE);
    ihdr.push_back(PNG_FILTER_TYPE_BASE);
    ihdr.push_back(PNG_INTERLACE_NONE);
    write_png_chunk(fp, "IHDR", ihdr, output_path);

    std::vector<png_byte> animation_control;
    animation_control.reserve(8U);
    append_u32_be(animation_control, static_cast<uint32_t>(frame_paths.size()));
    append_u32_be(animation_control, 0U);
    write_png_chunk(fp, "acTL", animation_control, output_path);

    uint32_t sequence_number = 0;
    write_frame_control_chunk(fp, output_path, sequence_number++, first_frame.width,
                              first_frame.height);
    const std::vector<png_byte> first_image_data =
        compress_png_rgb_frame(first_frame, frame_paths[0]);
    write_png_chunk(fp, "IDAT", first_image_data, output_path);

    for (size_t frame_index = 1; frame_index < frame_paths.size(); ++frame_index) {
        const RgbFrame frame = read_png_rgb_frame(frame_paths[frame_index]);
        if (frame.width != first_frame.width || frame.height != first_frame.height) {
            fclose(fp);
            fail("animated PNG frame dimensions do not match: " +
                 display_path(frame_paths[frame_index]));
        }

        write_frame_control_chunk(fp, output_path, sequence_number++, frame.width, frame.height);
        const std::vector<png_byte> compressed =
            compress_png_rgb_frame(frame, frame_paths[frame_index]);
        if (compressed.size() >
            static_cast<size_t>(std::numeric_limits<uint32_t>::max()) - 4U) {
            fclose(fp);
            fail("animated PNG frame chunk is too large: " + display_path(frame_paths[frame_index]));
        }

        std::vector<png_byte> frame_data;
        frame_data.reserve(compressed.size() + 4U);
        append_u32_be(frame_data, sequence_number++);
        frame_data.insert(frame_data.end(), compressed.begin(), compressed.end());
        write_png_chunk(fp, "fdAT", frame_data, output_path);
    }

    write_png_chunk(fp, "IEND", std::vector<png_byte>(), output_path);
    if (fclose(fp) != 0) {
        fail("failed to close animated PNG: " + display_path(output_path));
    }
}

void delete_files(const std::vector<fs::path>& paths)
{
    for (const fs::path& path : paths) {
        std::error_code ec;
        const bool removed = fs::remove(path, ec);
        if (!removed || ec) {
            fail("failed to delete frame: " + display_path(path));
        }
    }
}

} // namespace

int main(int argc, char* argv[])
{
    ensure_console_output();

    Params params = fill_params(argc, argv);
    std::vector<std::vector<fs::path>> animation_frames(kMapArtifactCount);
    std::vector<std::vector<fs::path>> animation_temp_files(kMapArtifactCount);
    std::vector<fs::path> animation_output_paths(kMapArtifactCount);
    std::vector<fs::path> step_outputs;
    std::vector<CheckpointRecord> checkpoints;
    BoundaryOverlayData last_boundary_state;
    fs::path initial_image_output_path;
    fs::path initial_material_output_path;
    fs::path initial_plate_output_path;
    fs::path initial_plate_class_output_path;
    fs::path initial_crust_age_output_path;
    fs::path initial_mineral_output_path;
    fs::path initial_tension_output_path;
    fs::path final_image_output_path;
    fs::path final_material_output_path;
    fs::path final_plate_output_path;
    fs::path final_plate_class_output_path;
    fs::path final_crust_age_output_path;
    fs::path final_mineral_output_path;
    fs::path final_tension_output_path;
    fs::path initial_boundary_output_path;
    fs::path final_boundary_output_path;
    fs::path run_manifest_output_path;
    std::string initial_image_label;
    std::string initial_material_label;
    std::string initial_plate_label;
    std::string initial_plate_class_label;
    std::string initial_crust_age_label;
    std::string initial_mineral_label;
    std::string initial_tension_label;
    std::string final_image_label;
    std::string final_material_label;
    std::string final_plate_label;
    std::string final_plate_class_label;
    std::string final_crust_age_label;
    std::string final_mineral_label;
    std::string final_tension_label;
    std::vector<std::string> animation_labels(kMapArtifactCount);
    std::string run_id;

    if (!params.aggregation_overlap_abs_explicit) {
        params.aggregation_overlap_abs =
            default_aggregation_overlap_abs(params.width, params.height);
    }

    ensure_output_directories();
    const std::time_t started_at = std::time(nullptr);
    run_id = make_run_id(params, started_at);
    initial_image_output_path = kDefaultInputDir / (run_id + "_HT.tiff");
    initial_material_output_path = kDefaultInputDir / (run_id + "_MT.png");
    initial_plate_output_path = kDefaultInputDir / (run_id + "_ID.png");
    initial_plate_class_output_path = kDefaultInputDir / (run_id + "_PT.png");
    initial_crust_age_output_path = kDefaultInputDir / (run_id + "_CA.png");
    initial_mineral_output_path = kDefaultInputDir / (run_id + "_MN.png");
    initial_tension_output_path = kDefaultInputDir / (run_id + "_TN.png");
    final_image_output_path = kDefaultOutputDir / (run_id + "_HT.tiff");
    final_material_output_path = kDefaultOutputDir / (run_id + "_MT.png");
    final_plate_output_path = kDefaultOutputDir / (run_id + "_ID.png");
    final_plate_class_output_path = kDefaultOutputDir / (run_id + "_PT.png");
    final_crust_age_output_path = kDefaultOutputDir / (run_id + "_CA.png");
    final_mineral_output_path = kDefaultOutputDir / (run_id + "_MN.png");
    final_tension_output_path = kDefaultOutputDir / (run_id + "_TN.png");
    initial_boundary_output_path = kDefaultInputDir / (run_id + "_B.png");
    final_boundary_output_path = kDefaultOutputDir / (run_id + "_B.png");
    run_manifest_output_path = kDefaultOutputDir / (run_id + ".run.json");
    for (size_t map_index = 0; map_index < kMapArtifactCount; ++map_index) {
        const MapArtifact artifact = static_cast<MapArtifact>(map_index);
        animation_output_paths[map_index] =
            kDefaultAnimationDir / animation_file_name(run_id, artifact);
        animation_labels[map_index] = display_path(animation_output_paths[map_index]);
    }
    initial_image_label = display_path(initial_image_output_path);
    initial_material_label = display_path(initial_material_output_path);
    initial_plate_label = display_path(initial_plate_output_path);
    initial_plate_class_label = display_path(initial_plate_class_output_path);
    initial_crust_age_label = display_path(initial_crust_age_output_path);
    initial_mineral_label = display_path(initial_mineral_output_path);
    initial_tension_label = display_path(initial_tension_output_path);
    final_image_label = display_path(final_image_output_path);
    final_material_label = display_path(final_material_output_path);
    final_plate_label = display_path(final_plate_output_path);
    final_plate_class_label = display_path(final_plate_class_output_path);
    final_crust_age_label = display_path(final_crust_age_output_path);
    final_mineral_label = display_path(final_mineral_output_path);
    final_tension_label = display_path(final_tension_output_path);

    const int32_t sea_level_override_m =
        params.has_sea_level_m ? static_cast<int32_t>(params.sea_level_m)
                               : TopographyCodec::kNoSeaLevelOverride;

    platec::scenario::Scenario scenario;
    scenario.seed = params.seed;
    scenario.width = params.width;
    scenario.height = params.height;
    scenario.sea_level = kSeaLevel;
    scenario.erosion_period = params.erosion_period;
    scenario.folding_ratio = params.folding_ratio;
    scenario.aggregation_overlap_abs = params.aggregation_overlap_abs;
    scenario.aggregation_overlap_rel = params.aggregation_overlap_rel;
    scenario.cycle_count = params.cycles;
    scenario.plate_count = params.plates;
    scenario.erosion_strength = params.erosion_strength;
    scenario.crust_rotation_strength = params.landmass_rotation;
    scenario.rotation_strength = params.rotation_strength;
    scenario.subduction_strength = params.subduction_strength;
    scenario.sea_level_m = sea_level_override_m;
    scenario.crust_type_boundary_m = params.crust_type_boundary_m;
    scenario.cycle_step_limit = params.cycle_steps;
    scenario.cycle_duration_myr = params.cycle_duration_myr;
    scenario.divergent_carve_strength = params.divergent_carve;
    scenario.delta_time_myr = params.delta_time_myr;
    scenario.gravity_mps2 = params.gravity_mps2;
    scenario.glacial_erosion_strength = params.glacial_erosion_strength;
    scenario.glacial_erosion_period = params.glacial_erosion_period;
    scenario.collision_uplift_ratio = params.collision_uplift_ratio;
    scenario.continental_compression_gain = params.continental_compression_gain;
    scenario.continental_boundary_fluidity = params.continental_boundary_fluidity;
    scenario.movement_energy = params.movement_energy;
    scenario.initial_max_height_m = params.max_initial_height_m;
    scenario.hf_noise_period = params.hf_noise_period;
    scenario.hf_noise_strength = params.hf_noise_strength;
    scenario.lf_noise_period = params.lf_noise_period;
    scenario.lf_noise_strength = params.lf_noise_strength;

    void* simulation = platec_api_create_from_scenario(&scenario);
    if (simulation == nullptr) {
        fail("failed to create simulation");
    }

    const uint16_t active_sea_level_m = platec_api_get_sea_level_m(simulation);
    const auto elapsed_updates = [&]() -> uint32_t {
        return static_cast<uint32_t>(std::llround(
            platec_api_get_time_myr(simulation) / params.delta_time_myr));
    };

    printf("Plate-tectonics simulation example\n");
    printf(" seed     : %u\n", params.seed);
    printf(" width    : %u\n", params.width);
    printf(" height   : %u\n", params.height);
    printf(" preview  : %s\n", params.colors ? "color" : "grayscale");
    printf(" tiff16   : grayscale uint16 metric elevation\n");
    printf(" cycles   : %u\n", params.cycles);
    printf(" cyclemax : %u\n", params.cycle_steps);
    printf(" plates   : %u\n", params.plates);
    printf(" aggregate: abs %u, rel %.3f\n", params.aggregation_overlap_abs,
           params.aggregation_overlap_rel);
    printf(" folding  : %.3f\n", params.folding_ratio);
    printf(" wind     : every %u updates, strength %.2f\n", params.erosion_period,
           params.erosion_strength);
    printf(" rotation : motion %.2f, landmass %.2f\n", params.rotation_strength,
           params.landmass_rotation);
    printf(" subduct  : %.2f\n", params.subduction_strength);
    printf(" diverge  : %.3f\n", params.divergent_carve);
    printf(" delta t  : %.3f Myr/update\n", params.delta_time_myr);
    if (params.cycle_duration_myr > 0.0) {
        printf(" cycle dur: %.3f Myr (%u steps)\n", params.cycle_duration_myr, params.cycle_steps);
    }
    printf(" gravity  : %.5f m/s^2\n", params.gravity_mps2);
    printf(" glacial  : every %u updates, strength %.2f\n", params.glacial_erosion_period,
           params.glacial_erosion_strength);
    printf(" collide q: %.3f\n", params.collision_uplift_ratio);
    printf(" cc gain  : %.3f\n", params.continental_compression_gain);
    printf(" cc fluid : %.3f\n", params.continental_boundary_fluidity);
    printf(" motion e : %.3f\n", params.movement_energy);
    if (params.crust_type_boundary_m != platec::scenario::kAutoCrustTypeBoundaryMeters) {
        printf(" crust ct : %d m\n", params.crust_type_boundary_m);
    }
    printf(" init max : %u\n", params.max_initial_height_m);
    printf(" hf noise : every %u, strength %.2f\n", params.hf_noise_period,
           params.hf_noise_strength);
    printf(" lf noise : every %u, strength %.2f\n", params.lf_noise_period,
           params.lf_noise_strength);
    printf(" mineral  : T0 %.1fC, dT %.1fC/km, dP %.2fkbar/km, hydro %.2f, fluids %.2f, meta %.2f, noise %.2f\n",
           params.mineral_surface_temp_c, params.mineral_geothermal_gradient_c_per_km,
           params.mineral_pressure_gradient_kbar_per_km, params.mineral_hydrothermal_activity,
           params.mineral_fluid_flux, params.mineral_metamorphism_gain, params.mineral_noise);
    printf(" sea lvl  : %u m\n", active_sea_level_m);
    printf(" init tiff16: %s\n", initial_image_label.c_str());
    printf(" init material: %s\n", initial_material_label.c_str());
    printf(" init plate ids: %s\n", initial_plate_label.c_str());
    printf(" init crust types: %s\n", initial_plate_class_label.c_str());
    printf(" init crust age: %s\n", initial_crust_age_label.c_str());
    printf(" init mineral : %s\n", initial_mineral_label.c_str());
    printf(" init tension : %s\n", initial_tension_label.c_str());
    printf(" init boundaries: %s\n", display_path(initial_boundary_output_path).c_str());
    printf(" final tiff16: %s\n", final_image_label.c_str());
    printf(" final material: %s\n", final_material_label.c_str());
    printf(" final plate ids: %s\n", final_plate_label.c_str());
    printf(" final crust types: %s\n", final_plate_class_label.c_str());
    printf(" final crust age: %s\n", final_crust_age_label.c_str());
    printf(" final mineral : %s\n", final_mineral_label.c_str());
    printf(" final tension : %s\n", final_tension_label.c_str());
    printf(" final boundaries: %s\n", display_path(final_boundary_output_path).c_str());
    if (params.animated_step != 0) {
        printf(" anim step: %u\n", params.animated_step);
        for (size_t map_index = 0; map_index < kMapArtifactCount; ++map_index) {
            const MapArtifact artifact = static_cast<MapArtifact>(map_index);
            printf(" anim %-3s: %s\n", map_suffix(artifact),
                   animation_labels[map_index].c_str());
        }
    }
    if (params.show_boundaries) {
        printf(" bounds   : convergent=red divergent=blue transform=yellow passive=brown\n");
    }
    if (params.step == 0) {
        printf(" step     : no\n");
    } else {
        printf(" step     : %u\n", params.step);
    }
    if (params.step != 0 || params.animated_step != 0) {
        printf(" frames   : %s\n", display_path(kDefaultFramesDir).c_str());
    }
    printf("\n");

    capture_boundary_overlay(simulation, static_cast<int>(params.width),
                             static_cast<int>(params.height), last_boundary_state);

    const bool overlay_boundaries_on_exports =
        params.show_boundaries && params.animated_step != 0;

    export_metric_heightmap(platec_api_get_heightmap(simulation), static_cast<int>(params.width),
                            static_cast<int>(params.height), active_sea_level_m,
                            initial_image_output_path);
    export_material_map(platec_api_get_material_map(simulation), static_cast<int>(params.width),
                        static_cast<int>(params.height), initial_material_output_path);
    export_plate_id_map(platec_api_get_platesmap(simulation), static_cast<int>(params.width),
                        static_cast<int>(params.height), initial_plate_output_path);
    export_crust_type_map(simulation, static_cast<int>(params.width), static_cast<int>(params.height),
                          initial_plate_class_output_path, overlay_boundaries_on_exports);
    export_crust_age_map(platec_api_get_crust_age_myr_map(simulation), static_cast<int>(params.width),
                         static_cast<int>(params.height), initial_crust_age_output_path);
    export_mineral_map(simulation, params, static_cast<int>(params.width),
                       static_cast<int>(params.height), initial_mineral_output_path);
    export_tension_map(platec_api_get_tension_map(simulation), static_cast<int>(params.width),
                       static_cast<int>(params.height), initial_tension_output_path);
    export_boundary_map(simulation, static_cast<int>(params.width), static_cast<int>(params.height),
                        initial_boundary_output_path, &last_boundary_state);
    checkpoints.push_back(
        CheckpointRecord{"initial", elapsed_updates(), platec_api_get_iteration_count(simulation),
                         platec_api_get_cycle_count(simulation),
                         platec_api_get_time_myr(simulation), ""});

    uint32_t saved_step_frame_count = 0;
    uint32_t saved_animation_frame_count = 0;
    const auto save_frame_set = [&](const std::string& prefix, uint32_t& frame_index,
                                    const BoundaryOverlayData* fallback_boundaries,
                                    std::vector<std::vector<fs::path>>* animation_frame_store =
                                        nullptr,
                                    std::vector<std::vector<fs::path>>* animation_temp_store =
                                        nullptr) {
        fs::path height_frame_path;
        const uint32_t current_frame = frame_index;
        for (size_t map_index = 0; map_index < kMapArtifactCount; ++map_index) {
            const MapArtifact artifact = static_cast<MapArtifact>(map_index);
            const fs::path frame_path =
                kDefaultFramesDir / frame_file_name(prefix, current_frame, artifact);
            save_map_png(simulation, params, artifact, frame_path, static_cast<int>(params.width),
                         static_cast<int>(params.height), params.colors,
                         overlay_boundaries_on_exports,
                         fallback_boundaries);
            if (artifact == MapArtifact::Height) {
                height_frame_path = frame_path;
            }
            if (animation_frame_store != nullptr) {
                (*animation_frame_store)[map_index].push_back(frame_path);
            }
            if (animation_temp_store != nullptr) {
                (*animation_temp_store)[map_index].push_back(frame_path);
            }
        }
        ++frame_index;
        return height_frame_path;
    };

    if (params.animated_step != 0) {
        save_frame_set(run_id + "_anim", saved_animation_frame_count, &last_boundary_state,
                       &animation_frames, &animation_temp_files);
    }

    uint32_t step = 0;
    while (platec_api_is_finished(simulation) == 0) {
        ++step;
        const BoundaryOverlayData boundary_before_step = last_boundary_state;
        platec_api_step(simulation);
        const bool finished_after_step = platec_api_is_finished(simulation) != 0;

        if (finished_after_step) {
            break;
        }

        if (params.animated_step != 0 && step % params.animated_step == 0) {
            save_frame_set(run_id + "_anim", saved_animation_frame_count, &boundary_before_step,
                           &animation_frames, &animation_temp_files);
        }

        if (params.step != 0 && step % params.step == 0) {
            const fs::path step_output_path =
                save_frame_set(run_id, saved_step_frame_count, &boundary_before_step);
            step_outputs.push_back(step_output_path);
            checkpoints.push_back(
                CheckpointRecord{"sample", elapsed_updates(),
                                 platec_api_get_iteration_count(simulation),
                                 platec_api_get_cycle_count(simulation),
                                 platec_api_get_time_myr(simulation),
                                 step_output_path.filename().string()});
            printf(" * frame %u (step %u, filename %s)\n",
                   saved_step_frame_count - 1, step, display_path(step_output_path).c_str());
        }

        capture_boundary_overlay(simulation, static_cast<int>(params.width),
                                 static_cast<int>(params.height), last_boundary_state);
    }

    export_metric_heightmap(platec_api_get_heightmap(simulation), static_cast<int>(params.width),
                            static_cast<int>(params.height), active_sea_level_m,
                            final_image_output_path);
    export_material_map(platec_api_get_material_map(simulation), static_cast<int>(params.width),
                        static_cast<int>(params.height), final_material_output_path);
    export_plate_id_map(platec_api_get_platesmap(simulation), static_cast<int>(params.width),
                        static_cast<int>(params.height), final_plate_output_path);
    export_crust_type_map(simulation, static_cast<int>(params.width), static_cast<int>(params.height),
                          final_plate_class_output_path, overlay_boundaries_on_exports);
    export_crust_age_map(platec_api_get_crust_age_myr_map(simulation), static_cast<int>(params.width),
                         static_cast<int>(params.height), final_crust_age_output_path);
    export_mineral_map(simulation, params, static_cast<int>(params.width),
                       static_cast<int>(params.height), final_mineral_output_path);
    export_tension_map(platec_api_get_tension_map(simulation), static_cast<int>(params.width),
                       static_cast<int>(params.height), final_tension_output_path);
    export_boundary_map(simulation, static_cast<int>(params.width), static_cast<int>(params.height),
                        final_boundary_output_path, &last_boundary_state);

    if (params.animated_step != 0) {
        if (animation_frames[0].empty() || step % params.animated_step != 0) {
            save_frame_set(run_id + "_anim", saved_animation_frame_count, &last_boundary_state,
                           &animation_frames, &animation_temp_files);
        }
        for (size_t map_index = 0; map_index < kMapArtifactCount; ++map_index) {
            create_animated_png(animation_output_paths[map_index], animation_frames[map_index]);
            delete_files(animation_temp_files[map_index]);
        }
    }

    checkpoints.push_back(
        CheckpointRecord{"final", elapsed_updates(),
                         platec_api_get_iteration_count(simulation),
                         platec_api_get_cycle_count(simulation),
                         platec_api_get_time_myr(simulation), ""});

    write_run_manifest(run_manifest_output_path, params, run_id,
                       active_sea_level_m, platec_api_get_time_origin_step(simulation),
                       platec_api_get_iteration_count(simulation),
                       platec_api_get_cycle_count(simulation),
                       platec_api_get_time_myr(simulation), checkpoints,
                       initial_image_output_path,
                       initial_material_output_path,
                       initial_plate_output_path,
                       initial_plate_class_output_path,
                       &initial_boundary_output_path,
                       final_image_output_path,
                       final_material_output_path,
                       final_plate_output_path,
                       final_plate_class_output_path,
                       &final_boundary_output_path,
                       params.animated_step != 0 ? &animation_output_paths[0] : nullptr);

    printf(" * initial metric TIFF16 exported (filename %s)\n", initial_image_label.c_str());
    printf(" * initial material map exported (filename %s)\n", initial_material_label.c_str());
    printf(" * initial plate ID map exported (filename %s)\n", initial_plate_label.c_str());
    printf(" * initial crust type map exported (filename %s)\n", initial_plate_class_label.c_str());
    printf(" * initial crust age map exported (filename %s)\n", initial_crust_age_label.c_str());
    printf(" * initial mineral map exported (filename %s)\n", initial_mineral_label.c_str());
    printf(" * initial tension map exported (filename %s)\n", initial_tension_label.c_str());
    printf(" * initial boundary map exported (filename %s)\n",
           display_path(initial_boundary_output_path).c_str());
    printf(" * final metric TIFF16 exported (filename %s)\n", final_image_label.c_str());
    printf(" * final material map exported (filename %s)\n", final_material_label.c_str());
    printf(" * final plate ID map exported (filename %s)\n", final_plate_label.c_str());
    printf(" * final crust type map exported (filename %s)\n", final_plate_class_label.c_str());
    printf(" * final crust age map exported (filename %s)\n", final_crust_age_label.c_str());
    printf(" * final mineral map exported (filename %s)\n", final_mineral_label.c_str());
    printf(" * final tension map exported (filename %s)\n", final_tension_label.c_str());
    printf(" * final boundary map exported (filename %s)\n",
           display_path(final_boundary_output_path).c_str());
    printf(" * run manifest exported (filename %s)\n",
           display_path(run_manifest_output_path).c_str());
    printf(" * elapsed geological time %.3f Myr\n", platec_api_get_time_myr(simulation));
    if (params.animated_step != 0) {
        for (size_t map_index = 0; map_index < kMapArtifactCount; ++map_index) {
            const MapArtifact artifact = static_cast<MapArtifact>(map_index);
            printf(" * %s animated PNG created (filename %s)\n", map_suffix(artifact),
                   animation_labels[map_index].c_str());
        }
    }

    platec_api_destroy(simulation);
    return 0;
}
