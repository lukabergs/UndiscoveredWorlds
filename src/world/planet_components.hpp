#pragma once

// Planet state grouped by world component (docs/architecture/world-model.md).
// Data only: planet's accessors and the simulation modules hold the behaviour.
// Grids are worldgrid<T> sized to the world; climate state is climateworld.

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "core/classes.hpp"
#include "simulations/climate/climate_biomes.hpp"
#include "simulations/geology/mineral_systems/mineral_deposit.hpp"
#include "simulations/society/social_generation.hpp"
#include "core/world/world_grid.hpp"

#define ARRAYWIDTH 2048
#define ARRAYHEIGHT 1024

#define MAXCRATERS 100000

constexpr int CLIMATEMAPSEACOLOURCOUNT = 3;
constexpr int CLIMATEMAPCOLOURCOUNT = 32;
constexpr int RIVERMAPCOLOURCOUNT = 10;
constexpr int RIVERMAPFEATURECOUNT = 7;
constexpr int MAPGRADIENTTYPECOUNT = 4;
constexpr int MAPGRADIENTMAXSTOPS = 12;

enum ClimateMapSeaColourSlot
{
    climateopensea = 0,
    climateseasonalseaice = 1,
    climatepermanentseaice = 2
};

enum RiverMapColourSlot
{
    rivermapbackground = 0,
    rivermaplowflow = 1,
    rivermaphighflow = 2,
    rivermaplake = 3,
    rivermapsaltpan = 4,
    rivermapwetlands = 5,
    rivermapmud = 6,
    rivermapsand = 7,
    rivermapshingle = 8,
    rivermapvolcano = 9
};

enum RiverMapFeatureSlot
{
    rivermapshowlakes = 0,
    rivermapshowsaltpans = 1,
    rivermapshowwetlands = 2,
    rivermapshowmud = 3,
    rivermapshowsand = 4,
    rivermapshowshingle = 5,
    rivermapshowvolcanoes = 6
};

enum MapGradientType
{
    mapgradientelevation = 0,
    mapgradienttemperature = 1,
    mapgradientprecipitation = 2,
    mapgradientriverflow = 3
};

enum class GeologicRegime : std::uint8_t
{
    stable = 0,
    convergent_arc = 1,
    continent_collision = 2,
    divergent_rift = 3,
    transform = 4,
    passive_margin = 5,
    mid_ocean_ridge = 6,
    trench_adjacent = 7
};

enum class CrustClass : std::uint8_t
{
    none = 0,
    oceanic = 1,
    transitional = 2,
    continental = 3
};

enum class BoundaryType : std::uint8_t
{
    none = 0,
    convergent = 1,
    divergent = 2,
    transform = 3,
    passive_margin = 4
};

enum class DeformingRegionType : std::uint8_t
{
    none = 0,
    continental_rift = 1,
    diffuse_collision = 2
};

enum class BasinClass : std::uint8_t
{
    none = 0,
    exorheic = 1,
    endorheic = 2,
    coastal = 3
};

struct TectonicBoundarySegment
{
    int id = 0;
    int leftPlateId = 0;
    int rightPlateId = 0;
    int cellCount = 0;
    int persistenceSteps = 0;
    float centroidX = 0.0f;
    float centroidY = 0.0f;
    float lengthCells = 0.0f;
    float averageNormalMotion = 0.0f;
    float averageShearMotion = 0.0f;
    int averageConvergenceScore = 0;
    int averageDivergenceScore = 0;
    int averageShearScore = 0;
    BoundaryType boundaryType = BoundaryType::none;
    GeologicRegime geologicRegime = GeologicRegime::stable;
    double ageMyr = 0.0;
};

struct TectonicDeformingRegion
{
    int id = 0;
    int boundarySegmentId = 0;
    int primaryPlateId = 0;
    int secondaryPlateId = 0;
    int cellCount = 0;
    int persistenceSteps = 0;
    float centroidX = 0.0f;
    float centroidY = 0.0f;
    float averageDeformationRate = 0.0f;
    float averageInterpolatedVelocityX = 0.0f;
    float averageInterpolatedVelocityY = 0.0f;
    float averageNormalMotion = 0.0f;
    float averageShearMotion = 0.0f;
    DeformingRegionType type = DeformingRegionType::none;
    double ageMyr = 0.0;
};

using namespace std;

// Grid geometry, orbit and rotation: constants of one planet.
struct planetbodystate
{
    int itssize; // Size of planet. 1=tiny; 2=small; 3=Earthlike.
    float itsradius_m = 0.0f; // Planet radius in metres (world.body.radius_km).
    int itswidth;   // width of global map
    int itsheight;  // height of global map
    long itsseed;    // seed number of this world
    bool itsrotation;    // 1 (true) like the Earth, 0 (false) the other way
    float itstilt; // Axial tilt (affects seasonal differences). Earthlike = 22.5
    float itseccentricity; // How elliptical the orbit is. Earthlike = 0.0167.
    short itsperihelion; // Closest point to the sun. Earthlike = 0 (Jan)
    float itsgravity; // Strength of gravity on the surface. Earthlike = 1.0
    float itslunar; // Strength of lunar attraction. Earthlike = 1.0
};

// Elevation, crust, landforms, plates, boundaries and mineral deposits.
struct geospherestate
{
    worldgrid<bool> shelfmap; // continental shelves: shallow margins of land
    int itstectonictimeoriginstep;
    float itstectonictimemyr;
    float itstectonicdeltatimemyr;
    int itstectonicplatecount;
    int itstectonicsealevelm;
    float itsmountainreduce;    // Amount to reduce the heights by, to make the global map match the regional map a little better
    int itsmaxheight;   // Maximum elevation
    int itssealevel;    // Sea level
    int itslandtotal; // Total cells of land
    int itsseatotal; // Total cells of sea
    int itscraterno; // Number of craters
    worldgrid<float> roughnessmap;
    worldgrid<short> mountainridges;
    worldgrid<float> mountainheights;
    worldgrid<float> craterrims;
    worldgrid<float> cratercentres;
    worldgrid<float> mapnom;   // metres above the datum
    worldgrid<std::uint8_t> geologicregimemap;
    worldgrid<std::uint8_t> tectonicconvergencemap;
    worldgrid<std::uint8_t> tectonicdivergencemap;
    worldgrid<std::uint8_t> tectonicshearmap;
    worldgrid<float> tectoniccrustagemyrmap;
    worldgrid<float> tectoniccrustthicknessmap;
    worldgrid<std::uint8_t> tectoniccrustclassmap;
    worldgrid<float> tectonicuplifttendencymap;
    worldgrid<float> tectonicsubsidencetendencymap;
    worldgrid<float> tectonicaccumulatedstrainmap;
    worldgrid<std::uint8_t> tectonicboundarytypemap;
    worldgrid<unsigned short> tectonicboundarydistancemap;
    worldgrid<int> tectonicboundarysegmentidmap;
    worldgrid<int> tectonicnearestboundaryidmap;
    worldgrid<float> tectonicboundaryhistorymap;
    worldgrid<int> tectonicdeformingregionidmap;
    worldgrid<std::uint8_t> tectonicdeformingregiontypemap;
    worldgrid<float> tectonicdeformationratemap;
    worldgrid<float> tectonicdeformationvelocityxmap;
    worldgrid<float> tectonicdeformationvelocityymap;
    // The engine's own fields of its last update, for the tectonic reports; no later process reads them.
    worldgrid<float> tectonicheightmap;
    worldgrid<int> tectonicplateidmap;
    worldgrid<int> tectoniccrustagestepsmap;
    worldgrid<float> tectoniccontinentalcrustmap;
    worldgrid<float> tectonicoceaniccrustmap;
    worldgrid<float> tectonicsedimentmap;
    worldgrid<float> tectonicriftdamagemap;
    worldgrid<float> tectonicplumeinfluencemap;
    std::string tectonicsnapshotrecord;   // JSON: time, scenario, diagnostics, statistics, plates
    worldgrid<std::uint8_t> erosionpotentialmap;
    worldgrid<std::uint8_t> depositionpotentialmap;
    worldgrid<std::uint8_t> floodplainfertilitymap;
    worldgrid<float> erodibilitymap, erosionratemap, depositionratemap, denudationmap, aggradationmap;
    worldgrid<std::uint8_t> progradedmap;
    float itserosionduration = 0.0f, itserodedvolume = 0.0f, itslanddepositedvolume = 0.0f, itsmarinedepositedvolume = 0.0f;
    worldgrid<std::uint8_t> metalorereservemap;
    worldgrid<std::uint8_t> placerreservemap;
    worldgrid<std::uint8_t> evaporitereservemap;
    worldgrid<std::uint8_t> volcanicreservemap;
    worldgrid<float> extraelevmap;
    worldgrid<bool> islandmap;
    worldgrid<bool> mountainislandmap;
    worldgrid<short> noisemap;
    worldgrid<short> oceanridgemap;
    worldgrid<float> oceanridgeheightmap;
    worldgrid<short> oceanriftmap;
    worldgrid<short> oceanridgeoffsetmap;
    worldgrid<short> oceanridgeanglemap;
    worldgrid<float> volcanomap;
    worldgrid<bool> stratomap;
    fourshorts cratercentreslist[MAXCRATERS]; // This duplicates the information in the cratercentres array, but it tells us the order in which to draw the craters, and their radius, which is also important.
    std::vector<TectonicBoundarySegment> tectonicboundarysegmentlist;
    std::vector<TectonicDeformingRegion> tectonicdeformingregionlist;
    std::vector<mineralsystems::Deposit> mineraldepositlist;
    mineralsystems::MineralBudget mineralbudgetstate;
};

// Rivers, lakes, deltas, tides, coastal features and drainage basins.
struct hydrospherestate
{
    float itsriverfactor;   // Divide river flow by this to get it in cubic metres/second
    int itsriverlandreduce; // Rivers will carve channels at the base elevation of the tile, minus this
    int itsestuarylimit;    // Rivers this size or above might have estuaries
    int itsmaxriverflow;    // Largest river flow
    worldgrid<short> rivermapdir;
    worldgrid<int> rivermapjan;
    worldgrid<int> rivermapjul;
    worldgrid<float> lakemap;
    worldgrid<std::uint8_t> flowdirectionmap;
    worldgrid<float> drainageareamap, dischargemap, filledelevationmap;
    worldgrid<std::int32_t> depressionidmap;
    worldgrid<float> runoffmap, sedimentfluxmap;
    worldgrid<std::uint8_t> importedlakemap;
    worldgrid<short> tidalmap;
    worldgrid<float> riftlakemapsurface;
    worldgrid<float> riftlakemapbed;
    worldgrid<bool> lakestartmap;
    worldgrid<short> specials;
    worldgrid<std::uint8_t> basinclassmap;
    worldgrid<short> deltamapdir;
    worldgrid<int> deltamapjan;
    worldgrid<int> deltamapjul;
};

// Glacial thresholds; sea ice is climate state (climateworld).
struct cryospherestate
{
    int itsglacialtemp;     // Areas with this average temperature or lower will have glacial terrain
    int itsglaciertemp;     // Rivers in areas with this maximum temperature or lower will be glaciers
};

// Planet parameters of the climate; climate fields live in climateworld.
struct atmosphereparameters
{
    float itstempdecrease; // How much temperature decreases with elevation. Earthlike = 6.5/km
    int itsnorthpolaradjust; // Amount to adjust north pole temperature. Earthlike = +3
    int itssouthpolaradjust; // Amount to adjust south pole temperature. Earthlike = -3
    int itsaveragetemp; // Average global temperature. Earthlike = 14
    float itswaterpickup; // How much water to pick up over oceans. Default = 1.0
    int itsclimateno;   // Total number of climate types
};

// Biological resources beyond the climate biomes (climateworld).
struct biospherestate
{
    worldgrid<std::uint8_t> fisheryreservemap;
};

// Settlements, polities, routes, ownership and history.
struct anthropospherestate
{
    worldgrid<std::uint8_t> settlement_suitability;
    worldgrid<std::uint8_t> social_infrastructure;
    worldgrid<std::uint8_t> agricultural_capacity;
    worldgrid<std::uint8_t> route_traffic;
    worldgrid<std::uint8_t> river_access;
    worldgrid<std::uint8_t> harbor_score;
    worldgrid<int> owner_settlement_id;
    worldgrid<int> owner_polity_id;
    std::vector<Settlement> settlementlist;
    std::vector<Polity> politylist;
    std::vector<RouteEdge> routeedgelist;
    std::vector<HistoryEvent> historyeventlist;
};

// Map colours, shading and display options; never simulation state.
struct presentationsettings
{
    float itslandshading;
    float itslakeshading;
    float itsseashading; // Amount of shading to do on these areas.
    short itsshadingdir = 4; // Direction of light source for shading.
    short itssnowchange = 1; // Kind of snow border.
    short itsseaiceappearance = 1; // How much sea ice to display.
    bool itscolourcliffs = 1; // If this is 1, high grass colours will only be used on steep slopes.
    bool itsshowmapoutline = 1;
    float itslandmarbling;
    float itslakemarbling;
    float itsseamarbling; // Amount of marbling to do on these areas on the regional map.
    int itsminriverflowglobal; // Rivers greater than this will be shown on the global map.
    int itsminriverflowregional; // Rivers greater than this will be shown on the regional map.
    bool itsmangroves;
    int itsseaice1;
    int itsseaice2;
    int itsseaice3;
    int itsocean1;
    int itsocean2;
    int itsocean3;
    int itsdeepocean1;
    int itsdeepocean2;
    int itsdeepocean3;
    int itsbase1;
    int itsbase2;
    int itsbase3;
    int itsbasetemp1;
    int itsbasetemp2;
    int itsbasetemp3;
    int itshighbase1;
    int itshighbase2;
    int itshighbase3;
    int itsdesert1;
    int itsdesert2;
    int itsdesert3;
    int itshighdesert1;
    int itshighdesert2;
    int itshighdesert3;
    int itscolddesert1;
    int itscolddesert2;
    int itscolddesert3;
    int itsgrass1;
    int itsgrass2;
    int itsgrass3;
    int itscold1;
    int itscold2;
    int itscold3;
    int itstundra1;
    int itstundra2;
    int itstundra3;
    int itseqtundra1;
    int itseqtundra2;
    int itseqtundra3;
    int itssaltpan1;
    int itssaltpan2;
    int itssaltpan3;
    int itserg1;
    int itserg2;
    int itserg3;
    int itswetlands1;
    int itswetlands2;
    int itswetlands3;
    int itslake1;
    int itslake2;
    int itslake3;
    int itsriver1;
    int itsriver2;
    int itsriver3;
    int itsglacier1;
    int itsglacier2;
    int itsglacier3;
    int itssand1;
    int itssand2;
    int itssand3;
    int itsmud1;
    int itsmud2;
    int itsmud3;
    int itsshingle1;
    int itsshingle2;
    int itsshingle3;
    int itsmangrove1;
    int itsmangrove2;
    int itsmangrove3;
    int itshighlight1;
    int itshighlight2;
    int itshighlight3;
    int itsoutline1;
    int itsoutline2;
    int itsoutline3;
    int itselevationlow1;
    int itselevationlow2;
    int itselevationlow3;
    int itselevationhigh1;
    int itselevationhigh2;
    int itselevationhigh3;
    int itstemperaturecold1;
    int itstemperaturecold2;
    int itstemperaturecold3;
    int itstemperaturetemperate1;
    int itstemperaturetemperate2;
    int itstemperaturetemperate3;
    int itstemperaturehot1;
    int itstemperaturehot2;
    int itstemperaturehot3;
    int itsprecipitationdry1;
    int itsprecipitationdry2;
    int itsprecipitationdry3;
    int itsprecipitationwet1;
    int itsprecipitationwet2;
    int itsprecipitationwet3;
    std::array<std::array<int, 3>, CLIMATEMAPSEACOLOURCOUNT> itsclimatemapseacolours;
    std::array<std::array<int, 3>, CLIMATEMAPCOLOURCOUNT> itsclimatemapcolours;
    std::array<std::array<int, 3>, BIOMEMAPCOLOURCOUNT> itsbiomemapcolours;
    std::array<std::array<int, 3>, RIVERMAPCOLOURCOUNT> itsrivermapcolours;
    std::array<bool, RIVERMAPFEATURECOUNT> itsshowrivermapfeatures;
    std::array<int, MAPGRADIENTTYPECOUNT> itsmapgradientstopcounts{};
    std::array<bool, MAPGRADIENTTYPECOUNT> itsmapgradientdiscrete{};
    std::array<std::array<int, MAPGRADIENTMAXSTOPS>, MAPGRADIENTTYPECOUNT> itsmapgradientpositions{};
    std::array<std::array<std::array<int, 3>, MAPGRADIENTMAXSTOPS>, MAPGRADIENTTYPECOUNT> itsmapgradientcolours{};
    worldgrid<bool> noshademap;
    worldgrid<int> testmap;
};
