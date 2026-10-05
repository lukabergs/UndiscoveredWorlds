#pragma once

// Global and regional map views, cached map images and drawing entry points.
// Desktop only: these use SFML images and textures. Split from functions.hpp
// so that simulation and pipeline code compile without SFML.

#include <array>

#include <SFML/Graphics.hpp>

#include "world/functions.hpp"

constexpr int GLOBALMAPTYPES = 22;
constexpr int DISPLAYMAPSIZEX = 1024;
constexpr int DISPLAYMAPSIZEY = 512;

enum mapviewenum
{
    elevation,
    temperature,
    precipitation,
    climate,
    biomes,
    rivers,
    relief,
    geology,
    tectonic_uplift,
    tectonic_boundaries,
    basins,
    erosion,
    deposition,
    fertility,
    resources,
    suitability,
    settlements_map,
    population,
    infrastructure_map,
    polities_map,
    trade_map,
    recent_conflict_map
};
enum mapdatakindenum { indexedmapdata, gradientmapdata };
enum mapindexedstyleenum { noindexedstyle, reliefindexedstyle, climateindexedstyle, biomeindexedstyle };
enum mapgradientstyleenum { nogradientstyle, standardgradientstyle, riversgradientstyle };

struct AppearanceSettings;
struct maplayer;

void drawglobalelevationmapimage(planet& world, maplayer& layer);
void drawglobaltemperaturemapimage(planet& world, maplayer& layer);
void drawglobalprecipitationmapimage(planet& world, maplayer& layer);
void drawglobalclimatemapimage(planet& world, maplayer& layer);
void drawglobalbiomemapimage(planet& world, maplayer& layer);
void drawglobalriversmapimage(planet& world, maplayer& layer);
void drawglobalreliefmapimage(planet& world, maplayer& layer);
void drawglobalgeologymapimage(planet& world, maplayer& layer);
void drawglobaltectonicupliftmapimage(planet& world, maplayer& layer);
void drawglobaltectonicboundariesmapimage(planet& world, maplayer& layer);
void drawglobalbasinsmapimage(planet& world, maplayer& layer);
void drawglobalerosionmapimage(planet& world, maplayer& layer);
void drawglobaldepositionmapimage(planet& world, maplayer& layer);
void drawglobalfertilitymapimage(planet& world, maplayer& layer);
void drawglobalresourcesmapimage(planet& world, maplayer& layer);
void drawglobalsuitabilitymapimage(planet& world, maplayer& layer);
void drawglobalsettlementsmapimage(planet& world, maplayer& layer);
void drawglobalpopulationmapimage(planet& world, maplayer& layer);
void drawglobalinfrastructuremapimage(planet& world, maplayer& layer);
void drawglobalpolitiesmapimage(planet& world, maplayer& layer);
void drawglobaltrademapimage(planet& world, maplayer& layer);
void drawglobalrecentconflictmapimage(planet& world, maplayer& layer);
void drawregionalelevationmapimage(planet& world, region& region, maplayer& layer);
void drawregionaltemperaturemapimage(planet& world, region& region, maplayer& layer);
void drawregionalprecipitationmapimage(planet& world, region& region, maplayer& layer);
void drawregionalclimatemapimage(planet& world, region& region, maplayer& layer);
void drawregionalbiomemapimage(planet& world, region& region, maplayer& layer);
void drawregionalriversmapimage(planet& world, region& region, maplayer& layer);
void drawregionalreliefmapimage(planet& world, region& region, maplayer& layer);
void drawregionalgeologymapimage(planet& world, region& region, maplayer& layer);
void drawregionaltectonicupliftmapimage(planet& world, region& region, maplayer& layer);
void drawregionaltectonicboundariesmapimage(planet& world, region& region, maplayer& layer);
void drawregionalbasinsmapimage(planet& world, region& region, maplayer& layer);
void drawregionalerosionmapimage(planet& world, region& region, maplayer& layer);
void drawregionaldepositionmapimage(planet& world, region& region, maplayer& layer);
void drawregionalfertilitymapimage(planet& world, region& region, maplayer& layer);
void drawregionalresourcesmapimage(planet& world, region& region, maplayer& layer);
void drawregionalsuitabilitymapimage(planet& world, region& region, maplayer& layer);
void drawregionalsettlementsmapimage(planet& world, region& region, maplayer& layer);
void drawregionalpopulationmapimage(planet& world, region& region, maplayer& layer);
void drawregionalinfrastructuremapimage(planet& world, region& region, maplayer& layer);
void drawregionalpolitiesmapimage(planet& world, region& region, maplayer& layer);
void drawregionaltrademapimage(planet& world, region& region, maplayer& layer);
void drawregionalrecentconflictmapimage(planet& world, region& region, maplayer& layer);
void drawgradientmapappearance(const struct mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, std::array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int colouralign, int otheralign);
void drawindexedmapappearance(const struct mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, std::array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int colouralign, int otheralign);
void drawstaticmapappearance(const struct mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, std::array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int colouralign, int otheralign);

struct mapviewdefinition
{
    mapviewenum view;
    const char* label;
    const char* exportstem;
    mapdatakindenum datakind;
    mapindexedstyleenum indexedstyle;
    mapgradientstyleenum gradientstyle;
    int gradientindex;
    const char* description;
    void (*drawglobal)(planet& world, maplayer& layer);
    void (*drawregional)(planet& world, region& region, maplayer& layer);
    void (*drawappearance)(const mapviewdefinition& definition, planet& world, AppearanceSettings& appearance, std::array<int, MAPGRADIENTTYPECOUNT>& selectedgradientstops, int colouralign, int otheralign);
};

constexpr std::array<mapviewenum, GLOBALMAPTYPES> allmapviews = { elevation, temperature, precipitation, climate, biomes, rivers, relief, geology, tectonic_uplift, tectonic_boundaries, basins, erosion, deposition, fertility, resources, suitability, settlements_map, population, infrastructure_map, polities_map, trade_map, recent_conflict_map };
constexpr std::array<mapviewdefinition, GLOBALMAPTYPES> allmapviewdefinitions =
{ {
    { relief, "Relief", "Relief", indexedmapdata, reliefindexedstyle, nogradientstyle, -1, "", drawglobalreliefmapimage, drawregionalreliefmapimage, drawindexedmapappearance },
    { elevation, "Elevation", "Elevation", gradientmapdata, noindexedstyle, standardgradientstyle, mapgradientelevation, "0 m is sea level. Negative values are below sea level.", drawglobalelevationmapimage, drawregionalelevationmapimage, drawgradientmapappearance },
    { temperature, "Temperature", "Temperature", gradientmapdata, noindexedstyle, standardgradientstyle, mapgradienttemperature, "", drawglobaltemperaturemapimage, drawregionaltemperaturemapimage, drawgradientmapappearance },
    { precipitation, "Precipitation", "Precipitation", gradientmapdata, noindexedstyle, standardgradientstyle, mapgradientprecipitation, "", drawglobalprecipitationmapimage, drawregionalprecipitationmapimage, drawgradientmapappearance },
    { climate, "Climate", "Climate", indexedmapdata, climateindexedstyle, nogradientstyle, -1, "", drawglobalclimatemapimage, drawregionalclimatemapimage, drawindexedmapappearance },
    { biomes, "Biomes", "Biomes", indexedmapdata, biomeindexedstyle, nogradientstyle, -1, "Sea uses the Climate tab sea palette.", drawglobalbiomemapimage, drawregionalbiomemapimage, drawindexedmapappearance },
    { rivers, "Rivers", "Rivers", gradientmapdata, noindexedstyle, riversgradientstyle, mapgradientriverflow, "", drawglobalriversmapimage, drawregionalriversmapimage, drawgradientmapappearance },
    { geology, "Geology", "Geology", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Derived geologic regimes from retained tectonic signals.", drawglobalgeologymapimage, drawregionalgeologymapimage, drawstaticmapappearance },
    { tectonic_uplift, "Tectonic Uplift", "TectonicUplift", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Native uplift/subsidence tendency with accumulated strain memory.", drawglobaltectonicupliftmapimage, drawregionaltectonicupliftmapimage, drawstaticmapappearance },
    { tectonic_boundaries, "Tectonic Boundaries", "TectonicBoundaries", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Native boundary class, proximity, lifecycle history, and deforming-region overlay.", drawglobaltectonicboundariesmapimage, drawregionaltectonicboundariesmapimage, drawstaticmapappearance },
    { basins, "Basins", "Basins", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Drainage basin classes derived from final hydrology.", drawglobalbasinsmapimage, drawregionalbasinsmapimage, drawstaticmapappearance },
    { erosion, "Erosion", "Erosion", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Relative erosion potential, normalized to 0..100.", drawglobalerosionmapimage, drawregionalerosionmapimage, drawstaticmapappearance },
    { deposition, "Deposition", "Deposition", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Relative deposition potential, normalized to 0..100.", drawglobaldepositionmapimage, drawregionaldepositionmapimage, drawstaticmapappearance },
    { fertility, "Fertility", "Fertility", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Floodplain and lowland fertility proxy, normalized to 0..100.", drawglobalfertilitymapimage, drawregionalfertilitymapimage, drawstaticmapappearance },
    { resources, "Resources", "Resources", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Dominant reserve-style physical resource potential. Gold tints mark surface-accessible lode and placer gold.", drawglobalresourcesmapimage, drawregionalresourcesmapimage, drawstaticmapappearance },
    { suitability, "Suitability", "Suitability", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Settlement suitability score (0..100).", drawglobalsuitabilitymapimage, drawregionalsuitabilitymapimage, drawstaticmapappearance },
    { settlements_map, "Settlements", "Settlements", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Settlement locations and urban size.", drawglobalsettlementsmapimage, drawregionalsettlementsmapimage, drawstaticmapappearance },
    { population, "Population", "Population", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Combined rural and urban population intensity.", drawglobalpopulationmapimage, drawregionalpopulationmapimage, drawstaticmapappearance },
    { infrastructure_map, "Infrastructure", "Infrastructure", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Infrastructure diffusion and hinterland reach.", drawglobalinfrastructuremapimage, drawregionalinfrastructuremapimage, drawstaticmapappearance },
    { polities_map, "Polities", "Polities", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Current polity ownership.", drawglobalpolitiesmapimage, drawregionalpolitiesmapimage, drawstaticmapappearance },
    { trade_map, "Trade", "Trade", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Route traffic and trade intensity.", drawglobaltrademapimage, drawregionaltrademapimage, drawstaticmapappearance },
    { recent_conflict_map, "Recent Conflict", "RecentConflict", indexedmapdata, noindexedstyle, nogradientstyle, -1, "Border friction and recent conflict pressure.", drawglobalrecentconflictmapimage, drawregionalrecentconflictmapimage, drawstaticmapappearance },
} };

constexpr const mapviewdefinition& getmapviewdefinition(mapviewenum view)
{
    for (const mapviewdefinition& definition : allmapviewdefinitions)
    {
        if (definition.view == view)
            return definition;
    }

    return allmapviewdefinitions[0];
}

struct maplayer
{
    sf::Image image;
    sf::Image displayimage;
    bool created = false;
};

struct mapcache
{
    std::array<maplayer, GLOBALMAPTYPES> layers;
};

constexpr int mapviewindex(mapviewenum mapview)
{
    return static_cast<int>(mapview);
}

void resetmapcache(mapcache& maps);
maplayer& getmaplayer(mapcache& maps, mapviewenum mapview);
const maplayer& getmaplayer(const mapcache& maps, mapviewenum mapview);
sf::Image& getmapimage(mapcache& maps, mapviewenum mapview);
const sf::Image& getmapimage(const mapcache& maps, mapviewenum mapview);
sf::Image& getdisplaymapimage(mapcache& maps, mapviewenum mapview);
const sf::Image& getdisplaymapimage(const mapcache& maps, mapviewenum mapview);
void updateTextureFromImage(sf::Texture& texture, const sf::Image& image);
void adjustforsize(planet& world, sf::Vector2i& globaltexturesize, mapcache& globalmaps, sf::Image& highlightimage, int highlightsize, sf::Image& minihighlightimage, int& minihighlightsize);
void drawhighlightobjects(planet& world, sf::Image& highlightimage, int highlightsize, sf::Image& minihighlightimage, int& minihighlightsize);

bool standardbutton(const char* label);

void drawglobalmapimage(mapviewenum mapview, planet& world, mapcache& maps);
void drawallglobalmapimages(planet& world, mapcache& maps);
void applyglobalmapview(mapviewenum mapview, planet& world, mapcache& maps, sf::Texture& texture, sf::Sprite& sprite, sf::Sprite* minimap = nullptr);
void drawglobalelevationmapimage(planet& world, maplayer& layer);
void drawglobaltemperaturemapimage(planet& world, maplayer& layer);
void drawglobalprecipitationmapimage(planet& world, maplayer& layer);
void drawglobalclimatemapimage(planet& world, maplayer& layer);
void drawglobalbiomemapimage(planet& world, maplayer& layer);
void drawglobalriversmapimage(planet& world, maplayer& layer);
void drawglobalreliefmapimage(planet& world, maplayer& layer);
void drawglobalgeologymapimage(planet& world, maplayer& layer);
void drawglobaltectonicupliftmapimage(planet& world, maplayer& layer);
void drawglobaltectonicboundariesmapimage(planet& world, maplayer& layer);
void drawglobalbasinsmapimage(planet& world, maplayer& layer);
void drawglobalerosionmapimage(planet& world, maplayer& layer);
void drawglobaldepositionmapimage(planet& world, maplayer& layer);
void drawglobalfertilitymapimage(planet& world, maplayer& layer);
void drawglobalresourcesmapimage(planet& world, maplayer& layer);
void drawglobalsuitabilitymapimage(planet& world, maplayer& layer);
void drawglobalsettlementsmapimage(planet& world, maplayer& layer);
void drawglobalpopulationmapimage(planet& world, maplayer& layer);
void drawglobalinfrastructuremapimage(planet& world, maplayer& layer);
void drawglobalpolitiesmapimage(planet& world, maplayer& layer);
void drawglobaltrademapimage(planet& world, maplayer& layer);
void drawglobalrecentconflictmapimage(planet& world, maplayer& layer);
sf::Color getclimatecolours(const planet& world, short climate);
void drawregionalmapimage(mapviewenum mapview, planet& world, region& region, mapcache& maps);
void drawallregionalmapimages(planet& world, region& region, mapcache& maps);
void applyregionalmapview(mapviewenum mapview, planet& world, region& region, mapcache& maps, sf::Texture& texture, sf::Sprite& sprite);
void drawregionalelevationmapimage(planet& world, region& region, maplayer& layer);
void drawregionaltemperaturemapimage(planet& world, region& region, maplayer& layer);
void drawregionalprecipitationmapimage(planet& world, region& region, maplayer& layer);
void drawregionalclimatemapimage(planet& world, region& region, maplayer& layer);
void drawregionalbiomemapimage(planet& world, region& region, maplayer& layer);
void drawregionalriversmapimage(planet& world, region& region, maplayer& layer);
void drawregionalreliefmapimage(planet& world, region& region, maplayer& layer);
void drawregionalgeologymapimage(planet& world, region& region, maplayer& layer);
void drawregionaltectonicupliftmapimage(planet& world, region& region, maplayer& layer);
void drawregionaltectonicboundariesmapimage(planet& world, region& region, maplayer& layer);
void drawregionalbasinsmapimage(planet& world, region& region, maplayer& layer);
void drawregionalerosionmapimage(planet& world, region& region, maplayer& layer);
void drawregionaldepositionmapimage(planet& world, region& region, maplayer& layer);
void drawregionalfertilitymapimage(planet& world, region& region, maplayer& layer);
void drawregionalresourcesmapimage(planet& world, region& region, maplayer& layer);
void drawregionalsuitabilitymapimage(planet& world, region& region, maplayer& layer);
void drawregionalsettlementsmapimage(planet& world, region& region, maplayer& layer);
void drawregionalpopulationmapimage(planet& world, region& region, maplayer& layer);
void drawregionalinfrastructuremapimage(planet& world, region& region, maplayer& layer);
void drawregionalpolitiesmapimage(planet& world, region& region, maplayer& layer);
void drawregionaltrademapimage(planet& world, region& region, maplayer& layer);
void drawregionalrecentconflictmapimage(planet& world, region& region, maplayer& layer);
