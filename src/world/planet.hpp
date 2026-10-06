//
//  planet.hpp
//  Undiscovered Worlds
//
//  Created by Jonathan Hill on 22/07/2019.
//
//  Please see functions.hpp for notes.

#ifndef planet_hpp
#define planet_hpp

#include <stdio.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "core/classes.hpp"
#include "core/climate_seasons.hpp"
#include "simulations/climate/climate_context.hpp"
#include "simulations/geology/mineral_systems/mineral_deposit.hpp"
#include "simulations/society/social_generation.hpp"
#include "core/world/field_access.hpp"
#include "simulations/climate/climate_biomes.hpp"
#include "simulations/climate/climate_world.hpp"

#include "world/planet_components.hpp"

class planet
{
    // Catalogue access for run bundles and access audits (planet_fields.cpp).
    friend struct planetfieldaccess;

public:

    planet(); // constructor
    ~planet();  // destructor

    climatecontext::Context& climatesimulation() { return climate_.climatesimulation(); }
    const climatecontext::Context& climatesimulation() const { return climate_.climatesimulation(); }
    // Climate module state with this world's terrain bound as its inputs.
    climateworld& climate();
    const climateworld& climate() const;

    // accessor functions

    int saveversion() const; // save format
    void setsaveversion(int amount);

    int settingssaveversion() const; // save format for settings
    void setsettingssaveversion(int amount);

    int tectonictimeoriginstep() const;
    void settectonictimeoriginstep(int amount);

    float tectonictimemyr() const;
    void settectonictimemyr(float amount);

    float tectonicdeltatimemyr() const;
    void settectonicdeltatimemyr(float amount);


    int tectonicplatecount() const;
    void settectonicplatecount(int amount);

    int tectonicsealevelm() const;
    void settectonicsealevelm(int amount);

    long seed() const;   // seed
    void setseed(long amount);

    int size() const;   // size
    void setsize(int amount);

    int width() const;   // global width
    void setwidth(int amount);

    int height() const;  // global height
    void setheight(int amount);

    bool rotation() const;   // global rotation
    void setrotation(bool amount);

    float tilt() const; // axial tilt
    void settilt(float amount);

    float eccentricity() const; // orbital eccentricity
    void seteccentricity(float amount);

    int perihelion() const; // perihelion
    void setperihelion (int amount);

    float gravity() const; // surface gravity
    void setgravity(float amount);

    float lunar() const; // strength of lunar gravity
    void setlunar(float amount);

    float tempdecrease() const; // amount temperatures decrease by per 1000 metres
    void settempdecrease(float amount);

    int northpolaradjust() const; // adjustment to north polar temperature
    void setnorthpolaradjust(int amount);

    int southpolaradjust() const; // adjustment to south polar temperature
    void setsouthpolaradjust(int amount);

    int averagetemp() const; // average world temperature
    void setaveragetemp(int amount);

    int northpolartemp() const; // north polar temperature
    void setnorthpolartemp(int amount);

    int southpolartemp() const; // south polar temperature
    void setsouthpolartemp(int amount);

    int eqtemp() const; // equatorial temperature
    void seteqtemp(int amount);

    float waterpickup() const; // amount of water to pick up over oceans
    void setwaterpickup(float amount);

    float riverfactor() const;   // global factor for calculating flow in cubic metres/second
    void setriverfactor(float amount);

    int riverlandreduce() const; // global factor for carving river valleys
    void setriverlandreduce(int amount);

    int estuarylimit() const; // global minimum river flow for estuaries
    void setestuarylimit(int amount);

    int glacialtemp() const; // global maximum temperature for glacially carved terrain
    void setglacialtemp(int amount);

    int glaciertemp() const; // global maximum temperature for rivers to become glaciers
    void setglaciertemp(int amount);

    float mountainreduce() const;    // global factor to reduce mountain heights by a little
    void setmountainreduce(float amount);

    int climatenumber() const;   // total number of possible climates
    void setclimatenumber(int amount);

    int maxelevation() const;    // global maximum elevation
    void setmaxelevation(int amount);

    int sealevel() const;    // global sea level
    void setsealevel(int amount);

    float landshading() const;
    void setlandshading(float amount);

    float lakeshading() const;
    void setlakeshading(float amount);

    float seashading() const;
    void setseashading(float amount);

    int snowchange() const;
    void setsnowchange(int amount); // 1=abrupt; 2=speckled; 3=gradual

    int seaiceappearance() const;
    void setseaiceappearance(int amount);

    bool colourcliffs() const;
    void setcolourcliffs(bool amount);
    bool showmapoutline() const;
    void setshowmapoutline(bool amount);

    int shadingdir() const;
    void setshadingdir(int amount);

    float landmarbling() const;
    void setlandmarbling(float amount);

    float lakemarbling() const;
    void setlakemarbling(float amount);

    float seamarbling() const;
    void setseamarbling(float amount);

    int minriverflowglobal() const;
    void setminriverflowglobal(int amount);

    int minriverflowregional() const;
    void setminriverflowregional(int amount);

    bool showmangroves() const;
    void setshowmangroves(bool amount);

    int landtotal() const;
    void setlandtotal(int amount);

    int seatotal() const;
    void setseatotal(int amount);

    int craterno() const;
    void setcraterno(int amount);

    int maxriverflow() const;    // greatest river flow

    // These are the accessor functions for the map colours.

    int seaice1() const { return presentation_.itsseaice1; };
    void setseaice1(int amount) { presentation_.itsseaice1 = amount; };

    int seaice2() const { return presentation_.itsseaice2; };
    void setseaice2(int amount) { presentation_.itsseaice2 = amount; };

    int seaice3() const { return presentation_.itsseaice3; };
    void setseaice3(int amount) { presentation_.itsseaice3 = amount; };

    int ocean1() const { return presentation_.itsocean1; };
    void setocean1(int amount) { presentation_.itsocean1 = amount; };

    int ocean2() const { return presentation_.itsocean2; };
    void setocean2(int amount) { presentation_.itsocean2 = amount; };

    int ocean3() const { return presentation_.itsocean3; };
    void setocean3(int amount) { presentation_.itsocean3 = amount; };

    int deepocean1() const { return presentation_.itsdeepocean1; };
    void setdeepocean1(int amount) { presentation_.itsdeepocean1 = amount; };

    int deepocean2() const { return presentation_.itsdeepocean2; };
    void setdeepocean2(int amount) { presentation_.itsdeepocean2 = amount; };

    int deepocean3() const { return presentation_.itsdeepocean3; };
    void setdeepocean3(int amount) { presentation_.itsdeepocean3 = amount; };

    int base1() const { return presentation_.itsbase1; };
    void setbase1(int amount) { presentation_.itsbase1 = amount; };

    int base2() const { return presentation_.itsbase2; };
    void setbase2(int amount) { presentation_.itsbase2 = amount; };

    int base3() const { return presentation_.itsbase3; };
    void setbase3(int amount) { presentation_.itsbase3 = amount; };

    int basetemp1() const { return presentation_.itsbasetemp1; };
    void setbasetemp1(int amount) { presentation_.itsbasetemp1 = amount; };

    int basetemp2() const { return presentation_.itsbasetemp2; };
    void setbasetemp2(int amount) { presentation_.itsbasetemp2 = amount; };

    int basetemp3() const { return presentation_.itsbasetemp3; };
    void setbasetemp3(int amount) { presentation_.itsbasetemp3 = amount; };

    int highbase1() const { return presentation_.itshighbase1; };
    void sethighbase1(int amount) { presentation_.itshighbase1 = amount; };

    int highbase2() const { return presentation_.itshighbase2; };
    void sethighbase2(int amount) { presentation_.itshighbase2 = amount; };

    int highbase3() const { return presentation_.itshighbase3; };
    void sethighbase3(int amount) { presentation_.itshighbase3 = amount; };

    int desert1() const { return presentation_.itsdesert1; };
    void setdesert1(int amount) { presentation_.itsdesert1 = amount; };

    int desert2() const { return presentation_.itsdesert2; };
    void setdesert2(int amount) { presentation_.itsdesert2 = amount; };

    int desert3() const { return presentation_.itsdesert3; };
    void setdesert3(int amount) { presentation_.itsdesert3 = amount; };

    int highdesert1() const { return presentation_.itshighdesert1; };
    void sethighdesert1(int amount) { presentation_.itshighdesert1 = amount; };

    int highdesert2() const { return presentation_.itshighdesert2; };
    void sethighdesert2(int amount) { presentation_.itshighdesert2 = amount; };

    int highdesert3() const { return presentation_.itshighdesert3; };
    void sethighdesert3(int amount) { presentation_.itshighdesert3 = amount; };

    int colddesert1() const { return presentation_.itscolddesert1; };
    void setcolddesert1(int amount) { presentation_.itscolddesert1 = amount; };

    int colddesert2() const { return presentation_.itscolddesert2; };
    void setcolddesert2(int amount) { presentation_.itscolddesert2 = amount; };

    int colddesert3() const { return presentation_.itscolddesert3; };
    void setcolddesert3(int amount) { presentation_.itscolddesert3 = amount; };

    int grass1() const { return presentation_.itsgrass1; };
    void setgrass1(int amount) { presentation_.itsgrass1 = amount; };

    int grass2() const { return presentation_.itsgrass2; };
    void setgrass2(int amount) { presentation_.itsgrass2 = amount; };

    int grass3() const { return presentation_.itsgrass3; };
    void setgrass3(int amount) { presentation_.itsgrass3 = amount; };

    int cold1() const { return presentation_.itscold1; };
    void setcold1(int amount) { presentation_.itscold1 = amount; };

    int cold2() const { return presentation_.itscold2; };
    void setcold2(int amount) { presentation_.itscold2 = amount; };

    int cold3() const { return presentation_.itscold3; };
    void setcold3(int amount) { presentation_.itscold3 = amount; };

    int tundra1() const { return presentation_.itstundra1; };
    void settundra1(int amount) { presentation_.itstundra1 = amount; };

    int tundra2() const { return presentation_.itstundra2; };
    void settundra2(int amount) { presentation_.itstundra2 = amount; };

    int tundra3() const { return presentation_.itstundra3; };
    void settundra3(int amount) { presentation_.itstundra3 = amount; };

    int eqtundra1() const { return presentation_.itseqtundra1; };
    void seteqtundra1(int amount) { presentation_.itseqtundra1 = amount; };

    int eqtundra2() const { return presentation_.itseqtundra2; };
    void seteqtundra2(int amount) { presentation_.itseqtundra2 = amount; };

    int eqtundra3() const { return presentation_.itseqtundra3; };
    void seteqtundra3(int amount) { presentation_.itseqtundra3 = amount; };

    int saltpan1() const { return presentation_.itssaltpan1; };
    void setsaltpan1(int amount) { presentation_.itssaltpan1 = amount; };

    int saltpan2() const { return presentation_.itssaltpan2; };
    void setsaltpan2(int amount) { presentation_.itssaltpan2 = amount; };

    int saltpan3() const { return presentation_.itssaltpan3; };
    void setsaltpan3(int amount) { presentation_.itssaltpan3 = amount; };

    int erg1() const { return presentation_.itserg1; };
    void seterg1(int amount) { presentation_.itserg1 = amount; };

    int erg2() const { return presentation_.itserg2; };
    void seterg2(int amount) { presentation_.itserg2 = amount; };

    int erg3() const { return presentation_.itserg3; };
    void seterg3(int amount) { presentation_.itserg3 = amount; };

    int wetlands1() const { return presentation_.itswetlands1; };
    void setwetlands1(int amount) { presentation_.itswetlands1 = amount; };

    int wetlands2() const { return presentation_.itswetlands2; };
    void setwetlands2(int amount) { presentation_.itswetlands2 = amount; };

    int wetlands3() const { return presentation_.itswetlands3; };
    void setwetlands3(int amount) { presentation_.itswetlands3 = amount; };

    int lake1() const { return presentation_.itslake1; };
    void setlake1(int amount) { presentation_.itslake1 = amount; };

    int lake2() const { return presentation_.itslake2; };
    void setlake2(int amount) { presentation_.itslake2 = amount; };

    int lake3() const { return presentation_.itslake3; };
    void setlake3(int amount) { presentation_.itslake3 = amount; };

    int river1() const { return presentation_.itsriver1; };
    void setriver1(int amount) { presentation_.itsriver1 = amount; };

    int river2() const { return presentation_.itsriver2; };
    void setriver2(int amount) { presentation_.itsriver2 = amount; };

    int river3() const { return presentation_.itsriver3; };
    void setriver3(int amount) { presentation_.itsriver3 = amount; };

    int glacier1() const { return presentation_.itsglacier1; };
    void setglacier1(int amount) { presentation_.itsglacier1 = amount; };

    int glacier2() const { return presentation_.itsglacier2; };
    void setglacier2(int amount) { presentation_.itsglacier2 = amount; };

    int glacier3() const { return presentation_.itsglacier3; };
    void setglacier3(int amount) { presentation_.itsglacier3 = amount; };

    int sand1() const { return presentation_.itssand1; };
    void setsand1(int amount) { presentation_.itssand1 = amount; };

    int sand2() const { return presentation_.itssand2; };
    void setsand2(int amount) { presentation_.itssand2 = amount; };

    int sand3() const { return presentation_.itssand3; };
    void setsand3(int amount) { presentation_.itssand3 = amount; };

    int mud1() const { return presentation_.itsmud1; };
    void setmud1(int amount) { presentation_.itsmud1 = amount; };

    int mud2() const { return presentation_.itsmud2; };
    void setmud2(int amount) { presentation_.itsmud2 = amount; };

    int mud3() const { return presentation_.itsmud3; };
    void setmud3(int amount) { presentation_.itsmud3 = amount; };

    int shingle1() const { return presentation_.itsshingle1; };
    void setshingle1(int amount) { presentation_.itsshingle1 = amount; };

    int shingle2() const { return presentation_.itsshingle2; };
    void setshingle2(int amount) { presentation_.itsshingle2 = amount; };

    int shingle3() const { return presentation_.itsshingle3; };
    void setshingle3(int amount) { presentation_.itsshingle3 = amount; };

    int mangrove1() const { return presentation_.itsmangrove1; };
    void setmangrove1(int amount) { presentation_.itsmangrove1 = amount; };

    int mangrove2() const { return presentation_.itsmangrove2; };
    void setmangrove2(int amount) { presentation_.itsmangrove2 = amount; };

    int mangrove3() const { return presentation_.itsmangrove3; };
    void setmangrove3(int amount) { presentation_.itsmangrove3 = amount; };

    int highlight1() const { return presentation_.itshighlight1; };
    void sethighlight1(int amount) { presentation_.itshighlight1 = amount; };

    int highlight2() const { return presentation_.itshighlight2; };
    void sethighlight2(int amount) { presentation_.itshighlight2 = amount; };

    int highlight3() const { return presentation_.itshighlight3; };
    void sethighlight3(int amount) { presentation_.itshighlight3 = amount; };

    int outline1() const { return presentation_.itsoutline1; };
    void setoutline1(int amount) { presentation_.itsoutline1 = amount; };

    int outline2() const { return presentation_.itsoutline2; };
    void setoutline2(int amount) { presentation_.itsoutline2 = amount; };

    int outline3() const { return presentation_.itsoutline3; };
    void setoutline3(int amount) { presentation_.itsoutline3 = amount; };

    int elevationlow1() const { return presentation_.itselevationlow1; };
    void setelevationlow1(int amount) { presentation_.itselevationlow1 = amount; };

    int elevationlow2() const { return presentation_.itselevationlow2; };
    void setelevationlow2(int amount) { presentation_.itselevationlow2 = amount; };

    int elevationlow3() const { return presentation_.itselevationlow3; };
    void setelevationlow3(int amount) { presentation_.itselevationlow3 = amount; };

    int elevationhigh1() const { return presentation_.itselevationhigh1; };
    void setelevationhigh1(int amount) { presentation_.itselevationhigh1 = amount; };

    int elevationhigh2() const { return presentation_.itselevationhigh2; };
    void setelevationhigh2(int amount) { presentation_.itselevationhigh2 = amount; };

    int elevationhigh3() const { return presentation_.itselevationhigh3; };
    void setelevationhigh3(int amount) { presentation_.itselevationhigh3 = amount; };

    int temperaturecold1() const { return presentation_.itstemperaturecold1; };
    void settemperaturecold1(int amount) { presentation_.itstemperaturecold1 = amount; };

    int temperaturecold2() const { return presentation_.itstemperaturecold2; };
    void settemperaturecold2(int amount) { presentation_.itstemperaturecold2 = amount; };

    int temperaturecold3() const { return presentation_.itstemperaturecold3; };
    void settemperaturecold3(int amount) { presentation_.itstemperaturecold3 = amount; };

    int temperaturetemperate1() const { return presentation_.itstemperaturetemperate1; };
    void settemperaturetemperate1(int amount) { presentation_.itstemperaturetemperate1 = amount; };

    int temperaturetemperate2() const { return presentation_.itstemperaturetemperate2; };
    void settemperaturetemperate2(int amount) { presentation_.itstemperaturetemperate2 = amount; };

    int temperaturetemperate3() const { return presentation_.itstemperaturetemperate3; };
    void settemperaturetemperate3(int amount) { presentation_.itstemperaturetemperate3 = amount; };

    int temperaturehot1() const { return presentation_.itstemperaturehot1; };
    void settemperaturehot1(int amount) { presentation_.itstemperaturehot1 = amount; };

    int temperaturehot2() const { return presentation_.itstemperaturehot2; };
    void settemperaturehot2(int amount) { presentation_.itstemperaturehot2 = amount; };

    int temperaturehot3() const { return presentation_.itstemperaturehot3; };
    void settemperaturehot3(int amount) { presentation_.itstemperaturehot3 = amount; };

    int precipitationdry1() const { return presentation_.itsprecipitationdry1; };
    void setprecipitationdry1(int amount) { presentation_.itsprecipitationdry1 = amount; };

    int precipitationdry2() const { return presentation_.itsprecipitationdry2; };
    void setprecipitationdry2(int amount) { presentation_.itsprecipitationdry2 = amount; };

    int precipitationdry3() const { return presentation_.itsprecipitationdry3; };
    void setprecipitationdry3(int amount) { presentation_.itsprecipitationdry3 = amount; };

    int precipitationwet1() const { return presentation_.itsprecipitationwet1; };
    void setprecipitationwet1(int amount) { presentation_.itsprecipitationwet1 = amount; };

    int precipitationwet2() const { return presentation_.itsprecipitationwet2; };
    void setprecipitationwet2(int amount) { presentation_.itsprecipitationwet2 = amount; };

    int precipitationwet3() const { return presentation_.itsprecipitationwet3; };
    void setprecipitationwet3(int amount) { presentation_.itsprecipitationwet3 = amount; };

    int climatemapseacolour(int slot, int channel) const;
    void setclimatemapseacolour(int slot, int channel, int amount);

    int climatemapcolour(int slot, int channel) const;
    void setclimatemapcolour(int slot, int channel, int amount);

    int biomemapcolour(int slot, int channel) const;
    void setbiomemapcolour(int slot, int channel, int amount);

    int rivermapcolour(int slot, int channel) const;
    void setrivermapcolour(int slot, int channel, int amount);

    bool showrivermapfeature(int slot) const;
    void setshowrivermapfeature(int slot, bool amount);

    int mapgradientstopcount(int gradient) const;
    void setmapgradientstopcount(int gradient, int count);
    bool mapgradientdiscrete(int gradient) const;
    void setmapgradientdiscrete(int gradient, bool amount);
    int mapgradientposition(int gradient, int stop) const;
    void setmapgradientposition(int gradient, int stop, int amount);
    int mapgradientcolour(int gradient, int stop, int channel) const;
    void setmapgradientcolour(int gradient, int stop, int channel, int amount);

    // These accessor functions are for location-specific information.
    // They don't check that x and y are valid coordinates.

    int map(int x, int y) const; // total terrain elevation, whole metres
    float mapfloat(int x, int y) const; // total terrain elevation

    bool sea(int x, int y) const; // whether this is sea or not

    bool outline(int x, int y) const;   // whether this is coast, next to sea
    // Continental-shelf mask (geosphere.continental_shelf), edited in place by terrain passes.
    worldgrid<bool>& continentalshelves() { UW_FIELD_READ(continental_shelf); return geosphere_.shelfmap; }
    const worldgrid<bool>& continentalshelves() const { UW_FIELD_READ(continental_shelf); return geosphere_.shelfmap; }

    // Routing and fluvial erosion fields (erosion plan, Phase A), written whole by their processes.
    const worldgrid<std::uint8_t>& flowdirection() const { UW_FIELD_READ(flow_direction); return hydrosphere_.flowdirectionmap; }
    worldgrid<std::uint8_t>& flowdirection() { return hydrosphere_.flowdirectionmap; }
    const worldgrid<float>& drainagearea() const { UW_FIELD_READ(drainage_area); return hydrosphere_.drainageareamap; }
    worldgrid<float>& drainagearea() { return hydrosphere_.drainageareamap; }
    const worldgrid<float>& discharge() const { UW_FIELD_READ(discharge); return hydrosphere_.dischargemap; }
    worldgrid<float>& discharge() { return hydrosphere_.dischargemap; }
    const worldgrid<float>& filledelevation() const { UW_FIELD_READ(filled_elevation); return hydrosphere_.filledelevationmap; }
    worldgrid<float>& filledelevation() { return hydrosphere_.filledelevationmap; }
    const worldgrid<std::int32_t>& depressionid() const { UW_FIELD_READ(depression_id); return hydrosphere_.depressionidmap; }
    worldgrid<std::int32_t>& depressionid() { return hydrosphere_.depressionidmap; }
    const worldgrid<float>& erodibility() const { UW_FIELD_READ(erodibility); return geosphere_.erodibilitymap; }
    worldgrid<float>& erodibility() { return geosphere_.erodibilitymap; }
    const worldgrid<float>& erosionrate() const { UW_FIELD_READ(erosion_rate); return geosphere_.erosionratemap; }
    worldgrid<float>& erosionrate() { return geosphere_.erosionratemap; }
    const worldgrid<float>& depositionrate() const { UW_FIELD_READ(deposition_rate); return geosphere_.depositionratemap; }
    worldgrid<float>& depositionrate() { return geosphere_.depositionratemap; }
    const worldgrid<float>& denudation() const { UW_FIELD_READ(denudation); return geosphere_.denudationmap; }
    worldgrid<float>& denudation() { return geosphere_.denudationmap; }
    const worldgrid<float>& aggradation() const { UW_FIELD_READ(aggradation); return geosphere_.aggradationmap; }
    worldgrid<float>& aggradation() { return geosphere_.aggradationmap; }
    float radius() const { UW_FIELD_READ(radius); return body_.itsradius_m; }
    void setradius(float amount) { body_.itsradius_m = amount; }
    float erosionduration() const { UW_FIELD_READ(erosion_duration); return geosphere_.itserosionduration; }
    void seterosionduration(float amount) { geosphere_.itserosionduration = amount; }
    float erodedvolume() const { UW_FIELD_READ(eroded_volume); return geosphere_.itserodedvolume; }
    void seterodedvolume(float amount) { geosphere_.itserodedvolume = amount; }
    float landdepositedvolume() const { UW_FIELD_READ(land_deposited_volume); return geosphere_.itslanddepositedvolume; }
    void setlanddepositedvolume(float amount) { geosphere_.itslanddepositedvolume = amount; }
    float marinedepositedvolume() const { UW_FIELD_READ(marine_deposited_volume); return geosphere_.itsmarinedepositedvolume; }
    void setmarinedepositedvolume(float amount) { geosphere_.itsmarinedepositedvolume = amount; }

    bool coast(int x, int y) const; // whether this is coast, next to land

    void longitude(int x, int& degrees, int& minutes, int& seconds, bool& negative) const;   // returns the longitude of this point

    void latitude(int y, int& degrees, int& minutes, int& seconds, bool& negative) const;   // returns the latitude of this point

    int reverselatitude(int lat) const; // returns the y coordinate of a latitude

    int nom(int x, int y) const; // no-mountains terrain elevation, whole metres
    void setnom(int x, int y, int amount);
    float nomfloat(int x, int y) const; // no-mountains terrain elevation
    void setnomfloat(int x, int y, float amount);
    // Land and sea follow the elevation rounded to whole metres (the rule of sea(), the renderers and the
    // save format): at or below sea level by that rule, for the base elevation and for the whole surface.
    bool basebelowsea(int x, int y) const;
    bool surfacebelowsea(int x, int y) const;

    //int nomridge(int x, int y) const; // no-mountains, with undersea ridge terrain elevation

    int oceanridges(int x, int y) const; // undersea ridges
    void setoceanridges(int x, int y, int amount);

    int oceanridgeheights(int x, int y) const; // undersea ridges
    float oceanridgeheightsfloat(int x, int y) const;
    void setoceanridgeheights(int x, int y, float amount);

    int oceanrifts(int x, int y) const; // undersea rifts
    void setoceanrifts(int x, int y, int amount);

    int oceanridgeoffset(int x, int y) const; // offset map for the ridges
    void setoceanridgeoffset(int x, int y, int amount);

    int oceanridgeangle(int x, int y) const; // the angle perpendicular to the line of the ridge
    void setoceanridgeangle(int x, int y, int amount);

    int volcano(int x, int y) const; // isolated peaks (negative value means it's extinct)
    float volcanofloat(int x, int y) const;
    void setvolcano(int x, int y, float amount);

    bool strato(int x, int y) const; // stratovolcanoes
    void setstrato(int x, int y, bool amount);

    int extraelev(int x, int y) const;   // extra terrain elevation
    float extraelevfloat(int x, int y) const;
    void setextraelev(int x, int y, float amount);

    int maxtemp(int x, int y) const; // maximum temperature
    void setmaxtemp(int x, int y, int amount);

    int mintemp(int x, int y) const; // minimum temperature
    void setmintemp(int x, int y, int amount);

    int jantemp(int x, int y) const; // January temperature
    void setjantemp(int x, int y, int amount);

    int jultemp(int x, int y) const; // July temperature
    void setjultemp(int x, int y, int amount);

    int aprtemp(int x, int y) const; // April temperature
    int octtemp(int x, int y) const; // October temperature

    int avetemp(int x, int y) const; // average temperature

    int seasonaltemp(int season, int x, int y) const;
    void setseasonaltemp(int season, int x, int y, int amount);

    int summerrain(int x, int y) const;  // summer precipitation
    void setsummerrain(int x, int y, int amount);

    int winterrain(int x, int y) const;  // winter precipitation
    void setwinterrain(int x, int y, int amount);

    int janrain(int x, int y) const; // January precipitation
    void setjanrain(int x, int y, int amount);

    int julrain(int x, int y) const; // July precipitation
    void setjulrain(int x, int y, int amount);

    int aprrain(int x, int y) const; // April precipitation
    int octrain(int x, int y) const; // October precipitation

    int averain(int x, int y) const; // average precipitation
    float averainfloat(int x, int y) const;
    // The same values as the whole-number accessors above, without rounding.
    float mintempfloat(int x, int y) const;
    float maxtempfloat(int x, int y) const;
    float avetempfloat(int x, int y) const;
    float summerrainfloat(int x, int y) const;
    float winterrainfloat(int x, int y) const;
    float jantempfloat(int x, int y) const;
    float jultempfloat(int x, int y) const;
    float janrainfloat(int x, int y) const;
    float julrainfloat(int x, int y) const;
    void setsummerrainfloat(int x, int y, float amount);
    void setwinterrainfloat(int x, int y, float amount);
    void setmintempfloat(int x, int y, float amount);
    void setmaxtempfloat(int x, int y, float amount);
    void setjantempfloat(int x, int y, float amount);
    void setjultempfloat(int x, int y, float amount);
    void setjanrainfloat(int x, int y, float amount);
    void setjulrainfloat(int x, int y, float amount);

    int seasonalrain(int season, int x, int y) const;
    void setseasonalrain(int season, int x, int y, int amount);
    float seasonalrainfloat(int season, int x, int y) const;
    float seasonaltempfloat(int season, int x, int y) const;
    void setseasonaltempfloat(int season, int x, int y, float amount);
    void setseasonalrainfloat(int season, int x, int y, float amount);

    int seasonalpressure(int season, int x, int y) const;
    void setseasonalpressure(int season, int x, int y, int amount);

    int seasonaluwind(int season, int x, int y) const;
    void setseasonaluwind(int season, int x, int y, int amount);

    int seasonalvwind(int season, int x, int y) const;
    void setseasonalvwind(int season, int x, int y, int amount);

    int seasonalupperheight(int season, int x, int y) const;
    void setseasonalupperheight(int season, int x, int y, int amount);

    int seasonalupperuwind(int season, int x, int y) const;
    void setseasonalupperuwind(int season, int x, int y, int amount);

    int seasonaluppervwind(int season, int x, int y) const;
    void setseasonaluppervwind(int season, int x, int y, int amount);

    int seasonalverticalvelocity(int season, int x, int y) const;
    void setseasonalverticalvelocity(int season, int x, int y, int amount);

    int seasonalcurrentu(int season, int x, int y) const;
    void setseasonalcurrentu(int season, int x, int y, int amount);

    int seasonalcurrentv(int season, int x, int y) const;
    void setseasonalcurrentv(int season, int x, int y, int amount);

    int seasonalsst(int season, int x, int y) const;
    void setseasonalsst(int season, int x, int y, int amount);
    float seasonalsstfloat(int season, int x, int y) const;

    int seasonalevaporation(int season, int x, int y) const;
    void setseasonalevaporation(int season, int x, int y, int amount);

    int seasonalmaritimeinfluence(int season, int x, int y) const;
    void setseasonalmaritimeinfluence(int season, int x, int y, int amount);

    int seasonalmaritimethermalanomaly(int season, int x, int y) const;
    void setseasonalmaritimethermalanomaly(int season, int x, int y, int amount);

    int seasonalmaritimefetch(int season, int x, int y) const;
    void setseasonalmaritimefetch(int season, int x, int y, int amount);

    int seasonalmoisture(int season, int x, int y) const;
    void setseasonalmoisture(int season, int x, int y, int amount);

    int seasonalconvergence(int season, int x, int y) const;
    void setseasonalconvergence(int season, int x, int y, int amount);

    int seasonaluplift(int season, int x, int y) const;
    void setseasonaluplift(int season, int x, int y, int amount);

    int seasonalsubsidence(int season, int x, int y) const;
    void setseasonalsubsidence(int season, int x, int y, int amount);

    int janmountainrain(int x, int y) const;  // jan precipitation on mountains
    void setjanmountainrain(int x, int y, int amount);

    int julmountainrain(int x, int y) const;  // jul precipitation on mountains
    void setjulmountainrain(int x, int y, int amount);

    int wintermountainrain(int x, int y) const;  // winter precipitation on mountains
    void setwintermountainrain(int x, int y, int amount);

    int summermountainrain(int x, int y) const;  // winter precipitation on mountains
    void setsummermountainrain(int x, int y, int amount);

    int janmountainraindir(int x, int y) const;  // direction of jan precipitation on mountains
    void setjanmountainraindir(int x, int y, int amount);

    int julmountainraindir(int x, int y) const;  // direction of jul precipitation on mountains
    void setjulmountainraindir(int x, int y, int amount);

    int wintermountainraindir(int x, int y) const;  // direction of winter precipitation on mountains
    void setwintermountainraindir(int x, int y, int amount);

    int summermountainraindir(int x, int y) const;  // direction of winter precipitation on mountains
    void setsummermountainraindir(int x, int y, int amount);

    int climate(int x, int y) const; // climate type
    void setclimate(int x, int y, int amount);

    int biome(int x, int y) const; // Holdridge biome slot
    void setbiome(int x, int y, int amount);

    int seaice(int x, int y) const;  // sea ice
    void setseaice(int x, int y, int amount);

    int riverdir(int x, int y) const;    // river flow direction
    void setriverdir(int x, int y, int amount);

    int riverjan(int x, int y) const;    // January river flow volume
    void setriverjan(int x, int y, int amount);

    int riverjul(int x, int y) const;    // July river flow volume
    void setriverjul(int x, int y, int amount);

    int riveraveflow(int x, int y) const; // average river flow

    int wind(int x, int y) const;    // wind direction
    void setwind(int x, int y, int amount);

    int lakesurface(int x, int y) const; // lake surface elevation
    float lakesurfacefloat(int x, int y) const;
    void setlakesurface(int x, int y, float amount);

    int truelake(int x, int y) const; // whether this is a true lake

    float roughness(int x, int y) const;   // roughness
    void setroughness(int x, int y, float amount);

    int mountainridge(int x, int y) const;   // mountain ridge directions
    void setmountainridge(int x, int y, int amount);

    int mountainheight(int x, int y) const;  // mountain elevation
    float mountainheightfloat(int x, int y) const;
    void setmountainheight(int x, int y, float amount);

    int craterrim(int x, int y) const;  // crater rims
    float craterrimfloat(int x, int y) const;
    void setcraterrim(int x, int y, float amount);

    int cratercentre(int x, int y) const;  // crater centre heights
    float cratercentrefloat(int x, int y) const;
    void setcratercentre(int x, int y, float amount);

    int craterx(int n) const; // crater centre x
    void setcraterx(int n, int amount);

    int cratery(int n) const; // crater centre y
    void setcratery(int n, int amount);

    int craterelev(int n) const; // crater peak elevation
    void setcraterelev(int n, int amount);

    int craterradius(int n) const; // crater radius
    void setcraterradius(int n, int amount);

    int tide(int x, int y) const;    // tidal strength
    void settide(int x, int y, int amount);

    int riftlakesurface(int x, int y) const; // rift lake surface elevation
    float riftlakesurfacefloat(int x, int y) const;
    void setriftlakesurface(int x, int y, float amount);

    int riftlakebed(int x, int y) const; // rift lake bed elevation
    float riftlakebedfloat(int x, int y) const;
    void setriftlakebed(int x, int y, float amount);

    int special(int x, int y) const; // special features. 100: salt lake. 110: salt pan. 120: dunes. 130: fresh wetlands. 131: brackish wetlands. 132: salt wetlands.
    void setspecial(int x, int y, int amount);

    GeologicRegime geologicregime(int x, int y) const;
    void setgeologicregime(int x, int y, GeologicRegime amount);

    int tectonicconvergence(int x, int y) const;
    void settectonicconvergence(int x, int y, int amount);

    int tectonicdivergence(int x, int y) const;
    void settectonicdivergence(int x, int y, int amount);

    int tectonicshear(int x, int y) const;
    void settectonicshear(int x, int y, int amount);

    float tectoniccrustagemyr(int x, int y) const;
    void settectoniccrustagemyr(int x, int y, float amount);

    float tectoniccrustthickness(int x, int y) const;
    void settectoniccrustthickness(int x, int y, float amount);

    CrustClass tectoniccrustclass(int x, int y) const;
    void settectoniccrustclass(int x, int y, CrustClass amount);

    float tectonicuplifttendency(int x, int y) const;
    void settectonicuplifttendency(int x, int y, float amount);

    float tectonicsubsidencetendency(int x, int y) const;
    void settectonicsubsidencetendency(int x, int y, float amount);

    float tectonicaccumulatedstrain(int x, int y) const;
    void settectonicaccumulatedstrain(int x, int y, float amount);

    BoundaryType tectonicboundarytype(int x, int y) const;
    void settectonicboundarytype(int x, int y, BoundaryType amount);

    int tectonicboundarydistance(int x, int y) const;
    void settectonicboundarydistance(int x, int y, int amount);

    int tectonicboundarysegmentid(int x, int y) const;
    void settectonicboundarysegmentid(int x, int y, int amount);

    int tectonicnearestboundaryid(int x, int y) const;
    void settectonicnearestboundaryid(int x, int y, int amount);

    float tectonicboundaryhistory(int x, int y) const;
    void settectonicboundaryhistory(int x, int y, float amount);

    int tectonicdeformingregionid(int x, int y) const;
    void settectonicdeformingregionid(int x, int y, int amount);

    DeformingRegionType tectonicdeformingregiontype(int x, int y) const;
    void settectonicdeformingregiontype(int x, int y, DeformingRegionType amount);

    float tectonicdeformationrate(int x, int y) const;
    void settectonicdeformationrate(int x, int y, float amount);

    float tectonicdeformationvelocityx(int x, int y) const;
    void settectonicdeformationvelocityx(int x, int y, float amount);

    float tectonicdeformationvelocityy(int x, int y) const;
    void settectonicdeformationvelocityy(int x, int y, float amount);

    const std::vector<TectonicBoundarySegment>& tectonicboundarysegments() const;
    std::vector<TectonicBoundarySegment>& tectonicboundarysegments();
    int tectonicboundarysegmentcount() const;
    // The engine's own fields of one cell and the record of its last update (geosphere.tectonic_snapshot).
    void settectonicsnapshotcell(int x, int y, float height, int plate, int crustagesteps, float continentalcrust,
        float oceaniccrust, float sediment, float riftdamage, float plumeinfluence);
    void settectonicsnapshotrecord(std::string record);
    void settectonicboundarysegments(std::vector<TectonicBoundarySegment> segments);
    template<typename SegmentT> void settectonicboundarysegments(const std::vector<SegmentT>& segments);

    const std::vector<TectonicDeformingRegion>& tectonicdeformingregions() const;
    std::vector<TectonicDeformingRegion>& tectonicdeformingregions();
    int tectonicdeformingregioncount() const;
    void settectonicdeformingregions(std::vector<TectonicDeformingRegion> regions);
    template<typename RegionT> void settectonicdeformingregions(const std::vector<RegionT>& regions);

    const TectonicBoundarySegment* findtectonicboundarysegment(int id) const;
    const TectonicDeformingRegion* findtectonicdeformingregion(int id) const;

    // Deposits recorded through the tectonic run, with placers added on the
    // final world. Positions are fractions of the global map.
    const std::vector<mineralsystems::Deposit>& mineraldeposits() const;
    std::vector<mineralsystems::Deposit>& mineraldeposits();
    void setmineraldeposits(std::vector<mineralsystems::Deposit> deposits);
    const mineralsystems::MineralBudget& mineralbudget() const;
    void setmineralbudget(const mineralsystems::MineralBudget& budget);

    BasinClass basinclass(int x, int y) const;
    void setbasinclass(int x, int y, BasinClass amount);

    int erosionpotential(int x, int y) const;
    void seterosionpotential(int x, int y, int amount);

    int depositionpotential(int x, int y) const;
    void setdepositionpotential(int x, int y, int amount);

    int floodplainfertility(int x, int y) const;
    void setfloodplainfertility(int x, int y, int amount);

    int metalorereserve(int x, int y) const;
    void setmetalorereserve(int x, int y, int amount);

    int placerreserve(int x, int y) const;
    void setplacerreserve(int x, int y, int amount);

    int evaporitereserve(int x, int y) const;
    void setevaporitereserve(int x, int y, int amount);

    int volcanicreserve(int x, int y) const;
    void setvolcanicreserve(int x, int y, int amount);

    int fisheryreserve(int x, int y) const;
    void setfisheryreserve(int x, int y, int amount);

    int settlementsuitability(int x, int y) const;
    void setsettlementsuitability(int x, int y, int amount);

    int infrastructure(int x, int y) const;
    void setinfrastructure(int x, int y, int amount);

    int agriculturalcapacity(int x, int y) const;
    void setagriculturalcapacity(int x, int y, int amount);

    int routetraffic(int x, int y) const;
    void setroutetraffic(int x, int y, int amount);

    int riveraccess(int x, int y) const;
    void setriveraccess(int x, int y, int amount);

    int harborscore(int x, int y) const;
    void setharborscore(int x, int y, int amount);

    int ownersettlementid(int x, int y) const;
    void setownersettlementid(int x, int y, int amount);

    int ownerpolityid(int x, int y) const;
    void setownerpolityid(int x, int y, int amount);

    const std::vector<Settlement>& settlements() const;
    std::vector<Settlement>& settlements();

    const std::vector<Polity>& polities() const;
    std::vector<Polity>& polities();

    const std::vector<RouteEdge>& routeedges() const;
    std::vector<RouteEdge>& routeedges();

    const std::vector<HistoryEvent>& historyevents() const;
    std::vector<HistoryEvent>& historyevents();

    void clearsocialstate();

    int deltadir(int x, int y) const;    // delta branch flow direction (reversed)
    int deltadir_no_bounds_check(int x, int y) const { return hydrosphere_.deltamapdir[x][y]; }
    void setdeltadir(int x, int y, int amount);

    int deltajan(int x, int y) const;    // January delta branch flow volume
    void setdeltajan(int x, int y, int amount);

    int deltajul(int x, int y) const;    // July delta branch flow volume
    void setdeltajul(int x, int y, int amount);

    bool lakestart(int x, int y) const; //  Whether a rift lake starts here
    void setlakestart(int x, int y, int amount);

    bool island(int x, int y) const; // Whether this is a one-tile island
    void setisland(int x, int y, bool amount);

    bool mountainisland(int x, int y) const; // Whether this is an island that's basically a mountain in the sea (if it is, the elevation of what would be the sea bed here is recorded)
    void setmountainisland(int x, int y, bool amount);

    bool noshade(int x, int y) const; // Whether no shading should be applied here on the global map.
    void setnoshade(int x, int y, bool amount);

    int noise(int x, int y) const; // The noise map - used for adding the terraced decoration to parts of the regional map.
    void setnoise(int x, int y, int amount);

    int horse(int x, int y) const;   // horse latitudes
    void sethorse(int x, int y, int amount);

    int test(int x, int y) const;    // Test array.
    void settest(int x, int y, int amount);

    // Now we have wrapped versions of all of those.
    // For these, the function will automatically wrap the X coordinate and clip the Y coordinate if need be.

    bool seawrap(int x, int y) const;
    bool outlinewrap(int x, int y) const;
    int mountainheightwrap(int x, int y) const;

    // Raw grid access for accelerated renderer backends.
    std::vector<int> wholemapnom() const;   // base elevation in whole metres, storage order
    const float* rawlakemap() const;
    const float* rawoceanridgeheightmap() const;
    const float* rawmountainheights() const;
    const float* rawvolcanomap() const;
    const float* rawextraelevmap() const;
    const float* rawcraterrims() const;
    const float* rawcratercentres() const;
    const short* rawoceanridgemap() const;
    const short* rawoceanriftmap() const;
    std::vector<short> legacyjantempmap() const;
    std::vector<short> legacyjultempmap() const;
    std::vector<short> legacyjanrainmap() const;
    std::vector<short> legacyjulrainmap() const;
    const short* rawseaicemap() const;
    const short* rawclimatemap() const;
    const short* rawbiomemap() const;
    const int* rawrivermapjan() const;
    const int* rawrivermapjul() const;
    const short* rawspecials() const;
    const short* rawdeltamapdir() const;
    const int* rawdeltamapjan() const;
    const int* rawdeltamapjul() const;
    const float* rawriftlakemapsurface() const;
    const bool* rawnoshademap() const;

    // Other public functions.

    void clear();   // Clears all of the maps.
    void smoothnom(int amount); // Smoothes the no-mountain map to a given amount.
    void smoothextraelev(int amount); // Smoothes the extra elevation map to a given amount.
    void shiftterrain(int offset); // Shifts the physical terrain by a given amount.
    void setmaxriverflow(); // Calculates the largest river flow on the map.
    void syncseasonalclimatefromlegacy(); // Populates explicit seasonal fields from the legacy Jan/Jul model.
    void cleartectonicprovenance();
    const std::vector<std::uint8_t>& weatheranomalystate() const;
    void setweatheranomalystate(std::vector<std::uint8_t> state);
    // Advances only the replayable anomaly state, never the seasonal maps.
    // A zero resolution preserves the saved LOD; failed steps leave it intact.
    bool advanceweather(float elapsedSeconds, int horizontalCells = 0);

    void saveworld(string filename); // Saves the world.
    bool loadworld(string filename); // Loads the world.


private:

    // Private variables.

    int itssaveversion; // The save version. Only saves of this version can be loaded in.
    int itssettingssaveversion; // As above, but for settings files.

    // State by world component; climate state is climate_ (climateworld).
    planetbodystate body_;
    geospherestate geosphere_;
    hydrospherestate hydrosphere_;
    cryospherestate cryosphere_;
    atmosphereparameters atmosphere_;
    biospherestate biosphere_;
    anthropospherestate anthroposphere_;
    presentationsettings presentation_;

    // reused temporary state
    string line_for_file_read;

    // Private functions.

    int wrapx(int x) const;
    int clipy(int y) const;
    void smoothoverland(worldgrid<int>& arr, int amount, bool vary);
    void smoothoverland(worldgrid<float>& arr, int amount, bool vary);

    template<typename T> void shift(worldgrid<T>& arr, int offset);
    template<typename T> void writevariable(ofstream& outfile, T val);
    template<typename T> void readvariable(ifstream& infile, T& val);
    void writeshortvectordata(ofstream& outfile, const std::vector<short>& arr);
    void readshortvectordata(ifstream& infile, std::vector<short>& arr);
    void readlegacyvectordata(ifstream& infile, std::vector<float>& arr, float scale = 1.0f);
    void writelegacydata(ofstream& outfile, const worldgrid<float>& arr);
    void readlegacydata(ifstream& infile, worldgrid<float>& arr);
    void readwholedata(ifstream& infile, worldgrid<float>& arr);
    void readwholeshortdata(ifstream& infile, worldgrid<float>& arr);
    float addedelevation(int x, int y) const;
    void writemineraldeposits(ofstream& outfile);
    bool readmineraldeposits(ifstream& infile);
    int seasonalclimateindex(int x, int y) const;
    bool validseasonindex(int season) const;
    void resizegrids();
    void cleartectonicprovenanceinternal();
    climateworld climate_; // atmosphere, ocean surface, sea ice and biome state
    climateinputs climateinputsview() const;
    template<typename T> void writedata(ofstream& outfile, const worldgrid<T>& arr);
    template<typename T> void readdata(ifstream& infile, worldgrid<T>& arr);
};


inline int planet::saveversion() const { return itssaveversion; }
inline void planet::setsaveversion(int amount) { itssaveversion = amount; }

inline int planet::settingssaveversion() const { return itssettingssaveversion; }
inline void planet::setsettingssaveversion(int amount) { itssettingssaveversion = amount; }

inline int planet::tectonictimeoriginstep() const { UW_FIELD_READ(tectonic_time_origin_step); return geosphere_.itstectonictimeoriginstep; }
inline void planet::settectonictimeoriginstep(int amount) { geosphere_.itstectonictimeoriginstep = (std::max)(0, amount); }

inline float planet::tectonictimemyr() const { UW_FIELD_READ(tectonic_time); return geosphere_.itstectonictimemyr; }
inline void planet::settectonictimemyr(float amount) { geosphere_.itstectonictimemyr = (std::max)(0.0f, amount); }

inline float planet::tectonicdeltatimemyr() const { UW_FIELD_READ(tectonic_time_step); return geosphere_.itstectonicdeltatimemyr; }
inline void planet::settectonicdeltatimemyr(float amount) { geosphere_.itstectonicdeltatimemyr = (std::max)(0.0f, amount); }


inline int planet::tectonicplatecount() const { UW_FIELD_READ(tectonic_plate_count); return geosphere_.itstectonicplatecount; }
inline void planet::settectonicplatecount(int amount) { geosphere_.itstectonicplatecount = (std::max)(0, amount); }

inline int planet::tectonicsealevelm() const { UW_FIELD_READ(tectonic_sea_level); return geosphere_.itstectonicsealevelm; }
inline void planet::settectonicsealevelm(int amount) { geosphere_.itstectonicsealevelm = (std::max)(0, amount); }

inline long planet::seed() const { UW_FIELD_READ(seed); return body_.itsseed; }
inline void planet::setseed(long amount) { body_.itsseed = amount; }

inline int planet::size() const { UW_FIELD_READ(size_class); return body_.itssize; }
inline void planet::setsize(int amount) { body_.itssize = amount; }

inline int planet::width() const { return body_.itswidth; }
inline void planet::setwidth(int amount)
{
    const int columns = std::clamp(amount + 1, 2, ARRAYWIDTH) & ~1;
    body_.itswidth = columns - 1;
    body_.itsheight = columns / 2 - 1;
    resizegrids();
}

inline int planet::height() const { return body_.itsheight; }
inline void planet::setheight(int amount)
{
    if (amount == (body_.itswidth + 1) / 2 - 1)
        body_.itsheight = amount;
    resizegrids();
}

inline bool planet::rotation() const { UW_FIELD_READ(prograde_rotation); return body_.itsrotation; }
inline void planet::setrotation(bool amount) { body_.itsrotation = amount; }

inline float planet::tilt() const { UW_FIELD_READ(axial_tilt); return body_.itstilt; }
inline void planet::settilt(float amount) { body_.itstilt = amount; }

inline float planet::eccentricity() const { UW_FIELD_READ(eccentricity); return body_.itseccentricity; }
inline void planet::seteccentricity(float amount) { body_.itseccentricity = amount; }

inline int planet::perihelion() const { UW_FIELD_READ(perihelion); return (int)body_.itsperihelion; }
inline void planet::setperihelion(int amount) { body_.itsperihelion = (short)amount; }

inline float planet::gravity() const { UW_FIELD_READ(gravity); return body_.itsgravity; }
inline void planet::setgravity(float amount) { body_.itsgravity = amount; }

inline float planet::lunar() const { UW_FIELD_READ(lunar_strength); return body_.itslunar; }
inline void planet::setlunar(float amount) { body_.itslunar = amount; }

inline float planet::tempdecrease() const { UW_FIELD_READ(lapse_rate); return atmosphere_.itstempdecrease; }
inline void planet::settempdecrease(float amount) { atmosphere_.itstempdecrease = amount; }

inline int planet::northpolaradjust() const { UW_FIELD_READ(north_polar_adjust); return atmosphere_.itsnorthpolaradjust; }
inline void planet::setnorthpolaradjust(int amount) { atmosphere_.itsnorthpolaradjust = amount; }

inline int planet::southpolaradjust() const { UW_FIELD_READ(south_polar_adjust); return atmosphere_.itssouthpolaradjust; }
inline void planet::setsouthpolaradjust(int amount) { atmosphere_.itssouthpolaradjust = amount; }

inline int planet::averagetemp() const { UW_FIELD_READ(average_temperature); return atmosphere_.itsaveragetemp; }
inline void planet::setaveragetemp(int amount) { atmosphere_.itsaveragetemp = amount; }

inline int planet::northpolartemp() const { return climate_.northpolartemp(); }
inline void planet::setnorthpolartemp(int amount) { climate_.setnorthpolartemp(amount); }

inline int planet::southpolartemp() const { return climate_.southpolartemp(); }
inline void planet::setsouthpolartemp(int amount) { climate_.setsouthpolartemp(amount); }

inline int planet::eqtemp() const { return climate_.eqtemp(); }
inline void planet::seteqtemp(int amount) { climate_.seteqtemp(amount); }

inline float planet::waterpickup() const { UW_FIELD_READ(water_pickup); return atmosphere_.itswaterpickup; }
inline void planet::setwaterpickup(float amount) { atmosphere_.itswaterpickup = amount; }

inline float planet::riverfactor() const { UW_FIELD_READ(river_factor); return hydrosphere_.itsriverfactor; }
inline void planet::setriverfactor(float amount) { hydrosphere_.itsriverfactor = amount; }

inline int planet::riverlandreduce() const { UW_FIELD_READ(river_land_reduce); return hydrosphere_.itsriverlandreduce; }
inline void planet::setriverlandreduce(int amount) { hydrosphere_.itsriverlandreduce = amount; }

inline int planet::estuarylimit() const { UW_FIELD_READ(estuary_limit); return hydrosphere_.itsestuarylimit; }
inline void planet::setestuarylimit(int amount) { hydrosphere_.itsestuarylimit = amount; }

inline int planet::glacialtemp() const { UW_FIELD_READ(glacial_temperature); return cryosphere_.itsglacialtemp; }
inline void planet::setglacialtemp(int amount) { cryosphere_.itsglacialtemp = amount; }

inline int planet::glaciertemp() const { UW_FIELD_READ(glacier_temperature); return cryosphere_.itsglaciertemp; }
inline void planet::setglaciertemp(int amount) { cryosphere_.itsglaciertemp = amount; }

inline float planet::mountainreduce() const { UW_FIELD_READ(mountain_reduce); return geosphere_.itsmountainreduce; }
inline void planet::setmountainreduce(float amount) { geosphere_.itsmountainreduce = amount; }

inline int planet::climatenumber() const { UW_FIELD_READ(climate_type_count); return atmosphere_.itsclimateno; }
inline void planet::setclimatenumber(int amount) { atmosphere_.itsclimateno = amount; }

inline int planet::maxelevation() const { UW_FIELD_READ(max_elevation); return geosphere_.itsmaxheight; }
inline void planet::setmaxelevation(int amount) { geosphere_.itsmaxheight = amount; }

inline int planet::sealevel() const { UW_FIELD_READ(sea_level); return geosphere_.itssealevel; }
inline void planet::setsealevel(int amount) { geosphere_.itssealevel = amount; }

inline float planet::landshading() const { return presentation_.itslandshading; };
inline void planet::setlandshading(float amount) { presentation_.itslandshading = amount; };

inline float planet::lakeshading() const { return presentation_.itslakeshading; };
inline void planet::setlakeshading(float amount) { presentation_.itslakeshading = amount; };

inline float planet::seashading() const { return presentation_.itsseashading; };
inline void planet::setseashading(float amount) { presentation_.itsseashading = amount; };

inline int planet::shadingdir() const { return (int)presentation_.itsshadingdir; };
inline void planet::setshadingdir(int amount) { presentation_.itsshadingdir = (short)amount; };

inline int planet::snowchange() const { return (int)presentation_.itssnowchange; };
inline void planet::setsnowchange(int amount) { presentation_.itssnowchange = (short)amount; };

inline int planet::seaiceappearance() const { return (int)presentation_.itsseaiceappearance; };
inline void planet::setseaiceappearance(int amount) { presentation_.itsseaiceappearance = (short)amount; };

inline bool planet::colourcliffs() const { return presentation_.itscolourcliffs; };
inline void planet::setcolourcliffs(bool amount) { presentation_.itscolourcliffs = amount; };
inline bool planet::showmapoutline() const { return presentation_.itsshowmapoutline; };
inline void planet::setshowmapoutline(bool amount) { presentation_.itsshowmapoutline = amount; };

inline float planet::landmarbling() const { return presentation_.itslandmarbling; };
inline void planet::setlandmarbling(float amount) { presentation_.itslandmarbling = amount; };

inline float planet::lakemarbling() const { return presentation_.itslakemarbling; };
inline void planet::setlakemarbling(float amount) { presentation_.itslakemarbling = amount; };

inline float planet::seamarbling() const { return presentation_.itsseamarbling; };
inline void planet::setseamarbling(float amount) { presentation_.itsseamarbling = amount; };

inline int planet::minriverflowglobal() const { return presentation_.itsminriverflowglobal; };
inline void planet::setminriverflowglobal(int amount) { presentation_.itsminriverflowglobal = amount; };

inline int planet::minriverflowregional() const { return presentation_.itsminriverflowregional; };
inline void planet::setminriverflowregional(int amount) { presentation_.itsminriverflowregional = amount; };

inline bool planet::showmangroves() const { return presentation_.itsmangroves; };
inline void planet::setshowmangroves(bool amount) { presentation_.itsmangroves = amount; };

inline int planet::climatemapseacolour(int slot, int channel) const
{
    if (slot < 0 || slot >= CLIMATEMAPSEACOLOURCOUNT || channel < 0 || channel > 2)
        return 0;

    return presentation_.itsclimatemapseacolours[slot][channel];
}

inline void planet::setclimatemapseacolour(int slot, int channel, int amount)
{
    if (slot < 0 || slot >= CLIMATEMAPSEACOLOURCOUNT || channel < 0 || channel > 2)
        return;

    presentation_.itsclimatemapseacolours[slot][channel] = amount;
}

inline int planet::climatemapcolour(int slot, int channel) const
{
    if (slot < 0 || slot >= CLIMATEMAPCOLOURCOUNT || channel < 0 || channel > 2)
        return 0;

    return presentation_.itsclimatemapcolours[slot][channel];
}

inline void planet::setclimatemapcolour(int slot, int channel, int amount)
{
    if (slot < 0 || slot >= CLIMATEMAPCOLOURCOUNT || channel < 0 || channel > 2)
        return;

    presentation_.itsclimatemapcolours[slot][channel] = amount;
}

inline int planet::biomemapcolour(int slot, int channel) const
{
    if (slot < 0 || slot >= BIOMEMAPCOLOURCOUNT || channel < 0 || channel > 2)
        return 0;

    return presentation_.itsbiomemapcolours[slot][channel];
}

inline void planet::setbiomemapcolour(int slot, int channel, int amount)
{
    if (slot < 0 || slot >= BIOMEMAPCOLOURCOUNT || channel < 0 || channel > 2)
        return;

    presentation_.itsbiomemapcolours[slot][channel] = amount;
}

inline int planet::rivermapcolour(int slot, int channel) const
{
    if (slot < 0 || slot >= RIVERMAPCOLOURCOUNT || channel < 0 || channel > 2)
        return 0;

    return presentation_.itsrivermapcolours[slot][channel];
}

inline void planet::setrivermapcolour(int slot, int channel, int amount)
{
    if (slot < 0 || slot >= RIVERMAPCOLOURCOUNT || channel < 0 || channel > 2)
        return;

    presentation_.itsrivermapcolours[slot][channel] = amount;
}

inline bool planet::showrivermapfeature(int slot) const
{
    if (slot < 0 || slot >= RIVERMAPFEATURECOUNT)
        return false;

    return presentation_.itsshowrivermapfeatures[slot];
}

inline void planet::setshowrivermapfeature(int slot, bool amount)
{
    if (slot < 0 || slot >= RIVERMAPFEATURECOUNT)
        return;

    presentation_.itsshowrivermapfeatures[slot] = amount;
}

inline int planet::mapgradientstopcount(int gradient) const
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT)
        return 0;

    return presentation_.itsmapgradientstopcounts[gradient];
}

inline void planet::setmapgradientstopcount(int gradient, int count)
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT)
        return;

    if (count < 0)
        count = 0;
    if (count > MAPGRADIENTMAXSTOPS)
        count = MAPGRADIENTMAXSTOPS;

    presentation_.itsmapgradientstopcounts[gradient] = count;
}

inline bool planet::mapgradientdiscrete(int gradient) const
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT)
        return false;

    return presentation_.itsmapgradientdiscrete[gradient];
}

inline void planet::setmapgradientdiscrete(int gradient, bool amount)
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT)
        return;

    presentation_.itsmapgradientdiscrete[gradient] = amount;
}

inline int planet::mapgradientposition(int gradient, int stop) const
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT || stop < 0 || stop >= MAPGRADIENTMAXSTOPS)
        return 0;

    return presentation_.itsmapgradientpositions[gradient][stop];
}

inline void planet::setmapgradientposition(int gradient, int stop, int amount)
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT || stop < 0 || stop >= MAPGRADIENTMAXSTOPS)
        return;

    presentation_.itsmapgradientpositions[gradient][stop] = amount;
}

inline int planet::mapgradientcolour(int gradient, int stop, int channel) const
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT || stop < 0 || stop >= MAPGRADIENTMAXSTOPS || channel < 0 || channel > 2)
        return 0;

    return presentation_.itsmapgradientcolours[gradient][stop][channel];
}

inline void planet::setmapgradientcolour(int gradient, int stop, int channel, int amount)
{
    if (gradient < 0 || gradient >= MAPGRADIENTTYPECOUNT || stop < 0 || stop >= MAPGRADIENTMAXSTOPS || channel < 0 || channel > 2)
        return;

    presentation_.itsmapgradientcolours[gradient][stop][channel] = amount;
}

inline const float* planet::rawlakemap() const { UW_FIELD_READ(lake_surface); return &hydrosphere_.lakemap[0][0]; }
inline const float* planet::rawoceanridgeheightmap() const { UW_FIELD_READ(ocean_ridge_height); return &geosphere_.oceanridgeheightmap[0][0]; }
inline const float* planet::rawmountainheights() const { UW_FIELD_READ(mountain_height); return &geosphere_.mountainheights[0][0]; }
inline const float* planet::rawvolcanomap() const { UW_FIELD_READ(volcano); return &geosphere_.volcanomap[0][0]; }
inline const float* planet::rawextraelevmap() const { UW_FIELD_READ(extra_elevation); return &geosphere_.extraelevmap[0][0]; }
inline const float* planet::rawcraterrims() const { UW_FIELD_READ(crater_rims); return &geosphere_.craterrims[0][0]; }
inline const float* planet::rawcratercentres() const { UW_FIELD_READ(crater_centres); return &geosphere_.cratercentres[0][0]; }
inline const short* planet::rawoceanridgemap() const { UW_FIELD_READ(ocean_ridge); return &geosphere_.oceanridgemap[0][0]; }
inline const short* planet::rawoceanriftmap() const { UW_FIELD_READ(ocean_rift); return &geosphere_.oceanriftmap[0][0]; }
inline std::vector<short> planet::legacyjantempmap() const { return climate_.legacyjantempmap(); }
inline std::vector<short> planet::legacyjultempmap() const { return climate_.legacyjultempmap(); }
inline std::vector<short> planet::legacyjanrainmap() const { return climate_.legacyjanrainmap(); }
inline std::vector<short> planet::legacyjulrainmap() const { return climate_.legacyjulrainmap(); }
inline const short* planet::rawseaicemap() const { return climate_.rawseaicemap(); }
inline const short* planet::rawclimatemap() const { return climate_.rawclimatemap(); }
inline const short* planet::rawbiomemap() const { return climate_.rawbiomemap(); }
inline const int* planet::rawrivermapjan() const { UW_FIELD_READ(river_flow_january); return &hydrosphere_.rivermapjan[0][0]; }
inline const int* planet::rawrivermapjul() const { UW_FIELD_READ(river_flow_july); return &hydrosphere_.rivermapjul[0][0]; }
inline const short* planet::rawspecials() const { UW_FIELD_READ(special_features); return &hydrosphere_.specials[0][0]; }
inline const short* planet::rawdeltamapdir() const { UW_FIELD_READ(delta_direction); return &hydrosphere_.deltamapdir[0][0]; }
inline const int* planet::rawdeltamapjan() const { UW_FIELD_READ(delta_flow_january); return &hydrosphere_.deltamapjan[0][0]; }
inline const int* planet::rawdeltamapjul() const { UW_FIELD_READ(delta_flow_july); return &hydrosphere_.deltamapjul[0][0]; }
inline const float* planet::rawriftlakemapsurface() const { UW_FIELD_READ(rift_lake_surface); return &hydrosphere_.riftlakemapsurface[0][0]; }
inline const bool* planet::rawnoshademap() const { UW_FIELD_READ(no_shade); return &presentation_.noshademap[0][0]; }

inline int planet::seasonalclimateindex(int x, int y) const
{
    return y * (body_.itswidth + 1) + x;
}

inline bool planet::validseasonindex(int season) const
{
    return season >= 0 && season < CLIMATESEASONCOUNT;
}

inline int planet::landtotal() const { return geosphere_.itslandtotal; };
inline void planet::setlandtotal(int amount) { geosphere_.itslandtotal = amount; };

inline int planet::seatotal() const { return geosphere_.itsseatotal; };
inline void planet::setseatotal(int amount) { geosphere_.itsseatotal = amount; };

inline int planet::craterno() const { return geosphere_.itscraterno; };
inline void planet::setcraterno(int amount) { geosphere_.itscraterno = amount; };

inline int planet::maxriverflow() const { UW_FIELD_READ(max_river_flow); return hydrosphere_.itsmaxriverflow; }

inline int planet::map(int x, int y) const
{
    UW_FIELD_READ(volcano);
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(ocean_ridge_height);
    UW_FIELD_READ(mountain_height);
    UW_FIELD_READ(extra_elevation);
    UW_FIELD_READ(crater_rims);
    UW_FIELD_READ(crater_centres);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.mapnom[x][y]) + wholemetres(addedelevation(x, y));
}

inline float planet::mapfloat(int x, int y) const
{
    UW_FIELD_READ(volcano);
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(ocean_ridge_height);
    UW_FIELD_READ(mountain_height);
    UW_FIELD_READ(extra_elevation);
    UW_FIELD_READ(crater_rims);
    UW_FIELD_READ(crater_centres);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.mapnom[x][y] + addedelevation(x, y);
}

// Ridges, volcanoes, mountains, canyon uplift and craters on top of the base
// elevation. Land and sea follow the base elevation in whole metres.
inline float planet::addedelevation(int x, int y) const
{
    const float thisvolcano = std::fabs(geosphere_.volcanomap[x][y]);

    if (wholemetres(geosphere_.mapnom[x][y]) <= geosphere_.itssealevel)
        return geosphere_.oceanridgeheightmap[x][y] + thisvolcano;

    if (thisvolcano > geosphere_.mountainheights[x][y])
        return thisvolcano + geosphere_.extraelevmap[x][y] + geosphere_.craterrims[x][y] + geosphere_.cratercentres[x][y];
    else
        return geosphere_.mountainheights[x][y] + geosphere_.extraelevmap[x][y] + geosphere_.craterrims[x][y] + geosphere_.cratercentres[x][y];
}

inline int planet::nom(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.mapnom[x][y]);
}

inline void planet::setnom(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.mapnom[x][y] = static_cast<float>(amount);
}

inline bool planet::basebelowsea(int x, int y) const
{
    UW_FIELD_READ(sea_level);
    return nom(x, y) <= geosphere_.itssealevel;
}

inline bool planet::surfacebelowsea(int x, int y) const
{
    return map(x, y) <= geosphere_.itssealevel;
}

inline float planet::nomfloat(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.mapnom[x][y];
}

inline void planet::setnomfloat(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.mapnom[x][y] = amount;
}

inline int planet::oceanridges(int x, int y) const
{
    UW_FIELD_READ(ocean_ridge);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)geosphere_.oceanridgemap[x][y];
}

inline void planet::setoceanridges(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.oceanridgemap[x][y] = (short)amount;
}

inline int planet::oceanridgeheights(int x, int y) const
{
    UW_FIELD_READ(ocean_ridge_height);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.oceanridgeheightmap[x][y]);
}

inline float planet::oceanridgeheightsfloat(int x, int y) const
{
    UW_FIELD_READ(ocean_ridge_height);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.oceanridgeheightmap[x][y];
}

inline void planet::setoceanridgeheights(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.oceanridgeheightmap[x][y] = amount;
}

inline int planet::oceanrifts(int x, int y) const
{
    UW_FIELD_READ(ocean_rift);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)geosphere_.oceanriftmap[x][y];
}

inline void planet::setoceanrifts(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.oceanriftmap[x][y] = (short)amount;
}

inline int planet::oceanridgeoffset(int x, int y) const
{
    UW_FIELD_READ(ocean_ridge_offset);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)geosphere_.oceanridgeoffsetmap[x][y];
}

inline void planet::setoceanridgeoffset(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.oceanridgeoffsetmap[x][y] = (short)amount;
}

inline int planet::oceanridgeangle(int x, int y) const
{
    UW_FIELD_READ(ocean_ridge_angle);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)geosphere_.oceanridgeanglemap[x][y];
}

inline void planet::setoceanridgeangle(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.oceanridgeanglemap[x][y] = (short)amount;
}

inline int planet::volcano(int x, int y) const
{
    UW_FIELD_READ(volcano);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.volcanomap[x][y]);
}

inline float planet::volcanofloat(int x, int y) const
{
    UW_FIELD_READ(volcano);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.volcanomap[x][y];
}

inline void planet::setvolcano(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.volcanomap[x][y] = amount;
}

inline bool planet::strato(int x, int y) const
{
    UW_FIELD_READ(stratovolcano);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.stratomap[x][y];
}

inline void planet::setstrato(int x, int y, bool amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.stratomap[x][y] = amount;
}

inline int planet::extraelev(int x, int y) const
{
    UW_FIELD_READ(extra_elevation);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.extraelevmap[x][y]);
}

inline float planet::extraelevfloat(int x, int y) const
{
    UW_FIELD_READ(extra_elevation);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.extraelevmap[x][y];
}

inline void planet::setextraelev(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.extraelevmap[x][y] = amount;
}

inline int planet::jantemp(int x, int y) const { return climate_.jantemp(x, y); }

inline void planet::setjantemp(int x, int y, int amount) { climate_.setjantemp(x, y, amount); }

inline int planet::jultemp(int x, int y) const { return climate_.jultemp(x, y); }

inline void planet::setjultemp(int x, int y, int amount) { climate_.setjultemp(x, y, amount); }

inline int planet::seasonaltemp(int season, int x, int y) const { return climate_.seasonaltemp(season, x, y); }

inline void planet::setseasonaltemp(int season, int x, int y, int amount) { climate_.setseasonaltemp(season, x, y, amount); }

inline int planet::aprtemp(int x, int y) const { return climate_.aprtemp(x, y); }

inline int planet::octtemp(int x, int y) const { return climate_.octtemp(x, y); }


inline int planet::avetemp(int x, int y) const { return climate_.avetemp(x, y); }

inline int planet::janrain(int x, int y) const { return climate_.janrain(x, y); }

inline void planet::setjanrain(int x, int y, int amount) { climate_.setjanrain(x, y, amount); }

inline int planet::julrain(int x, int y) const { return climate_.julrain(x, y); }

inline void planet::setjulrain(int x, int y, int amount) { climate_.setjulrain(x, y, amount); }

inline int planet::seasonalrain(int season, int x, int y) const { return climate_.seasonalrain(season, x, y); }

inline void planet::setseasonalrain(int season, int x, int y, int amount) { climate_.setseasonalrain(season, x, y, amount); }

inline float planet::seasonalrainfloat(int season, int x, int y) const { return climate_.seasonalrainfloat(season, x, y); }

inline float planet::seasonaltempfloat(int season, int x, int y) const { return climate_.seasonaltempfloat(season, x, y); }

inline void planet::setseasonaltempfloat(int season, int x, int y, float amount) { climate_.setseasonaltempfloat(season, x, y, amount); }

inline void planet::setseasonalrainfloat(int season, int x, int y, float amount) { climate_.setseasonalrainfloat(season, x, y, amount); }

inline int planet::aprrain(int x, int y) const { return climate_.aprrain(x, y); }

inline int planet::octrain(int x, int y) const { return climate_.octrain(x, y); }

inline int planet::averain(int x, int y) const { return climate_.averain(x, y); }

inline float planet::averainfloat(int x, int y) const { return climate_.averainfloat(x, y); }

inline int planet::seasonalpressure(int season, int x, int y) const { return climate_.seasonalpressure(season, x, y); }

inline void planet::setseasonalpressure(int season, int x, int y, int amount) { climate_.setseasonalpressure(season, x, y, amount); }

inline int planet::seasonaluwind(int season, int x, int y) const { return climate_.seasonaluwind(season, x, y); }

inline void planet::setseasonaluwind(int season, int x, int y, int amount) { climate_.setseasonaluwind(season, x, y, amount); }

inline int planet::seasonalvwind(int season, int x, int y) const { return climate_.seasonalvwind(season, x, y); }

inline void planet::setseasonalvwind(int season, int x, int y, int amount) { climate_.setseasonalvwind(season, x, y, amount); }

inline int planet::seasonalupperheight(int season, int x, int y) const { return climate_.seasonalupperheight(season, x, y); }

inline void planet::setseasonalupperheight(int season, int x, int y, int amount) { climate_.setseasonalupperheight(season, x, y, amount); }

inline int planet::seasonalupperuwind(int season, int x, int y) const { return climate_.seasonalupperuwind(season, x, y); }

inline void planet::setseasonalupperuwind(int season, int x, int y, int amount) { climate_.setseasonalupperuwind(season, x, y, amount); }

inline int planet::seasonaluppervwind(int season, int x, int y) const { return climate_.seasonaluppervwind(season, x, y); }

inline void planet::setseasonaluppervwind(int season, int x, int y, int amount) { climate_.setseasonaluppervwind(season, x, y, amount); }

inline int planet::seasonalverticalvelocity(int season, int x, int y) const { return climate_.seasonalverticalvelocity(season, x, y); }

inline void planet::setseasonalverticalvelocity(int season, int x, int y, int amount) { climate_.setseasonalverticalvelocity(season, x, y, amount); }

inline int planet::seasonalcurrentu(int season, int x, int y) const { return climate_.seasonalcurrentu(season, x, y); }

inline void planet::setseasonalcurrentu(int season, int x, int y, int amount) { climate_.setseasonalcurrentu(season, x, y, amount); }

inline int planet::seasonalcurrentv(int season, int x, int y) const { return climate_.seasonalcurrentv(season, x, y); }

inline void planet::setseasonalcurrentv(int season, int x, int y, int amount) { climate_.setseasonalcurrentv(season, x, y, amount); }

inline int planet::seasonalsst(int season, int x, int y) const { return climate_.seasonalsst(season, x, y); }

inline void planet::setseasonalsst(int season, int x, int y, int amount) { climate_.setseasonalsst(season, x, y, amount); }

inline float planet::seasonalsstfloat(int season, int x, int y) const { return climate_.seasonalsstfloat(season, x, y); }

inline int planet::seasonalevaporation(int season, int x, int y) const { return climate_.seasonalevaporation(season, x, y); }

inline void planet::setseasonalevaporation(int season, int x, int y, int amount) { climate_.setseasonalevaporation(season, x, y, amount); }

inline int planet::seasonalmaritimeinfluence(int season, int x, int y) const { return climate_.seasonalmaritimeinfluence(season, x, y); }

inline void planet::setseasonalmaritimeinfluence(int season, int x, int y, int amount) { climate_.setseasonalmaritimeinfluence(season, x, y, amount); }

inline int planet::seasonalmaritimethermalanomaly(int season, int x, int y) const { return climate_.seasonalmaritimethermalanomaly(season, x, y); }

inline void planet::setseasonalmaritimethermalanomaly(int season, int x, int y, int amount) { climate_.setseasonalmaritimethermalanomaly(season, x, y, amount); }

inline int planet::seasonalmaritimefetch(int season, int x, int y) const { return climate_.seasonalmaritimefetch(season, x, y); }

inline void planet::setseasonalmaritimefetch(int season, int x, int y, int amount) { climate_.setseasonalmaritimefetch(season, x, y, amount); }

inline int planet::seasonalmoisture(int season, int x, int y) const { return climate_.seasonalmoisture(season, x, y); }

inline void planet::setseasonalmoisture(int season, int x, int y, int amount) { climate_.setseasonalmoisture(season, x, y, amount); }

inline int planet::seasonalconvergence(int season, int x, int y) const { return climate_.seasonalconvergence(season, x, y); }

inline void planet::setseasonalconvergence(int season, int x, int y, int amount) { climate_.setseasonalconvergence(season, x, y, amount); }

inline int planet::seasonaluplift(int season, int x, int y) const { return climate_.seasonaluplift(season, x, y); }

inline void planet::setseasonaluplift(int season, int x, int y, int amount) { climate_.setseasonaluplift(season, x, y, amount); }

inline int planet::seasonalsubsidence(int season, int x, int y) const { return climate_.seasonalsubsidence(season, x, y); }

inline void planet::setseasonalsubsidence(int season, int x, int y, int amount) { climate_.setseasonalsubsidence(season, x, y, amount); }

inline int planet::janmountainrain(int x, int y) const { return climate_.janmountainrain(x, y); }

inline void planet::setjanmountainrain(int x, int y, int amount) { climate_.setjanmountainrain(x, y, amount); }

inline int planet::julmountainrain(int x, int y) const { return climate_.julmountainrain(x, y); }

inline void planet::setjulmountainrain(int x, int y, int amount) { climate_.setjulmountainrain(x, y, amount); }

inline int planet::wintermountainrain(int x, int y) const { return climate_.wintermountainrain(x, y); }

inline void planet::setwintermountainrain(int x, int y, int amount) { climate_.setwintermountainrain(x, y, amount); }

inline int planet::summermountainrain(int x, int y) const { return climate_.summermountainrain(x, y); }

inline void planet::setsummermountainrain(int x, int y, int amount) { climate_.setsummermountainrain(x, y, amount); }

inline int planet::janmountainraindir(int x, int y) const { return climate_.janmountainraindir(x, y); }

inline void planet::setjanmountainraindir(int x, int y, int amount) { climate_.setjanmountainraindir(x, y, amount); }

inline int planet::julmountainraindir(int x, int y) const { return climate_.julmountainraindir(x, y); }

inline void planet::setjulmountainraindir(int x, int y, int amount) { climate_.setjulmountainraindir(x, y, amount); }

inline int planet::wintermountainraindir(int x, int y) const { return climate_.wintermountainraindir(x, y); }

inline int planet::summermountainraindir(int x, int y) const { return climate_.summermountainraindir(x, y); }

inline void planet::setwintermountainraindir(int x, int y, int amount) { climate_.setwintermountainraindir(x, y, amount); }

inline void planet::setsummermountainraindir(int x, int y, int amount) { climate_.setsummermountainraindir(x, y, amount); }

inline int planet::climate(int x, int y) const { return climate_.climate(x, y); }

inline void planet::setclimate(int x, int y, int amount) { climate_.setclimate(x, y, amount); }

inline int planet::biome(int x, int y) const { return climate_.biome(x, y); }

inline void planet::setbiome(int x, int y, int amount) { climate_.setbiome(x, y, amount); }

inline int planet::seaice(int x, int y) const { return climate_.seaice(x, y); }

inline void planet::setseaice(int x, int y, int amount) { climate_.setseaice(x, y, amount); }

inline int planet::riverdir(int x, int y) const
{
    UW_FIELD_READ(river_direction);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (short)hydrosphere_.rivermapdir[x][y];
}

inline void planet::setriverdir(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.rivermapdir[x][y] = (short)amount;
}

inline int planet::riverjan(int x, int y) const
{
    UW_FIELD_READ(river_flow_january);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return hydrosphere_.rivermapjan[x][y];
}

inline void planet::setriverjan(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.rivermapjan[x][y] = amount;
}

inline int planet::riverjul(int x, int y) const
{
    UW_FIELD_READ(river_flow_july);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return hydrosphere_.rivermapjul[x][y];
}

inline void planet::setriverjul(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.rivermapjul[x][y] = amount;
}

inline int planet::riveraveflow(int x, int y) const
{
    UW_FIELD_READ(river_flow_january);
    UW_FIELD_READ(river_flow_july);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (hydrosphere_.rivermapjan[x][y] + hydrosphere_.rivermapjul[x][y]) / 2;
}

inline int planet::wind(int x, int y) const { return climate_.wind(x, y); }

inline void planet::setwind(int x, int y, int amount) { climate_.setwind(x, y, amount); }

inline int planet::lakesurface(int x, int y) const
{
    UW_FIELD_READ(lake_surface);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(hydrosphere_.lakemap[x][y]);
}

inline float planet::lakesurfacefloat(int x, int y) const
{
    UW_FIELD_READ(lake_surface);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return hydrosphere_.lakemap[x][y];
}

inline void planet::setlakesurface(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.lakemap[x][y] = amount;
}

inline float planet::roughness(int x, int y) const
{
    UW_FIELD_READ(roughness);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.roughnessmap[x][y];
}

inline void planet::setroughness(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.roughnessmap[x][y] = amount;
}

inline int planet::mountainridge(int x, int y) const
{
    UW_FIELD_READ(mountain_ridges);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)geosphere_.mountainridges[x][y];
}

inline void planet::setmountainridge(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.mountainridges[x][y] = (short)amount;
}

inline int planet::mountainheight(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.mountainheights[x][y]);
}

inline float planet::mountainheightfloat(int x, int y) const
{
    UW_FIELD_READ(mountain_height);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.mountainheights[x][y];
}

inline void planet::setmountainheight(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.mountainheights[x][y] = amount;
}

inline int planet::craterrim(int x, int y) const
{
    UW_FIELD_READ(crater_rims);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.craterrims[x][y]);
}

inline float planet::craterrimfloat(int x, int y) const
{
    UW_FIELD_READ(crater_rims);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.craterrims[x][y];
}

inline void planet::setcraterrim(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.craterrims[x][y] = amount;
}

inline int planet::cratercentre(int x, int y) const
{
    UW_FIELD_READ(crater_centres);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(geosphere_.cratercentres[x][y]);
}

inline float planet::cratercentrefloat(int x, int y) const
{
    UW_FIELD_READ(crater_centres);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.cratercentres[x][y];
}

inline void planet::setcratercentre(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.cratercentres[x][y] = amount;
}

inline int planet::craterx(int n) const
{
    UW_TABLE_READ(craters);
    if (n<0 || n>=MAXCRATERS)
        return 0;

    return (int)geosphere_.cratercentreslist[n].x;
}

inline void planet::setcraterx(int n, int amount)
{
    if (n < 0 || n >= MAXCRATERS)
        return;

    geosphere_.cratercentreslist[n].x=short(amount);
}

inline int planet::cratery(int n) const
{
    UW_TABLE_READ(craters);
    if (n < 0 || n >= MAXCRATERS)
        return 0;

    return (int)geosphere_.cratercentreslist[n].y;
}

inline void planet::setcratery(int n, int amount)
{
    if (n < 0 || n >= MAXCRATERS)
        return;

    geosphere_.cratercentreslist[n].y = short(amount);
}

inline int planet::craterelev(int n) const
{
    UW_TABLE_READ(craters);
    if (n < 0 || n >= MAXCRATERS)
        return 0;

    return (int)geosphere_.cratercentreslist[n].w;
}

inline void planet::setcraterelev(int n, int amount)
{
    if (n < 0 || n >= MAXCRATERS)
        return;

    geosphere_.cratercentreslist[n].w = short(amount);
}

    inline int planet::craterradius(int n) const
{
    UW_TABLE_READ(craters);
    if (n<0 || n>=MAXCRATERS)
        return 0;

    return (int)geosphere_.cratercentreslist[n].z;
}

inline void planet::setcraterradius(int n, int amount)
{
    if (n < 0 || n >= MAXCRATERS)
        return;

    geosphere_.cratercentreslist[n].z=short(amount);
}

inline int planet::tide(int x, int y) const
{
    UW_FIELD_READ(tide);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)hydrosphere_.tidalmap[x][y];
}

inline void planet::settide(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.tidalmap[x][y] = (short)amount;
}

inline int planet::riftlakesurface(int x, int y) const
{
    UW_FIELD_READ(rift_lake_surface);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(hydrosphere_.riftlakemapsurface[x][y]);
}

inline float planet::riftlakesurfacefloat(int x, int y) const
{
    UW_FIELD_READ(rift_lake_surface);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return hydrosphere_.riftlakemapsurface[x][y];
}

inline void planet::setriftlakesurface(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.riftlakemapsurface[x][y] = amount;
}

inline int planet::riftlakebed(int x, int y) const
{
    UW_FIELD_READ(rift_lake_bed);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return wholemetres(hydrosphere_.riftlakemapbed[x][y]);
}

inline float planet::riftlakebedfloat(int x, int y) const
{
    UW_FIELD_READ(rift_lake_bed);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return hydrosphere_.riftlakemapbed[x][y];
}

inline void planet::setriftlakebed(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.riftlakemapbed[x][y] = amount;
}

inline int planet::special(int x, int y) const
{
    UW_FIELD_READ(special_features);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return hydrosphere_.specials[x][y];
}

inline void planet::setspecial(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.specials[x][y] = (short)amount;
}

inline GeologicRegime planet::geologicregime(int x, int y) const
{
    UW_FIELD_READ(geologic_regime);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return GeologicRegime::stable;

    return static_cast<GeologicRegime>(geosphere_.geologicregimemap[x][y]);
}

inline void planet::setgeologicregime(int x, int y, GeologicRegime amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.geologicregimemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::tectonicconvergence(int x, int y) const
{
    UW_FIELD_READ(tectonic_convergence);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.tectonicconvergencemap[x][y]);
}

inline void planet::settectonicconvergence(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.tectonicconvergencemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::tectonicdivergence(int x, int y) const
{
    UW_FIELD_READ(tectonic_divergence);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.tectonicdivergencemap[x][y]);
}

inline void planet::settectonicdivergence(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.tectonicdivergencemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::tectonicshear(int x, int y) const
{
    UW_FIELD_READ(tectonic_shear);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.tectonicshearmap[x][y]);
}

inline void planet::settectonicshear(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.tectonicshearmap[x][y] = static_cast<std::uint8_t>(amount);
}

inline float planet::tectoniccrustagemyr(int x, int y) const
{
    UW_FIELD_READ(crust_age);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectoniccrustagemyrmap[x][y];
}

inline void planet::settectoniccrustagemyr(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectoniccrustagemyrmap[x][y] = (std::max)(0.0f, amount);
}

inline float planet::tectoniccrustthickness(int x, int y) const
{
    UW_FIELD_READ(crust_thickness);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectoniccrustthicknessmap[x][y];
}

inline void planet::settectoniccrustthickness(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectoniccrustthicknessmap[x][y] = (std::max)(0.0f, amount);
}

inline CrustClass planet::tectoniccrustclass(int x, int y) const
{
    UW_FIELD_READ(crust_class);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return CrustClass::none;

    return static_cast<CrustClass>(geosphere_.tectoniccrustclassmap[x][y]);
}

inline void planet::settectoniccrustclass(int x, int y, CrustClass amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectoniccrustclassmap[x][y] = static_cast<std::uint8_t>(amount);
}

inline float planet::tectonicuplifttendency(int x, int y) const
{
    UW_FIELD_READ(uplift_tendency);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicuplifttendencymap[x][y];
}

inline void planet::settectonicuplifttendency(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicuplifttendencymap[x][y] = std::clamp(amount, 0.0f, 1.0f);
}

inline float planet::tectonicsubsidencetendency(int x, int y) const
{
    UW_FIELD_READ(subsidence_tendency);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicsubsidencetendencymap[x][y];
}

inline void planet::settectonicsubsidencetendency(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicsubsidencetendencymap[x][y] = std::clamp(amount, 0.0f, 1.0f);
}

inline float planet::tectonicaccumulatedstrain(int x, int y) const
{
    UW_FIELD_READ(accumulated_strain);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicaccumulatedstrainmap[x][y];
}

inline void planet::settectonicaccumulatedstrain(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicaccumulatedstrainmap[x][y] = std::clamp(amount, 0.0f, 1.0f);
}

inline BoundaryType planet::tectonicboundarytype(int x, int y) const
{
    UW_FIELD_READ(boundary_type);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return BoundaryType::none;

    return static_cast<BoundaryType>(geosphere_.tectonicboundarytypemap[x][y]);
}

inline void planet::settectonicboundarytype(int x, int y, BoundaryType amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicboundarytypemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::tectonicboundarydistance(int x, int y) const
{
    UW_FIELD_READ(boundary_distance);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.tectonicboundarydistancemap[x][y]);
}

inline void planet::settectonicboundarydistance(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicboundarydistancemap[x][y] = static_cast<unsigned short>(std::clamp(amount, 0, 65535));
}

inline int planet::tectonicboundarysegmentid(int x, int y) const
{
    UW_FIELD_READ(boundary_segment);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.tectonicboundarysegmentidmap[x][y];
}

inline void planet::settectonicboundarysegmentid(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicboundarysegmentidmap[x][y] = (std::max)(0, amount);
}

inline int planet::tectonicnearestboundaryid(int x, int y) const
{
    UW_FIELD_READ(nearest_boundary);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.tectonicnearestboundaryidmap[x][y];
}

inline void planet::settectonicnearestboundaryid(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicnearestboundaryidmap[x][y] = (std::max)(0, amount);
}

inline float planet::tectonicboundaryhistory(int x, int y) const
{
    UW_FIELD_READ(boundary_history);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicboundaryhistorymap[x][y];
}

inline void planet::settectonicboundaryhistory(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicboundaryhistorymap[x][y] = std::clamp(amount, 0.0f, 1.0f);
}

inline int planet::tectonicdeformingregionid(int x, int y) const
{
    UW_FIELD_READ(deforming_region);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.tectonicdeformingregionidmap[x][y];
}

inline void planet::settectonicdeformingregionid(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicdeformingregionidmap[x][y] = (std::max)(0, amount);
}

inline DeformingRegionType planet::tectonicdeformingregiontype(int x, int y) const
{
    UW_FIELD_READ(deforming_region_type);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return DeformingRegionType::none;

    return static_cast<DeformingRegionType>(geosphere_.tectonicdeformingregiontypemap[x][y]);
}

inline void planet::settectonicdeformingregiontype(int x, int y, DeformingRegionType amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicdeformingregiontypemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline float planet::tectonicdeformationrate(int x, int y) const
{
    UW_FIELD_READ(deformation_rate);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicdeformationratemap[x][y];
}

inline void planet::settectonicdeformationrate(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicdeformationratemap[x][y] = std::clamp(amount, 0.0f, 1.0f);
}

inline float planet::tectonicdeformationvelocityx(int x, int y) const
{
    UW_FIELD_READ(deformation_velocity_x);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicdeformationvelocityxmap[x][y];
}

inline void planet::settectonicdeformationvelocityx(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicdeformationvelocityxmap[x][y] = amount;
}

inline float planet::tectonicdeformationvelocityy(int x, int y) const
{
    UW_FIELD_READ(deformation_velocity_y);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0.0f;

    return geosphere_.tectonicdeformationvelocityymap[x][y];
}

inline void planet::settectonicdeformationvelocityy(int x, int y, float amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicdeformationvelocityymap[x][y] = amount;
}

inline BasinClass planet::basinclass(int x, int y) const
{
    UW_FIELD_READ(basin_class);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return BasinClass::none;

    return static_cast<BasinClass>(hydrosphere_.basinclassmap[x][y]);
}

inline void planet::setbasinclass(int x, int y, BasinClass amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.basinclassmap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::erosionpotential(int x, int y) const
{
    UW_FIELD_READ(erosion_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.erosionpotentialmap[x][y]);
}

inline void planet::seterosionpotential(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.erosionpotentialmap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::depositionpotential(int x, int y) const
{
    UW_FIELD_READ(deposition_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.depositionpotentialmap[x][y]);
}

inline void planet::setdepositionpotential(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.depositionpotentialmap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::floodplainfertility(int x, int y) const
{
    UW_FIELD_READ(floodplain_fertility);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.floodplainfertilitymap[x][y]);
}

inline void planet::setfloodplainfertility(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.floodplainfertilitymap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::metalorereserve(int x, int y) const
{
    UW_FIELD_READ(metal_ore_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.metalorereservemap[x][y]);
}

inline void planet::setmetalorereserve(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.metalorereservemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::placerreserve(int x, int y) const
{
    UW_FIELD_READ(placer_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.placerreservemap[x][y]);
}

inline void planet::setplacerreserve(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.placerreservemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::evaporitereserve(int x, int y) const
{
    UW_FIELD_READ(evaporite_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.evaporitereservemap[x][y]);
}

inline void planet::setevaporitereserve(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.evaporitereservemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::volcanicreserve(int x, int y) const
{
    UW_FIELD_READ(volcanic_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(geosphere_.volcanicreservemap[x][y]);
}

inline void planet::setvolcanicreserve(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    geosphere_.volcanicreservemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::fisheryreserve(int x, int y) const
{
    UW_FIELD_READ(fishery_potential);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return static_cast<int>(biosphere_.fisheryreservemap[x][y]);
}

inline void planet::setfisheryreserve(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    amount = std::clamp(amount, 0, 100);
    biosphere_.fisheryreservemap[x][y] = static_cast<std::uint8_t>(amount);
}

inline int planet::settlementsuitability(int x, int y) const
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return 0;

    return static_cast<int>(anthroposphere_.settlement_suitability[x][y]);
}

inline void planet::setsettlementsuitability(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.settlement_suitability[x][y] = static_cast<std::uint8_t>(std::clamp(amount, 0, 100));
}

inline int planet::infrastructure(int x, int y) const
{
    UW_FIELD_READ(infrastructure);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return 0;

    return static_cast<int>(anthroposphere_.social_infrastructure[x][y]);
}

inline void planet::setinfrastructure(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.social_infrastructure[x][y] = static_cast<std::uint8_t>(std::clamp(amount, 0, 100));
}

inline int planet::agriculturalcapacity(int x, int y) const
{
    UW_FIELD_READ(agricultural_capacity);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return 0;

    return static_cast<int>(anthroposphere_.agricultural_capacity[x][y]);
}

inline void planet::setagriculturalcapacity(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.agricultural_capacity[x][y] = static_cast<std::uint8_t>(std::clamp(amount, 0, 100));
}

inline int planet::routetraffic(int x, int y) const
{
    UW_FIELD_READ(route_traffic);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return 0;

    return static_cast<int>(anthroposphere_.route_traffic[x][y]);
}

inline void planet::setroutetraffic(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.route_traffic[x][y] = static_cast<std::uint8_t>(std::clamp(amount, 0, 100));
}

inline int planet::riveraccess(int x, int y) const
{
    UW_FIELD_READ(river_access);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return 0;

    return static_cast<int>(anthroposphere_.river_access[x][y]);
}

inline void planet::setriveraccess(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.river_access[x][y] = static_cast<std::uint8_t>(std::clamp(amount, 0, 100));
}

inline int planet::harborscore(int x, int y) const
{
    UW_FIELD_READ(harbor_score);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return 0;

    return static_cast<int>(anthroposphere_.harbor_score[x][y]);
}

inline void planet::setharborscore(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.harbor_score[x][y] = static_cast<std::uint8_t>(std::clamp(amount, 0, 100));
}

inline int planet::ownersettlementid(int x, int y) const
{
    UW_FIELD_READ(owner_settlement);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return -1;

    return anthroposphere_.owner_settlement_id[x][y];
}

inline void planet::setownersettlementid(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.owner_settlement_id[x][y] = amount;
}

inline int planet::ownerpolityid(int x, int y) const
{
    UW_FIELD_READ(owner_polity);
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return -1;

    return anthroposphere_.owner_polity_id[x][y];
}

inline void planet::setownerpolityid(int x, int y, int amount)
{
    if (y < 0 || y > body_.itsheight || x < 0 || x > body_.itswidth)
        return;

    anthroposphere_.owner_polity_id[x][y] = amount;
}

inline int planet::deltadir(int x, int y) const
{
    UW_FIELD_READ(delta_direction);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)hydrosphere_.deltamapdir[x][y];
}

inline void planet::setdeltadir(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.deltamapdir[x][y] = (short)amount;
}

inline int planet::deltajan(int x, int y) const
{
    UW_FIELD_READ(delta_flow_january);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return hydrosphere_.deltamapjan[x][y];
}

inline void planet::setdeltajan(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.deltamapjan[x][y] = amount;
}

inline int planet::deltajul(int x, int y) const
{
    UW_FIELD_READ(delta_flow_july);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return hydrosphere_.deltamapjul[x][y];
}

inline void planet::setdeltajul(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.deltamapjul[x][y] = amount;
}

inline bool planet::lakestart(int x, int y) const
{
    UW_FIELD_READ(lake_start);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return hydrosphere_.lakestartmap[x][y];
}

inline void planet::setlakestart(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    hydrosphere_.lakestartmap[x][y] = amount;
}

inline bool planet::island(int x, int y) const
{
    UW_FIELD_READ(island);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.islandmap[x][y];
}

inline void planet::setisland(int x, int y, bool amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.islandmap[x][y] = amount;
}

inline bool planet::mountainisland(int x, int y) const
{
    UW_FIELD_READ(mountain_island);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return geosphere_.mountainislandmap[x][y];
}

inline void planet::setmountainisland(int x, int y, bool amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.mountainislandmap[x][y] = amount;
}

inline bool planet::noshade(int x, int y) const
{
    UW_FIELD_READ(no_shade);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return presentation_.noshademap[x][y];
}

inline void planet::setnoshade(int x, int y, bool amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    presentation_.noshademap[x][y] = amount;
}

inline int planet::noise(int x, int y) const
{
    UW_FIELD_READ(detail_noise);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return (int)geosphere_.noisemap[x][y];
}

inline void planet::setnoise(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.noisemap[x][y] = (short)amount;
}

inline int planet::horse(int x, int y) const { return climate_.horse(x, y); }

inline void planet::sethorse(int x, int y, int amount) { climate_.sethorse(x, y, amount); }

inline int planet::test(int x, int y) const
{
    UW_FIELD_READ(debug_test);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    return presentation_.testmap[x][y];
}

inline void planet::settest(int x, int y, int amount)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    presentation_.testmap[x][y] = amount;
}

inline const std::vector<Settlement>& planet::settlements() const { return anthroposphere_.settlementlist; }
inline std::vector<Settlement>& planet::settlements() { return anthroposphere_.settlementlist; }
inline const std::vector<Polity>& planet::polities() const { UW_TABLE_READ(polities); return anthroposphere_.politylist; }
inline std::vector<Polity>& planet::polities() { UW_TABLE_READ(polities); return anthroposphere_.politylist; }
inline const std::vector<RouteEdge>& planet::routeedges() const { UW_TABLE_READ(route_edges); return anthroposphere_.routeedgelist; }
inline std::vector<RouteEdge>& planet::routeedges() { UW_TABLE_READ(route_edges); return anthroposphere_.routeedgelist; }
inline const std::vector<HistoryEvent>& planet::historyevents() const { UW_TABLE_READ(history_events); return anthroposphere_.historyeventlist; }
inline std::vector<HistoryEvent>& planet::historyevents() { UW_TABLE_READ(history_events); return anthroposphere_.historyeventlist; }

inline const std::vector<mineralsystems::Deposit>& planet::mineraldeposits() const { UW_TABLE_READ(mineral_deposits); return geosphere_.mineraldepositlist; }
inline std::vector<mineralsystems::Deposit>& planet::mineraldeposits() { UW_TABLE_READ(mineral_deposits); return geosphere_.mineraldepositlist; }
inline void planet::setmineraldeposits(std::vector<mineralsystems::Deposit> deposits) { geosphere_.mineraldepositlist = std::move(deposits); }
inline const mineralsystems::MineralBudget& planet::mineralbudget() const { UW_TABLE_READ(mineral_budget); return geosphere_.mineralbudgetstate; }
inline void planet::setmineralbudget(const mineralsystems::MineralBudget& budget) { geosphere_.mineralbudgetstate = budget; }

inline const std::vector<TectonicBoundarySegment>& planet::tectonicboundarysegments() const { UW_TABLE_READ(tectonic_boundary_segments); return geosphere_.tectonicboundarysegmentlist; }
inline std::vector<TectonicBoundarySegment>& planet::tectonicboundarysegments() { UW_TABLE_READ(tectonic_boundary_segments); return geosphere_.tectonicboundarysegmentlist; }
inline int planet::tectonicboundarysegmentcount() const { UW_TABLE_READ(tectonic_boundary_segments); return static_cast<int>(geosphere_.tectonicboundarysegmentlist.size()); }
inline void planet::settectonicsnapshotcell(int x, int y, float height, int plate, int crustagesteps,
    float continentalcrust, float oceaniccrust, float sediment, float riftdamage, float plumeinfluence)
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return;

    geosphere_.tectonicheightmap[x][y] = height;
    geosphere_.tectonicplateidmap[x][y] = plate;
    geosphere_.tectoniccrustagestepsmap[x][y] = crustagesteps;
    geosphere_.tectoniccontinentalcrustmap[x][y] = continentalcrust;
    geosphere_.tectonicoceaniccrustmap[x][y] = oceaniccrust;
    geosphere_.tectonicsedimentmap[x][y] = sediment;
    geosphere_.tectonicriftdamagemap[x][y] = riftdamage;
    geosphere_.tectonicplumeinfluencemap[x][y] = plumeinfluence;
}

inline void planet::settectonicsnapshotrecord(std::string record) { geosphere_.tectonicsnapshotrecord = std::move(record); }

inline void planet::settectonicboundarysegments(std::vector<TectonicBoundarySegment> segments) { geosphere_.tectonicboundarysegmentlist = std::move(segments); }

template<typename SegmentT>
inline void planet::settectonicboundarysegments(const std::vector<SegmentT>& segments)
{
    geosphere_.tectonicboundarysegmentlist.resize(segments.size());

    for (size_t i = 0; i < segments.size(); i++)
    {
        const SegmentT& src = segments[i];
        TectonicBoundarySegment& dst = geosphere_.tectonicboundarysegmentlist[i];
        dst.id = static_cast<int>(src.id);
        dst.leftPlateId = static_cast<int>(src.left_plate_id);
        dst.rightPlateId = static_cast<int>(src.right_plate_id);
        dst.cellCount = static_cast<int>(src.cell_count);
        dst.persistenceSteps = static_cast<int>(src.persistence_steps);
        dst.centroidX = src.centroid_x;
        dst.centroidY = src.centroid_y;
        dst.lengthCells = src.length_cells;
        dst.averageNormalMotion = src.average_normal_motion;
        dst.averageShearMotion = src.average_shear_motion;
        dst.averageConvergenceScore = static_cast<int>(src.average_convergence_score);
        dst.averageDivergenceScore = static_cast<int>(src.average_divergence_score);
        dst.averageShearScore = static_cast<int>(src.average_shear_score);
        dst.boundaryType = static_cast<BoundaryType>(static_cast<std::uint8_t>(src.boundary_type));
        dst.geologicRegime = static_cast<GeologicRegime>(static_cast<std::uint8_t>(src.geologic_regime));
        dst.ageMyr = static_cast<double>(src.age_myr);
    }
}

inline const std::vector<TectonicDeformingRegion>& planet::tectonicdeformingregions() const { UW_TABLE_READ(tectonic_deforming_regions); return geosphere_.tectonicdeformingregionlist; }
inline std::vector<TectonicDeformingRegion>& planet::tectonicdeformingregions() { UW_TABLE_READ(tectonic_deforming_regions); return geosphere_.tectonicdeformingregionlist; }
inline int planet::tectonicdeformingregioncount() const { UW_TABLE_READ(tectonic_deforming_regions); return static_cast<int>(geosphere_.tectonicdeformingregionlist.size()); }
inline void planet::settectonicdeformingregions(std::vector<TectonicDeformingRegion> regions) { geosphere_.tectonicdeformingregionlist = std::move(regions); }

template<typename RegionT>
inline void planet::settectonicdeformingregions(const std::vector<RegionT>& regions)
{
    geosphere_.tectonicdeformingregionlist.resize(regions.size());

    for (size_t i = 0; i < regions.size(); i++)
    {
        const RegionT& src = regions[i];
        TectonicDeformingRegion& dst = geosphere_.tectonicdeformingregionlist[i];
        dst.id = static_cast<int>(src.id);
        dst.boundarySegmentId = static_cast<int>(src.boundary_segment_id);
        dst.primaryPlateId = static_cast<int>(src.primary_plate_id);
        dst.secondaryPlateId = static_cast<int>(src.secondary_plate_id);
        dst.cellCount = static_cast<int>(src.cell_count);
        dst.persistenceSteps = static_cast<int>(src.persistence_steps);
        dst.centroidX = src.centroid_x;
        dst.centroidY = src.centroid_y;
        dst.averageDeformationRate = src.average_deformation_rate;
        dst.averageInterpolatedVelocityX = src.average_interpolated_velocity_x;
        dst.averageInterpolatedVelocityY = src.average_interpolated_velocity_y;
        dst.averageNormalMotion = src.average_normal_motion;
        dst.averageShearMotion = src.average_shear_motion;
        dst.type = static_cast<DeformingRegionType>(static_cast<std::uint8_t>(src.type));
        dst.ageMyr = static_cast<double>(src.age_myr);
    }
}

inline const TectonicBoundarySegment* planet::findtectonicboundarysegment(int id) const
{
    UW_TABLE_READ(tectonic_boundary_segments);
    if (id <= 0)
        return nullptr;

    for (const TectonicBoundarySegment& segment : geosphere_.tectonicboundarysegmentlist)
    {
        if (segment.id == id)
            return &segment;
    }

    return nullptr;
}

inline const TectonicDeformingRegion* planet::findtectonicdeformingregion(int id) const
{
    UW_TABLE_READ(tectonic_deforming_regions);
    if (id <= 0)
        return nullptr;

    for (const TectonicDeformingRegion& region : geosphere_.tectonicdeformingregionlist)
    {
        if (region.id == id)
            return &region;
    }

    return nullptr;
}

inline void planet::clearsocialstate()
{
    anthroposphere_.settlementlist.clear();
    anthroposphere_.politylist.clear();
    anthroposphere_.routeedgelist.clear();
    anthroposphere_.historyeventlist.clear();
}

inline bool planet::sea(int x, int y) const
{
    UW_FIELD_READ(base_elevation);
    UW_FIELD_READ(sea_level);
    UW_FIELD_READ(lake_surface);
    UW_FIELD_READ(volcano);
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    if (wholemetres(geosphere_.mapnom[x][y]) <= geosphere_.itssealevel && hydrosphere_.lakemap[x][y] == 0)
        return 1;

    if (geosphere_.volcanomap[x][y] == 0)
        return 0;

    const float thisvolcano = std::fabs(geosphere_.volcanomap[x][y]);

    if (wholemetres(geosphere_.mapnom[x][y] + thisvolcano) <= geosphere_.itssealevel)
        return 1;
    else
        return 0;
}

inline int planet::mintemp(int x, int y) const { return climate_.mintemp(x, y); }

inline void planet::setmintemp(int x, int y, int amount) { climate_.setmintemp(x, y, amount); }

inline float planet::mintempfloat(int x, int y) const { return climate_.mintempfloat(x, y); }
inline float planet::maxtempfloat(int x, int y) const { return climate_.maxtempfloat(x, y); }
inline float planet::avetempfloat(int x, int y) const { return climate_.avetempfloat(x, y); }
inline float planet::summerrainfloat(int x, int y) const { return climate_.summerrainfloat(x, y); }
inline float planet::winterrainfloat(int x, int y) const { return climate_.winterrainfloat(x, y); }
inline float planet::jantempfloat(int x, int y) const { return climate_.jantempfloat(x, y); }
inline float planet::jultempfloat(int x, int y) const { return climate_.jultempfloat(x, y); }
inline float planet::janrainfloat(int x, int y) const { return climate_.janrainfloat(x, y); }
inline float planet::julrainfloat(int x, int y) const { return climate_.julrainfloat(x, y); }
inline void planet::setsummerrainfloat(int x, int y, float amount) { climate_.setsummerrainfloat(x, y, amount); }
inline void planet::setwinterrainfloat(int x, int y, float amount) { climate_.setwinterrainfloat(x, y, amount); }
inline void planet::setmintempfloat(int x, int y, float amount) { climate_.setmintempfloat(x, y, amount); }
inline void planet::setmaxtempfloat(int x, int y, float amount) { climate_.setmaxtempfloat(x, y, amount); }
inline void planet::setjantempfloat(int x, int y, float amount) { climate_.setjantempfloat(x, y, amount); }
inline void planet::setjultempfloat(int x, int y, float amount) { climate_.setjultempfloat(x, y, amount); }
inline void planet::setjanrainfloat(int x, int y, float amount) { climate_.setjanrainfloat(x, y, amount); }
inline void planet::setjulrainfloat(int x, int y, float amount) { climate_.setjulrainfloat(x, y, amount); }

inline int planet::maxtemp(int x, int y) const { return climate_.maxtemp(x, y); }

inline void planet::setmaxtemp(int x, int y, int amount) { climate_.setmaxtemp(x, y, amount); }

inline int planet::winterrain(int x, int y) const { return climate_.winterrain(x, y); }

inline void planet::setwinterrain(int x, int y, int amount) { climate_.setwinterrain(x, y, amount); }

inline int planet::summerrain(int x, int y) const { return climate_.summerrain(x, y); }

inline void planet::setsummerrain(int x, int y, int amount) { climate_.setsummerrain(x, y, amount); }

inline int planet::truelake(int x, int y) const
{
    if (y<0 || y>body_.itsheight || x<0 || x>body_.itswidth)
        return 0;

    if (lakesurface(x, y) != 0 && special(x, y) < 110)
        return (1);
    else
        return(0);
}


#endif /* planet_hpp */
