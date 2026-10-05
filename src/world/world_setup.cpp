// Initial state of a world and a region, default map appearance, land and sea totals.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "world/world_setup.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <queue>

#include "core/classes.hpp"
#include "world/world_parameters.hpp"
#include "world/planet.hpp"
#include "world/functions.hpp"

using namespace std;

// This function plugs the global variables into the world.

void initialiseworld(planet& world)
{
    world.clear();

    int saveversion = 18;        // Saves of this version (and version 17) can be loaded.
    int settingssaveversion = 5; // As above, but for settings files.
    int size = tuning::world::defaults::size;
    int width = tuning::world::defaults::width;
    int height = tuning::world::defaults::height;
    bool rotation = tuning::world::defaults::rotation;
    float tilt = tuning::world::defaults::tilt;
    float eccentricity = tuning::world::defaults::eccentricity;
    short perihelion = tuning::world::defaults::perihelion;
    float gravity = tuning::world::defaults::gravity;
    float lunar = tuning::world::defaults::lunar;
    float tempdecrease = tuning::world::defaults::tempDecrease;
    int northpolaradjust = tuning::world::defaults::northPolarAdjust;
    int southpolaradjust = tuning::world::defaults::southPolarAdjust;
    int averagetemp = tuning::world::defaults::averageTemp;
    float waterpickup = tuning::world::defaults::waterPickup;
    float riverfactor = tuning::world::defaults::riverFactor;
    int riverlandreduce = tuning::world::defaults::riverLandReduce;
    int estuarylimit = tuning::world::defaults::estuaryLimit;
    int glacialtemp = tuning::world::defaults::glacialTemp;
    int glaciertemp = tuning::world::defaults::glacierTemp;
    float mountainreduce = tuning::world::defaults::mountainReduce;
    int climatenumber = tuning::world::defaults::climateCount;
    int maxelevation = tuning::world::defaults::maxElevation;
    int sealevel = tuning::world::defaults::seaLevel;
    int craterno = tuning::world::defaults::craterCount;

    world.setsaveversion(saveversion);
    world.setsettingssaveversion(settingssaveversion);
    world.setsize(size);
    world.setwidth(width);
    world.setrotation(rotation);
    world.settilt(tilt);
    world.seteccentricity(eccentricity);
    world.setperihelion(perihelion);
    world.setgravity(gravity);
    world.setlunar(lunar);
    world.settempdecrease(tempdecrease);
    world.setnorthpolaradjust(northpolaradjust);
    world.setsouthpolaradjust(southpolaradjust);
    world.setaveragetemp(averagetemp);
    world.setwaterpickup(waterpickup);
    world.setriverfactor(riverfactor);
    world.setriverlandreduce(riverlandreduce);
    world.setestuarylimit(estuarylimit);
    world.setglacialtemp(glacialtemp);
    world.setglaciertemp(glaciertemp);
    world.setmountainreduce(mountainreduce);
    world.setclimatenumber(climatenumber);
    world.setmaxelevation(maxelevation);
    world.setsealevel(sealevel);
    world.settectonicplatecount(0);
    world.settectonicsealevelm(0);
    world.setlandtotal(0);
    world.setseatotal(0);
    world.setcraterno(craterno);

    world.setseatotal(0);
    world.setlandtotal(0);
}

void initialisegradientmapappearance(planet& world)
{
    const int elevationmin = 0 - world.sealevel();
    const int elevationmax = max(1, world.maxelevation() - world.sealevel());
    const float sealevelfactor = static_cast<float>(world.sealevel()) / static_cast<float>(max(1, world.maxelevation()));

    const int elevationmid1 = static_cast<int>(roundf(static_cast<float>(world.elevationlow1()) + static_cast<float>(world.elevationhigh1() - world.elevationlow1()) * sealevelfactor));
    const int elevationmid2 = static_cast<int>(roundf(static_cast<float>(world.elevationlow2()) + static_cast<float>(world.elevationhigh2() - world.elevationlow2()) * sealevelfactor));
    const int elevationmid3 = static_cast<int>(roundf(static_cast<float>(world.elevationlow3()) + static_cast<float>(world.elevationhigh3() - world.elevationlow3()) * sealevelfactor));
    const int riverflowmax = max(1000, world.maxriverflow());

    world.setmapgradientstopcount(mapgradientelevation, 3);
    world.setmapgradientdiscrete(mapgradientelevation, false);
    world.setmapgradientposition(mapgradientelevation, 0, elevationmin);
    world.setmapgradientposition(mapgradientelevation, 1, 0);
    world.setmapgradientposition(mapgradientelevation, 2, elevationmax);
    world.setmapgradientcolour(mapgradientelevation, 0, 0, world.elevationlow1());
    world.setmapgradientcolour(mapgradientelevation, 0, 1, world.elevationlow2());
    world.setmapgradientcolour(mapgradientelevation, 0, 2, world.elevationlow3());
    world.setmapgradientcolour(mapgradientelevation, 1, 0, elevationmid1);
    world.setmapgradientcolour(mapgradientelevation, 1, 1, elevationmid2);
    world.setmapgradientcolour(mapgradientelevation, 1, 2, elevationmid3);
    world.setmapgradientcolour(mapgradientelevation, 2, 0, world.elevationhigh1());
    world.setmapgradientcolour(mapgradientelevation, 2, 1, world.elevationhigh2());
    world.setmapgradientcolour(mapgradientelevation, 2, 2, world.elevationhigh3());

    world.setmapgradientstopcount(mapgradienttemperature, 3);
    world.setmapgradientdiscrete(mapgradienttemperature, false);
    world.setmapgradientposition(mapgradienttemperature, 0, -25);
    world.setmapgradientposition(mapgradienttemperature, 1, 0);
    world.setmapgradientposition(mapgradienttemperature, 2, 25);
    world.setmapgradientcolour(mapgradienttemperature, 0, 0, world.temperaturecold1());
    world.setmapgradientcolour(mapgradienttemperature, 0, 1, world.temperaturecold2());
    world.setmapgradientcolour(mapgradienttemperature, 0, 2, world.temperaturecold3());
    world.setmapgradientcolour(mapgradienttemperature, 1, 0, world.temperaturetemperate1());
    world.setmapgradientcolour(mapgradienttemperature, 1, 1, world.temperaturetemperate2());
    world.setmapgradientcolour(mapgradienttemperature, 1, 2, world.temperaturetemperate3());
    world.setmapgradientcolour(mapgradienttemperature, 2, 0, world.temperaturehot1());
    world.setmapgradientcolour(mapgradienttemperature, 2, 1, world.temperaturehot2());
    world.setmapgradientcolour(mapgradienttemperature, 2, 2, world.temperaturehot3());

    world.setmapgradientstopcount(mapgradientprecipitation, 2);
    world.setmapgradientdiscrete(mapgradientprecipitation, false);
    world.setmapgradientposition(mapgradientprecipitation, 0, 0);
    world.setmapgradientposition(mapgradientprecipitation, 1, 1020);
    world.setmapgradientcolour(mapgradientprecipitation, 0, 0, world.precipitationdry1());
    world.setmapgradientcolour(mapgradientprecipitation, 0, 1, world.precipitationdry2());
    world.setmapgradientcolour(mapgradientprecipitation, 0, 2, world.precipitationdry3());
    world.setmapgradientcolour(mapgradientprecipitation, 1, 0, world.precipitationwet1());
    world.setmapgradientcolour(mapgradientprecipitation, 1, 1, world.precipitationwet2());
    world.setmapgradientcolour(mapgradientprecipitation, 1, 2, world.precipitationwet3());

    world.setmapgradientstopcount(mapgradientriverflow, 2);
    world.setmapgradientdiscrete(mapgradientriverflow, false);
    world.setmapgradientposition(mapgradientriverflow, 0, 0);
    world.setmapgradientposition(mapgradientriverflow, 1, riverflowmax);
    world.setmapgradientcolour(mapgradientriverflow, 0, 0, world.rivermapcolour(rivermaplowflow, 0));
    world.setmapgradientcolour(mapgradientriverflow, 0, 1, world.rivermapcolour(rivermaplowflow, 1));
    world.setmapgradientcolour(mapgradientriverflow, 0, 2, world.rivermapcolour(rivermaplowflow, 2));
    world.setmapgradientcolour(mapgradientriverflow, 1, 0, world.rivermapcolour(rivermaphighflow, 0));
    world.setmapgradientcolour(mapgradientriverflow, 1, 1, world.rivermapcolour(rivermaphighflow, 1));
    world.setmapgradientcolour(mapgradientriverflow, 1, 2, world.rivermapcolour(rivermaphighflow, 2));

    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        const int stopcount = world.mapgradientstopcount(gradient);

        for (int stop = stopcount; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            world.setmapgradientposition(gradient, stop, world.mapgradientposition(gradient, stopcount - 1));
            world.setmapgradientcolour(gradient, stop, 0, world.mapgradientcolour(gradient, stopcount - 1, 0));
            world.setmapgradientcolour(gradient, stop, 1, world.mapgradientcolour(gradient, stopcount - 1, 1));
            world.setmapgradientcolour(gradient, stop, 2, world.mapgradientcolour(gradient, stopcount - 1, 2));
        }
    }
}

void setdefaultnonreliefmapappearance(planet& world)
{
    static const int defaultclimatecolours[CLIMATEMAPCOLOURCOUNT][3] =
    {
        { 0, 0, 0 },
        { 0, 0, 254 },
        { 1, 119, 255 },
        { 70, 169, 250 },
        { 70, 169, 250 },
        { 249, 15, 0 },
        { 251, 150, 149 },
        { 245, 163, 1 },
        { 254, 219, 99 },
        { 255, 255, 0 },
        { 198, 199, 1 },
        { 184, 184, 114 },
        { 138, 255, 162 },
        { 86, 199, 112 },
        { 30, 150, 66 },
        { 192, 254, 109 },
        { 76, 255, 93 },
        { 19, 203, 74 },
        { 255, 8, 245 },
        { 204, 3, 192 },
        { 154, 51, 144 },
        { 153, 100, 146 },
        { 172, 178, 249 },
        { 91, 121, 213 },
        { 78, 83, 175 },
        { 54, 3, 130 },
        { 0, 255, 245 },
        { 32, 200, 250 },
        { 0, 126, 125 },
        { 0, 69, 92 },
        { 178, 178, 178 },
        { 104, 104, 104 }
    };

    static const int defaultclimateseacolours[CLIMATEMAPSEACOLOURCOUNT][3] =
    {
        { 13, 49, 109 },
        { 228, 228, 255 },
        { 243, 243, 255 }
    };

    world.setshowmapoutline(1);

    world.setoutline1(0);
    world.setoutline2(0);
    world.setoutline3(0);

    world.setelevationlow1(0);
    world.setelevationlow2(0);
    world.setelevationlow3(0);
    world.setelevationhigh1(255);
    world.setelevationhigh2(255);
    world.setelevationhigh3(255);

    world.settemperaturecold1(0);
    world.settemperaturecold2(0);
    world.settemperaturecold3(250);
    world.settemperaturetemperate1(250);
    world.settemperaturetemperate2(250);
    world.settemperaturetemperate3(250);
    world.settemperaturehot1(250);
    world.settemperaturehot2(100);
    world.settemperaturehot3(0);

    world.setprecipitationdry1(255);
    world.setprecipitationdry2(255);
    world.setprecipitationdry3(255);
    world.setprecipitationwet1(0);
    world.setprecipitationwet2(64);
    world.setprecipitationwet3(255);

    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            world.setclimatemapseacolour(i, j, defaultclimateseacolours[i][j]);
    }

    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            world.setclimatemapcolour(i, j, defaultclimatecolours[i][j]);
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            world.setbiomemapcolour(i, j, defaultbiomemapcolours[i][j]);
    }

    world.setrivermapcolour(rivermapbackground, 0, 255);
    world.setrivermapcolour(rivermapbackground, 1, 255);
    world.setrivermapcolour(rivermapbackground, 2, 255);
    world.setrivermapcolour(rivermaplowflow, 0, 255);
    world.setrivermapcolour(rivermaplowflow, 1, 255);
    world.setrivermapcolour(rivermaplowflow, 2, 255);
    world.setrivermapcolour(rivermaphighflow, 0, world.river1());
    world.setrivermapcolour(rivermaphighflow, 1, world.river2());
    world.setrivermapcolour(rivermaphighflow, 2, world.river3());
    world.setrivermapcolour(rivermaplake, 0, world.lake1());
    world.setrivermapcolour(rivermaplake, 1, world.lake2());
    world.setrivermapcolour(rivermaplake, 2, world.lake3());
    world.setrivermapcolour(rivermapsaltpan, 0, world.saltpan1());
    world.setrivermapcolour(rivermapsaltpan, 1, world.saltpan2());
    world.setrivermapcolour(rivermapsaltpan, 2, world.saltpan3());
    world.setrivermapcolour(rivermapwetlands, 0, world.wetlands1());
    world.setrivermapcolour(rivermapwetlands, 1, world.wetlands2());
    world.setrivermapcolour(rivermapwetlands, 2, world.wetlands3());
    world.setrivermapcolour(rivermapmud, 0, world.mud1());
    world.setrivermapcolour(rivermapmud, 1, world.mud2());
    world.setrivermapcolour(rivermapmud, 2, world.mud3());
    world.setrivermapcolour(rivermapsand, 0, world.sand1());
    world.setrivermapcolour(rivermapsand, 1, world.sand2());
    world.setrivermapcolour(rivermapsand, 2, world.sand3());
    world.setrivermapcolour(rivermapshingle, 0, world.shingle1());
    world.setrivermapcolour(rivermapshingle, 1, world.shingle2());
    world.setrivermapcolour(rivermapshingle, 2, world.shingle3());
    world.setrivermapcolour(rivermapvolcano, 0, 240);
    world.setrivermapcolour(rivermapvolcano, 1, 0);
    world.setrivermapcolour(rivermapvolcano, 2, 0);

    initialisegradientmapappearance(world);

    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
        world.setshowrivermapfeature(i, 1);
}

// This function sets up the default colours (and other settings) for drawing relief maps.

void initialisemapcolours(planet& world)
{
    world.setlandshading(0.7f);
    world.setlakeshading(0.3f);
    world.setseashading(0.5f);

    world.setlandmarbling(0.5f);
    world.setlakemarbling(0.8f);
    world.setseamarbling(0.8f);

    world.setshadingdir(8);
    world.setsnowchange(1);
    world.setseaiceappearance(1);
    world.setminriverflowglobal(1000000);
    world.setminriverflowregional(150);
    world.setshowmangroves(1);
    world.setcolourcliffs(0);

    world.setseaice1(243);
    world.setseaice2(243);
    world.setseaice3(243); // Colours for sea ice.

    world.setocean1(62); // 62
    world.setocean2(93); // 93
    world.setocean3(149); // 149 // Colours for the shallow ocean.

    world.setdeepocean1(3); // 13
    world.setdeepocean2(24); // 34
    world.setdeepocean3(33); // 63 // Colours for the deep ocean.

    world.setbase1(48);
    world.setbase2(69);
    world.setbase3(13); // These are the basic colours for most land.

    world.setbasetemp1(108);
    world.setbasetemp2(111);
    world.setbasetemp3(51); // Basic colours for temperate regions.

    world.sethighbase1(171);
    world.sethighbase2(153);
    world.sethighbase3(135); // High base colours.

    world.setdesert1(225);
    world.setdesert2(175);
    world.setdesert3(118); // Desert colours.

    world.sethighdesert1(132);
    world.sethighdesert2(102);
    world.sethighdesert3(66); // High desert colours.

    world.setcolddesert1(111);
    world.setcolddesert2(93);
    world.setcolddesert3(72); // Cold desert colours.

    world.setgrass1(62);
    world.setgrass2(81);
    world.setgrass3(41); // Grass colours (in between base and desert).

    world.setcold1(207);
    world.setcold2(219);
    world.setcold3(231); // Cold colours.

    world.settundra1(117);
    world.settundra2(120);
    world.settundra3(62); // Tundra colours (in between base and cold).

    world.seteqtundra1(72);
    world.seteqtundra2(63);
    world.seteqtundra3(32); // Tundra at equatorial latitudes.

    world.setsaltpan1(200);
    world.setsaltpan2(200);
    world.setsaltpan3(200); // Salt pan colours.

    world.seterg1(244);
    world.seterg2(231);
    world.seterg3(90); // Erg colours.

    world.setwetlands1(54);
    world.setwetlands2(73);
    world.setwetlands3(75); // Wetlands colours.

    world.setlake1(93); //107;
    world.setlake2(118); //106;
    world.setlake3(178); //146; // Lake colours.

    world.setriver1(75);
    world.setriver2(104);
    world.setriver3(169); // River colours.

    world.setglacier1(209);
    world.setglacier2(222);
    world.setglacier3(251); // Glacier colours.

    world.setsand1(174);
    world.setsand2(168);
    world.setsand3(84); // Sand colours.

    world.setmud1(94);
    world.setmud2(85);
    world.setmud3(55); // Mud colours.

    world.setshingle1(142);
    world.setshingle2(113);
    world.setshingle3(58); // Shingle colours.

    world.setmangrove1(37);
    world.setmangrove2(70);
    world.setmangrove3(17); // Mangrove colours.

    world.sethighlight1(0);
    world.sethighlight2(255);
    world.sethighlight3(255); // Colours of highlights on the map.

    setdefaultnonreliefmapappearance(world);
}

// This randomises some of the world's properties. (Note: the possible changes are fairly modest.)

void changeworldproperties(planet& world)
{
    long seed = world.seed();
    fast_srand(seed);

    // A run may keep the body of the defaults instead (world.defaults.body_from_seed).
    if (!tuning::world::defaults::bodyFromSeed)
        return;

    // Size

    world.setsize(2); // Large world
    world.setgravity(1.0f);

    if (random(1, 16) == 1) // Medium-sized world
    {
        world.setsize(1);
        world.setgravity(0.4f);
    }

    if (random(1, 25) == 1) // Small world
    {
        world.setsize(0);
        world.setgravity(0.15f);
    }

    // Rotation

    if (random(1, 6) == 1)
        world.setrotation(0);
    else
        world.setrotation(1);

    // Gravity

    float oldgravity = world.gravity();

    int gravityvar = random(75, 125);

    world.setgravity(oldgravity * ((float)gravityvar / 100.0f));

    // Perihelion

    if (random(1, 6) == 1)
        world.setperihelion(1);

    // Obliquity

    float newtilt = 25.0f - (float)random(1, 50);
    newtilt = newtilt / 10.0f;

    newtilt = newtilt + 22.5f;

    world.settilt(newtilt);

    // Eccentricity

    float neweccentricity = (float)random(1, 200);
    neweccentricity = neweccentricity / 10000.0f + 0.005f;

    world.seteccentricity(neweccentricity);

    // Lunar pull

    float newlunar = (float)random(1, 40);

    newlunar = newlunar / 100.0f;
    newlunar = newlunar + 0.8f;

    world.setlunar(newlunar);

    // Average temperature

    int averagetemp  = random(1, 30);
    averagetemp = averagetemp - 1;

    world.setaveragetemp(averagetemp);

    // Polar temperature adjustments

    int npoleadjust = randomsign(random(0, 4));
    int spoleadjust = randomsign(random(0, 4));

    world.setnorthpolaradjust(npoleadjust);
    world.setsouthpolaradjust(spoleadjust);

    // Temperature decrease rate

    float newtempdecrease = (float)random(1, 100);
    newtempdecrease = newtempdecrease / 100.0f;

    newtempdecrease = newtempdecrease + 6.0f;

    world.settempdecrease(newtempdecrease);

    // Water pickup rate

    float newwater = 80.0f + (float)random(1, 40);

    newwater = newwater / 100.0f;

    world.setwaterpickup(newwater);

    // Glacial temperature

    int newglacial = -10 + random(1, 20);

    world.setglacialtemp(newglacial);

}

// This works out the total areas of land and sea in the world.

void getlandandseatotals(planet& world)
{
    int width = world.width();
    int height = world.height();
    int sealevel = world.sealevel();
    
    int seatotal = 0;
    int landtotal = 0;

    for (int i = 1; i < width; i++) // Avoid the margins, to avoid false positives.
    {
        for (int j = 1; j < height; j++)
        {
            if (!world.basebelowsea(i, j))
                landtotal++;
            else
                seatotal++;
        }
    }

    if (landtotal > 0)
    {
        if ((float)seatotal / (float)landtotal < 0.05) // If it's this close to 0, assume it is 0.
            seatotal = 0;
    }

    world.setlandtotal(landtotal);
    world.setseatotal(seatotal);
}
