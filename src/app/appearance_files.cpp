// Saving and loading the map appearance settings (.uws files) of the desktop application.
// Moved unchanged from misc.cpp (Jonathan Hill, 2022).

#include "app/appearance_files.hpp"

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <queue>

#include "core/classes.hpp"
#include "world/generation_tuning.hpp"
#include "world/planet.hpp"
#include "world/region.hpp"
#include "world/functions.hpp"

using namespace std;

// This function saves settings.

void savesettings(planet& world, string filename)
{
    ofstream outfile;
    outfile.open(filename, ios::out);

    outfile << world.settingssaveversion() << '\n';
    outfile << world.landshading() << '\n';
    outfile << world.lakeshading() << '\n';
    outfile << world.seashading() << '\n';
    outfile << world.shadingdir() << '\n';
    outfile << world.snowchange() << '\n';
    outfile << world.seaiceappearance() << '\n';
    outfile << world.landmarbling() << '\n';
    outfile << world.lakemarbling() << '\n';
    outfile << world.seamarbling() << '\n';
    outfile << world.minriverflowglobal() << '\n';
    outfile << world.minriverflowregional() << '\n';
    outfile << world.showmangroves() << '\n';
    outfile << world.colourcliffs() << '\n';
    outfile << world.seaice1() << '\n';
    outfile << world.seaice2() << '\n';
    outfile << world.seaice3() << '\n';
    outfile << world.ocean1() << '\n';
    outfile << world.ocean2() << '\n';
    outfile << world.ocean3() << '\n';
    outfile << world.deepocean1() << '\n';
    outfile << world.deepocean2() << '\n';
    outfile << world.deepocean3() << '\n';
    outfile << world.base1() << '\n';
    outfile << world.base2() << '\n';
    outfile << world.base3() << '\n';
    outfile << world.basetemp1() << '\n';
    outfile << world.basetemp2() << '\n';
    outfile << world.basetemp3() << '\n';
    outfile << world.highbase1() << '\n';
    outfile << world.highbase2() << '\n';
    outfile << world.highbase3() << '\n';
    outfile << world.desert1() << '\n';
    outfile << world.desert2() << '\n';
    outfile << world.desert3() << '\n';
    outfile << world.highdesert1() << '\n';
    outfile << world.highdesert2() << '\n';
    outfile << world.highdesert3() << '\n';
    outfile << world.colddesert1() << '\n';
    outfile << world.colddesert2() << '\n';
    outfile << world.colddesert3() << '\n';
    outfile << world.grass1() << '\n';
    outfile << world.grass2() << '\n';
    outfile << world.grass3() << '\n';
    outfile << world.cold1() << '\n';
    outfile << world.cold2() << '\n';
    outfile << world.cold3() << '\n';
    outfile << world.tundra1() << '\n';
    outfile << world.tundra2() << '\n';
    outfile << world.tundra3() << '\n';
    outfile << world.eqtundra1() << '\n';
    outfile << world.eqtundra2() << '\n';
    outfile << world.eqtundra3() << '\n';
    outfile << world.saltpan1() << '\n';
    outfile << world.saltpan2() << '\n';
    outfile << world.saltpan3() << '\n';
    outfile << world.erg1() << '\n';
    outfile << world.erg2() << '\n';
    outfile << world.erg3() << '\n';
    outfile << world.wetlands1() << '\n';
    outfile << world.wetlands2() << '\n';
    outfile << world.wetlands3() << '\n';
    outfile << world.lake1() << '\n';
    outfile << world.lake2() << '\n';
    outfile << world.lake3() << '\n';
    outfile << world.river1() << '\n';
    outfile << world.river2() << '\n';
    outfile << world.river3() << '\n';
    outfile << world.glacier1() << '\n';
    outfile << world.glacier2() << '\n';
    outfile << world.glacier3() << '\n';
    outfile << world.sand1() << '\n';
    outfile << world.sand2() << '\n';
    outfile << world.sand3() << '\n';
    outfile << world.mud1() << '\n';
    outfile << world.mud2() << '\n';
    outfile << world.mud3() << '\n';
    outfile << world.shingle1() << '\n';
    outfile << world.shingle2() << '\n';
    outfile << world.shingle3() << '\n';
    outfile << world.mangrove1() << '\n';
    outfile << world.mangrove2() << '\n';
    outfile << world.mangrove3() << '\n';
    outfile << world.highlight1() << '\n';
    outfile << world.highlight2() << '\n';
    outfile << world.highlight3() << '\n';
    outfile << world.showmapoutline() << '\n';
    outfile << world.outline1() << '\n';
    outfile << world.outline2() << '\n';
    outfile << world.outline3() << '\n';
    outfile << world.elevationlow1() << '\n';
    outfile << world.elevationlow2() << '\n';
    outfile << world.elevationlow3() << '\n';
    outfile << world.elevationhigh1() << '\n';
    outfile << world.elevationhigh2() << '\n';
    outfile << world.elevationhigh3() << '\n';
    outfile << world.temperaturecold1() << '\n';
    outfile << world.temperaturecold2() << '\n';
    outfile << world.temperaturecold3() << '\n';
    outfile << world.temperaturetemperate1() << '\n';
    outfile << world.temperaturetemperate2() << '\n';
    outfile << world.temperaturetemperate3() << '\n';
    outfile << world.temperaturehot1() << '\n';
    outfile << world.temperaturehot2() << '\n';
    outfile << world.temperaturehot3() << '\n';
    outfile << world.precipitationdry1() << '\n';
    outfile << world.precipitationdry2() << '\n';
    outfile << world.precipitationdry3() << '\n';
    outfile << world.precipitationwet1() << '\n';
    outfile << world.precipitationwet2() << '\n';
    outfile << world.precipitationwet3() << '\n';

    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            outfile << world.climatemapseacolour(i, j) << '\n';
    }

    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            outfile << world.climatemapcolour(i, j) << '\n';
    }

    for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            outfile << world.rivermapcolour(i, j) << '\n';
    }

    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
        outfile << world.showrivermapfeature(i) << '\n';

    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        outfile << world.mapgradientstopcount(gradient) << '\n';
        outfile << world.mapgradientdiscrete(gradient) << '\n';

        for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            outfile << world.mapgradientposition(gradient, stop) << '\n';

            for (int channel = 0; channel < 3; channel++)
                outfile << world.mapgradientcolour(gradient, stop, channel) << '\n';
        }
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        for (int j = 0; j < 3; j++)
            outfile << world.biomemapcolour(i, j) << '\n';
    }

    outfile.close();
}

// This function loads settings.

bool loadsettings(planet& world, string filename)
{
    ifstream infile;
    infile.open(filename, ios::in);

    string line;
    int fileversion;
    int val;
    float fval;
    bool bval;

    getline(infile, line);
    fileversion = stoi(line);

    if (fileversion < 1 || fileversion > world.settingssaveversion()) // Incompatible file format!
        return 0;

    getline(infile, line);
    fval = stof(line);
    world.setlandshading(fval);

    getline(infile, line);
    fval = stof(line);
    world.setlakeshading(fval);

    getline(infile, line);
    fval = stof(line);
    world.setseashading(fval);

    getline(infile, line);
    val = stoi(line);
    world.setshadingdir(val);

    getline(infile, line);
    val = stoi(line);
    world.setsnowchange(val);

    getline(infile, line);
    val = stoi(line);
    world.setseaiceappearance(val);

    getline(infile, line);
    fval = stof(line);
    world.setlandmarbling(fval);

    getline(infile, line);
    fval = stof(line);
    world.setlakemarbling(fval);

    getline(infile, line);
    fval = stof(line);
    world.setseamarbling(fval);

    getline(infile, line);
    val = stoi(line);
    world.setminriverflowglobal(val);

    getline(infile, line);
    val = stoi(line);
    world.setminriverflowregional(val);

    getline(infile, line);
    bval = stob(line);
    world.setshowmangroves(bval);

    getline(infile, line);
    bval = stob(line);
    world.setcolourcliffs(bval);

    getline(infile, line);
    val = stoi(line);
    world.setseaice1(val);

    getline(infile, line);
    val = stoi(line);
    world.setseaice2(val);

    getline(infile, line);
    val = stoi(line);
    world.setseaice3(val);

    getline(infile, line);
    val = stoi(line);
    world.setocean1(val);

    getline(infile, line);
    val = stoi(line);
    world.setocean2(val);

    getline(infile, line);
    val = stoi(line);
    world.setocean3(val);

    getline(infile, line);
    val = stoi(line);
    world.setdeepocean1(val);

    getline(infile, line);
    val = stoi(line);
    world.setdeepocean2(val);

    getline(infile, line);
    val = stoi(line);
    world.setdeepocean3(val);

    getline(infile, line);
    val = stoi(line);
    world.setbase1(val);

    getline(infile, line);
    val = stoi(line);
    world.setbase2(val);

    getline(infile, line);
    val = stoi(line);
    world.setbase3(val);

    getline(infile, line);
    val = stoi(line);
    world.setbasetemp1(val);

    getline(infile, line);
    val = stoi(line);
    world.setbasetemp2(val);

    getline(infile, line);
    val = stoi(line);
    world.setbasetemp3(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighbase1(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighbase2(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighbase3(val);

    getline(infile, line);
    val = stoi(line);
    world.setdesert1(val);

    getline(infile, line);
    val = stoi(line);
    world.setdesert2(val);

    getline(infile, line);
    val = stoi(line);
    world.setdesert3(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighdesert1(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighdesert2(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighdesert3(val);

    getline(infile, line);
    val = stoi(line);
    world.setcolddesert1(val);

    getline(infile, line);
    val = stoi(line);
    world.setcolddesert2(val);

    getline(infile, line);
    val = stoi(line);
    world.setcolddesert3(val);

    getline(infile, line);
    val = stoi(line);
    world.setgrass1(val);

    getline(infile, line);
    val = stoi(line);
    world.setgrass2(val);

    getline(infile, line);
    val = stoi(line);
    world.setgrass3(val);

    getline(infile, line);
    val = stoi(line);
    world.setcold1(val);

    getline(infile, line);
    val = stoi(line);
    world.setcold2(val);

    getline(infile, line);
    val = stoi(line);
    world.setcold3(val);

    getline(infile, line);
    val = stoi(line);
    world.settundra1(val);

    getline(infile, line);
    val = stoi(line);
    world.settundra2(val);

    getline(infile, line);
    val = stoi(line);
    world.settundra3(val);

    getline(infile, line);
    val = stoi(line);
    world.seteqtundra1(val);

    getline(infile, line);
    val = stoi(line);
    world.seteqtundra2(val);

    getline(infile, line);
    val = stoi(line);
    world.seteqtundra3(val);

    getline(infile, line);
    val = stoi(line);
    world.setsaltpan1(val);

    getline(infile, line);
    val = stoi(line);
    world.setsaltpan2(val);

    getline(infile, line);
    val = stoi(line);
    world.setsaltpan3(val);

    getline(infile, line);
    val = stoi(line);
    world.seterg1(val);

    getline(infile, line);
    val = stoi(line);
    world.seterg2(val);

    getline(infile, line);
    val = stoi(line);
    world.seterg3(val);

    getline(infile, line);
    val = stoi(line);
    world.setwetlands1(val);

    getline(infile, line);
    val = stoi(line);
    world.setwetlands2(val);

    getline(infile, line);
    val = stoi(line);
    world.setwetlands3(val);

    getline(infile, line);
    val = stoi(line);
    world.setlake1(val);

    getline(infile, line);
    val = stoi(line);
    world.setlake2(val);

    getline(infile, line);
    val = stoi(line);
    world.setlake3(val);

    getline(infile, line);
    val = stoi(line);
    world.setriver1(val);

    getline(infile, line);
    val = stoi(line);
    world.setriver2(val);

    getline(infile, line);
    val = stoi(line);
    world.setriver3(val);

    getline(infile, line);
    val = stoi(line);
    world.setglacier1(val);

    getline(infile, line);
    val = stoi(line);
    world.setglacier2(val);

    getline(infile, line);
    val = stoi(line);
    world.setglacier3(val);

    getline(infile, line);
    val = stoi(line);
    world.setsand1(val);

    getline(infile, line);
    val = stoi(line);
    world.setsand2(val);

    getline(infile, line);
    val = stoi(line);
    world.setsand3(val);

    getline(infile, line);
    val = stoi(line);
    world.setmud1(val);

    getline(infile, line);
    val = stoi(line);
    world.setmud2(val);

    getline(infile, line);
    val = stoi(line);
    world.setmud3(val);

    getline(infile, line);
    val = stoi(line);
    world.setshingle1(val);

    getline(infile, line);
    val = stoi(line);
    world.setshingle2(val);

    getline(infile, line);
    val = stoi(line);
    world.setshingle3(val);

    getline(infile, line);
    val = stoi(line);
    world.setmangrove1(val);

    getline(infile, line);
    val = stoi(line);
    world.setmangrove2(val);

    getline(infile, line);
    val = stoi(line);
    world.setmangrove3(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighlight1(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighlight2(val);

    getline(infile, line);
    val = stoi(line);
    world.sethighlight3(val);

    if (fileversion >= 2)
    {
        getline(infile, line);
        bval = stob(line);
        world.setshowmapoutline(bval);

        getline(infile, line);
        val = stoi(line);
        world.setoutline1(val);

        getline(infile, line);
        val = stoi(line);
        world.setoutline2(val);

        getline(infile, line);
        val = stoi(line);
        world.setoutline3(val);

        getline(infile, line);
        val = stoi(line);
        world.setelevationlow1(val);

        getline(infile, line);
        val = stoi(line);
        world.setelevationlow2(val);

        getline(infile, line);
        val = stoi(line);
        world.setelevationlow3(val);

        getline(infile, line);
        val = stoi(line);
        world.setelevationhigh1(val);

        getline(infile, line);
        val = stoi(line);
        world.setelevationhigh2(val);

        getline(infile, line);
        val = stoi(line);
        world.setelevationhigh3(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturecold1(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturecold2(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturecold3(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturetemperate1(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturetemperate2(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturetemperate3(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturehot1(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturehot2(val);

        getline(infile, line);
        val = stoi(line);
        world.settemperaturehot3(val);

        getline(infile, line);
        val = stoi(line);
        world.setprecipitationdry1(val);

        getline(infile, line);
        val = stoi(line);
        world.setprecipitationdry2(val);

        getline(infile, line);
        val = stoi(line);
        world.setprecipitationdry3(val);

        getline(infile, line);
        val = stoi(line);
        world.setprecipitationwet1(val);

        getline(infile, line);
        val = stoi(line);
        world.setprecipitationwet2(val);

        getline(infile, line);
        val = stoi(line);
        world.setprecipitationwet3(val);

        if (fileversion >= 3)
        {
            for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    getline(infile, line);
                    world.setclimatemapseacolour(i, j, stoi(line));
                }
            }

            for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    getline(infile, line);
                    world.setclimatemapcolour(i, j, stoi(line));
                }
            }

            for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    getline(infile, line);
                    world.setrivermapcolour(i, j, stoi(line));
                }
            }

            for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
            {
                getline(infile, line);
                world.setshowrivermapfeature(i, stob(line));
            }

            if (fileversion >= 4)
            {
                for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
                {
                    getline(infile, line);
                    world.setmapgradientstopcount(gradient, stoi(line));
                    getline(infile, line);
                    world.setmapgradientdiscrete(gradient, stob(line));

                    for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
                    {
                        getline(infile, line);
                        world.setmapgradientposition(gradient, stop, stoi(line));

                        for (int channel = 0; channel < 3; channel++)
                        {
                            getline(infile, line);
                            world.setmapgradientcolour(gradient, stop, channel, stoi(line));
                        }
                    }
                }
            }
            else
                initialisegradientmapappearance(world);

            if (fileversion >= 5)
            {
                for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        getline(infile, line);
                        world.setbiomemapcolour(i, j, stoi(line));
                    }
                }
            }
            else
            {
                for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
                {
                    for (int j = 0; j < 3; j++)
                        world.setbiomemapcolour(i, j, defaultbiomemapcolours[i][j]);
                }
            }
        }
        else
        {
            static const int defaultclimatecolours[CLIMATEMAPCOLOURCOUNT][3] =
            {
                { 0, 0, 0 }, { 0, 0, 254 }, { 1, 119, 255 }, { 70, 169, 250 }, { 70, 169, 250 }, { 249, 15, 0 }, { 251, 150, 149 }, { 245, 163, 1 },
                { 254, 219, 99 }, { 255, 255, 0 }, { 198, 199, 1 }, { 184, 184, 114 }, { 138, 255, 162 }, { 86, 199, 112 }, { 30, 150, 66 }, { 192, 254, 109 },
                { 76, 255, 93 }, { 19, 203, 74 }, { 255, 8, 245 }, { 204, 3, 192 }, { 154, 51, 144 }, { 153, 100, 146 }, { 172, 178, 249 }, { 91, 121, 213 },
                { 78, 83, 175 }, { 54, 3, 130 }, { 0, 255, 245 }, { 32, 200, 250 }, { 0, 126, 125 }, { 0, 69, 92 }, { 178, 178, 178 }, { 104, 104, 104 }
            };
            static const int defaultclimateseacolours[CLIMATEMAPSEACOLOURCOUNT][3] =
            {
                { 13, 49, 109 }, { 228, 228, 255 }, { 243, 243, 255 }
            };

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

            for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
                world.setshowrivermapfeature(i, 1);

            initialisegradientmapappearance(world);

            for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
            {
                for (int j = 0; j < 3; j++)
                    world.setbiomemapcolour(i, j, defaultbiomemapcolours[i][j]);
            }
        }
    }
    else
        setdefaultnonreliefmapappearance(world);

    infile.close();

    return 1;
}
