#include "appearance_settings.hpp"

#include <array>

namespace
{
using ColourGetter = int (planet::*)() const;
using ColourSetter = void (planet::*)(int);

struct ColourBinding
{
    ImVec4 AppearanceSettings::* colour;
    std::array<ColourGetter, 3> getters;
    std::array<ColourSetter, 3> setters;
    bool redrawMaps = true;
};

const std::array<ColourBinding, 32> colourbindings =
{ {
    { &AppearanceSettings::oceancolour,             { &planet::ocean1,                &planet::ocean2,                &planet::ocean3 },                { &planet::setocean1,                &planet::setocean2,                &planet::setocean3 },                true },
    { &AppearanceSettings::deepoceancolour,         { &planet::deepocean1,            &planet::deepocean2,            &planet::deepocean3 },            { &planet::setdeepocean1,            &planet::setdeepocean2,            &planet::setdeepocean3 },            true },
    { &AppearanceSettings::basecolour,              { &planet::base1,                 &planet::base2,                 &planet::base3 },                 { &planet::setbase1,                 &planet::setbase2,                 &planet::setbase3 },                 true },
    { &AppearanceSettings::grasscolour,             { &planet::grass1,                &planet::grass2,                &planet::grass3 },                { &planet::setgrass1,                &planet::setgrass2,                &planet::setgrass3 },                true },
    { &AppearanceSettings::basetempcolour,          { &planet::basetemp1,             &planet::basetemp2,             &planet::basetemp3 },             { &planet::setbasetemp1,             &planet::setbasetemp2,             &planet::setbasetemp3 },             true },
    { &AppearanceSettings::highbasecolour,          { &planet::highbase1,             &planet::highbase2,             &planet::highbase3 },             { &planet::sethighbase1,             &planet::sethighbase2,             &planet::sethighbase3 },             true },
    { &AppearanceSettings::desertcolour,            { &planet::desert1,               &planet::desert2,               &planet::desert3 },               { &planet::setdesert1,               &planet::setdesert2,               &planet::setdesert3 },               true },
    { &AppearanceSettings::highdesertcolour,        { &planet::highdesert1,           &planet::highdesert2,           &planet::highdesert3 },           { &planet::sethighdesert1,           &planet::sethighdesert2,           &planet::sethighdesert3 },           true },
    { &AppearanceSettings::colddesertcolour,        { &planet::colddesert1,           &planet::colddesert2,           &planet::colddesert3 },           { &planet::setcolddesert1,           &planet::setcolddesert2,           &planet::setcolddesert3 },           true },
    { &AppearanceSettings::eqtundracolour,          { &planet::eqtundra1,             &planet::eqtundra2,             &planet::eqtundra3 },             { &planet::seteqtundra1,             &planet::seteqtundra2,             &planet::seteqtundra3 },             true },
    { &AppearanceSettings::tundracolour,            { &planet::tundra1,               &planet::tundra2,               &planet::tundra3 },               { &planet::settundra1,               &planet::settundra2,               &planet::settundra3 },               true },
    { &AppearanceSettings::coldcolour,              { &planet::cold1,                 &planet::cold2,                 &planet::cold3 },                 { &planet::setcold1,                 &planet::setcold2,                 &planet::setcold3 },                 true },
    { &AppearanceSettings::seaicecolour,            { &planet::seaice1,               &planet::seaice2,               &planet::seaice3 },               { &planet::setseaice1,               &planet::setseaice2,               &planet::setseaice3 },               true },
    { &AppearanceSettings::glaciercolour,           { &planet::glacier1,              &planet::glacier2,              &planet::glacier3 },              { &planet::setglacier1,              &planet::setglacier2,              &planet::setglacier3 },              true },
    { &AppearanceSettings::saltpancolour,           { &planet::saltpan1,              &planet::saltpan2,              &planet::saltpan3 },              { &planet::setsaltpan1,              &planet::setsaltpan2,              &planet::setsaltpan3 },              true },
    { &AppearanceSettings::ergcolour,               { &planet::erg1,                  &planet::erg2,                  &planet::erg3 },                  { &planet::seterg1,                  &planet::seterg2,                  &planet::seterg3 },                  true },
    { &AppearanceSettings::wetlandscolour,          { &planet::wetlands1,             &planet::wetlands2,             &planet::wetlands3 },             { &planet::setwetlands1,             &planet::setwetlands2,             &planet::setwetlands3 },             true },
    { &AppearanceSettings::lakecolour,              { &planet::lake1,                 &planet::lake2,                 &planet::lake3 },                 { &planet::setlake1,                 &planet::setlake2,                 &planet::setlake3 },                 true },
    { &AppearanceSettings::rivercolour,             { &planet::river1,                &planet::river2,                &planet::river3 },                { &planet::setriver1,                &planet::setriver2,                &planet::setriver3 },                true },
    { &AppearanceSettings::sandcolour,              { &planet::sand1,                 &planet::sand2,                 &planet::sand3 },                 { &planet::setsand1,                 &planet::setsand2,                 &planet::setsand3 },                 true },
    { &AppearanceSettings::mudcolour,               { &planet::mud1,                  &planet::mud2,                  &planet::mud3 },                  { &planet::setmud1,                  &planet::setmud2,                  &planet::setmud3 },                  true },
    { &AppearanceSettings::shinglecolour,           { &planet::shingle1,              &planet::shingle2,              &planet::shingle3 },              { &planet::setshingle1,              &planet::setshingle2,              &planet::setshingle3 },              true },
    { &AppearanceSettings::mangrovecolour,          { &planet::mangrove1,             &planet::mangrove2,             &planet::mangrove3 },             { &planet::setmangrove1,             &planet::setmangrove2,             &planet::setmangrove3 },             true },
    { &AppearanceSettings::highlightcolour,         { &planet::highlight1,            &planet::highlight2,            &planet::highlight3 },            { &planet::sethighlight1,            &planet::sethighlight2,            &planet::sethighlight3 },            false },
    { &AppearanceSettings::outlinecolour,           { &planet::outline1,              &planet::outline2,              &planet::outline3 },              { &planet::setoutline1,              &planet::setoutline2,              &planet::setoutline3 },              true },
    { &AppearanceSettings::elevationlowcolour,      { &planet::elevationlow1,         &planet::elevationlow2,         &planet::elevationlow3 },         { &planet::setelevationlow1,         &planet::setelevationlow2,         &planet::setelevationlow3 },         true },
    { &AppearanceSettings::elevationhighcolour,     { &planet::elevationhigh1,        &planet::elevationhigh2,        &planet::elevationhigh3 },        { &planet::setelevationhigh1,        &planet::setelevationhigh2,        &planet::setelevationhigh3 },        true },
    { &AppearanceSettings::temperaturecoldcolour,   { &planet::temperaturecold1,      &planet::temperaturecold2,      &planet::temperaturecold3 },      { &planet::settemperaturecold1,      &planet::settemperaturecold2,      &planet::settemperaturecold3 },      true },
    { &AppearanceSettings::temperaturetemperatecolour, { &planet::temperaturetemperate1, &planet::temperaturetemperate2, &planet::temperaturetemperate3 }, { &planet::settemperaturetemperate1, &planet::settemperaturetemperate2, &planet::settemperaturetemperate3 }, true },
    { &AppearanceSettings::temperaturehotcolour,    { &planet::temperaturehot1,       &planet::temperaturehot2,       &planet::temperaturehot3 },       { &planet::settemperaturehot1,       &planet::settemperaturehot2,       &planet::settemperaturehot3 },       true },
    { &AppearanceSettings::precipitationdrycolour,  { &planet::precipitationdry1,     &planet::precipitationdry2,     &planet::precipitationdry3 },     { &planet::setprecipitationdry1,     &planet::setprecipitationdry2,     &planet::setprecipitationdry3 },     true },
    { &AppearanceSettings::precipitationwetcolour,  { &planet::precipitationwet1,     &planet::precipitationwet2,     &planet::precipitationwet3 },     { &planet::setprecipitationwet1,     &planet::setprecipitationwet2,     &planet::setprecipitationwet3 },     true },
} };

ImVec4 getcolourvalue(const planet& world, const ColourBinding& binding)
{
    return ImVec4(
        (float)(world.*binding.getters[0])() / 255.0f,
        (float)(world.*binding.getters[1])() / 255.0f,
        (float)(world.*binding.getters[2])() / 255.0f,
        1.0f);
}

void synccoloursetting(const planet& world, AppearanceSettings& settings, const ColourBinding& binding)
{
    settings.*(binding.colour) = getcolourvalue(world, binding);
}

bool colourchanged(const planet& world, const AppearanceSettings& settings, const ColourBinding& binding)
{
    const ImVec4& colour = settings.*(binding.colour);

    return (world.*binding.getters[0])() != colour.x * 255.0f
        || (world.*binding.getters[1])() != colour.y * 255.0f
        || (world.*binding.getters[2])() != colour.z * 255.0f;
}

void applycoloursetting(planet& world, const AppearanceSettings& settings, const ColourBinding& binding)
{
    const ImVec4& colour = settings.*(binding.colour);

    (world.*binding.setters[0])((int)(colour.x * 255.0f));
    (world.*binding.setters[1])((int)(colour.y * 255.0f));
    (world.*binding.setters[2])((int)(colour.z * 255.0f));
}

ImVec4 getindexedcolourvalue(int r, int g, int b)
{
    return ImVec4(
        static_cast<float>(r) / 255.0f,
        static_cast<float>(g) / 255.0f,
        static_cast<float>(b) / 255.0f,
        1.0f);
}

bool indexedcolourchanged(int r, int g, int b, const ImVec4& colour)
{
    return r != colour.x * 255.0f
        || g != colour.y * 255.0f
        || b != colour.z * 255.0f;
}

void syncmapgradients(const planet& world, AppearanceSettings& settings)
{
    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        MapGradientSettings& gradientsettings = settings.mapgradients[gradient];
        gradientsettings.stopcount = world.mapgradientstopcount(gradient);
        gradientsettings.discrete = world.mapgradientdiscrete(gradient);

        for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            gradientsettings.stops[stop].position = world.mapgradientposition(gradient, stop);
            gradientsettings.stops[stop].colour = getindexedcolourvalue(
                world.mapgradientcolour(gradient, stop, 0),
                world.mapgradientcolour(gradient, stop, 1),
                world.mapgradientcolour(gradient, stop, 2));
        }
    }
}

bool mapgradientschanged(const planet& world, const AppearanceSettings& settings)
{
    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        const MapGradientSettings& gradientsettings = settings.mapgradients[gradient];

        if (world.mapgradientstopcount(gradient) != gradientsettings.stopcount
            || world.mapgradientdiscrete(gradient) != gradientsettings.discrete)
        {
            return true;
        }

        for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            if (world.mapgradientposition(gradient, stop) != gradientsettings.stops[stop].position
                || indexedcolourchanged(
                    world.mapgradientcolour(gradient, stop, 0),
                    world.mapgradientcolour(gradient, stop, 1),
                    world.mapgradientcolour(gradient, stop, 2),
                    gradientsettings.stops[stop].colour))
            {
                return true;
            }
        }
    }

    return false;
}

void applymapgradients(planet& world, const AppearanceSettings& settings)
{
    for (int gradient = 0; gradient < MAPGRADIENTTYPECOUNT; gradient++)
    {
        const MapGradientSettings& gradientsettings = settings.mapgradients[gradient];
        world.setmapgradientstopcount(gradient, gradientsettings.stopcount);
        world.setmapgradientdiscrete(gradient, gradientsettings.discrete);

        for (int stop = 0; stop < MAPGRADIENTMAXSTOPS; stop++)
        {
            world.setmapgradientposition(gradient, stop, gradientsettings.stops[stop].position);
            world.setmapgradientcolour(gradient, stop, 0, static_cast<int>(gradientsettings.stops[stop].colour.x * 255.0f));
            world.setmapgradientcolour(gradient, stop, 1, static_cast<int>(gradientsettings.stops[stop].colour.y * 255.0f));
            world.setmapgradientcolour(gradient, stop, 2, static_cast<int>(gradientsettings.stops[stop].colour.z * 255.0f));
        }
    }
}

void syncindexedcolours(const planet& world, AppearanceSettings& settings)
{
    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        settings.climateseacolours[i] = getindexedcolourvalue(
            world.climatemapseacolour(i, 0),
            world.climatemapseacolour(i, 1),
            world.climatemapseacolour(i, 2));
    }

    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        settings.climatecolours[i] = getindexedcolourvalue(
            world.climatemapcolour(i, 0),
            world.climatemapcolour(i, 1),
            world.climatemapcolour(i, 2));
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        settings.biomecolours[i] = getindexedcolourvalue(
            world.biomemapcolour(i, 0),
            world.biomemapcolour(i, 1),
            world.biomemapcolour(i, 2));
    }

    for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
    {
        settings.rivermapcolours[i] = getindexedcolourvalue(
            world.rivermapcolour(i, 0),
            world.rivermapcolour(i, 1),
            world.rivermapcolour(i, 2));
    }

    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
        settings.showrivermapfeatures[i] = world.showrivermapfeature(i);
}

bool indexedcolourschanged(const planet& world, const AppearanceSettings& settings)
{
    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        if (indexedcolourchanged(world.climatemapseacolour(i, 0), world.climatemapseacolour(i, 1), world.climatemapseacolour(i, 2), settings.climateseacolours[i]))
            return true;
    }

    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        if (indexedcolourchanged(world.climatemapcolour(i, 0), world.climatemapcolour(i, 1), world.climatemapcolour(i, 2), settings.climatecolours[i]))
            return true;
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        if (indexedcolourchanged(world.biomemapcolour(i, 0), world.biomemapcolour(i, 1), world.biomemapcolour(i, 2), settings.biomecolours[i]))
            return true;
    }

    for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
    {
        if (indexedcolourchanged(world.rivermapcolour(i, 0), world.rivermapcolour(i, 1), world.rivermapcolour(i, 2), settings.rivermapcolours[i]))
            return true;
    }

    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
    {
        if (world.showrivermapfeature(i) != settings.showrivermapfeatures[i])
            return true;
    }

    return false;
}

void applyindexedcolours(planet& world, const AppearanceSettings& settings)
{
    for (int i = 0; i < CLIMATEMAPSEACOLOURCOUNT; i++)
    {
        world.setclimatemapseacolour(i, 0, static_cast<int>(settings.climateseacolours[i].x * 255.0f));
        world.setclimatemapseacolour(i, 1, static_cast<int>(settings.climateseacolours[i].y * 255.0f));
        world.setclimatemapseacolour(i, 2, static_cast<int>(settings.climateseacolours[i].z * 255.0f));
    }

    for (int i = 0; i < CLIMATEMAPCOLOURCOUNT; i++)
    {
        world.setclimatemapcolour(i, 0, static_cast<int>(settings.climatecolours[i].x * 255.0f));
        world.setclimatemapcolour(i, 1, static_cast<int>(settings.climatecolours[i].y * 255.0f));
        world.setclimatemapcolour(i, 2, static_cast<int>(settings.climatecolours[i].z * 255.0f));
    }

    for (int i = 0; i < BIOMEMAPCOLOURCOUNT; i++)
    {
        world.setbiomemapcolour(i, 0, static_cast<int>(settings.biomecolours[i].x * 255.0f));
        world.setbiomemapcolour(i, 1, static_cast<int>(settings.biomecolours[i].y * 255.0f));
        world.setbiomemapcolour(i, 2, static_cast<int>(settings.biomecolours[i].z * 255.0f));
    }

    for (int i = 0; i < RIVERMAPCOLOURCOUNT; i++)
    {
        world.setrivermapcolour(i, 0, static_cast<int>(settings.rivermapcolours[i].x * 255.0f));
        world.setrivermapcolour(i, 1, static_cast<int>(settings.rivermapcolours[i].y * 255.0f));
        world.setrivermapcolour(i, 2, static_cast<int>(settings.rivermapcolours[i].z * 255.0f));
    }

    for (int i = 0; i < RIVERMAPFEATURECOUNT; i++)
        world.setshowrivermapfeature(i, settings.showrivermapfeatures[i]);
}

int getuishadingdir(int shadingdir)
{
    switch (shadingdir)
    {
    case 6:
        return 1;

    case 2:
        return 2;

    case 8:
        return 3;

    case 4:
    default:
        return 0;
    }
}

int getworldshadingdir(int shadingdir)
{
    switch (shadingdir)
    {
    case 1:
        return 6;

    case 2:
        return 2;

    case 3:
        return 8;

    case 0:
    default:
        return 4;
    }
}

void syncscalarsettings(const planet& world, AppearanceSettings& settings)
{
    settings.shadingland = world.landshading();
    settings.shadinglake = world.lakeshading();
    settings.shadingsea = world.seashading();
    settings.marblingland = world.landmarbling();
    settings.marblinglake = world.lakemarbling();
    settings.marblingsea = world.seamarbling();
    settings.globalriversentry = world.minriverflowglobal();
    settings.regionalriversentry = world.minriverflowregional();
    settings.shadingdir = getuishadingdir(world.shadingdir());
    settings.snowchange = world.snowchange() - 1;
    settings.seaiceappearance = world.seaiceappearance() - 1;
    settings.showmapoutline = world.showmapoutline();
    settings.colourcliffs = world.colourcliffs();
    settings.mangroves = world.showmangroves();
}
}

AppearanceSettings makeappearancesettings(const planet& world)
{
    AppearanceSettings settings;
    syncappearancesettings(world, settings);
    return settings;
}

void syncappearancesettings(const planet& world, AppearanceSettings& settings)
{
    for (const ColourBinding& binding : colourbindings)
        synccoloursetting(world, settings, binding);

    syncindexedcolours(world, settings);
    syncmapgradients(world, settings);
    syncscalarsettings(world, settings);
}

AppearanceChangeFlags getappearancechanges(const planet& world, const AppearanceSettings& settings)
{
    AppearanceChangeFlags changes;

    for (const ColourBinding& binding : colourbindings)
    {
        if (!colourchanged(world, settings, binding))
            continue;

        if (binding.redrawMaps)
            changes.mapAppearanceChanged = true;
        else
            changes.highlightChanged = true;
    }

    if (indexedcolourschanged(world, settings))
        changes.mapAppearanceChanged = true;

    if (mapgradientschanged(world, settings))
        changes.mapAppearanceChanged = true;

    if (world.landshading() != settings.shadingland
        || world.lakeshading() != settings.shadinglake
        || world.seashading() != settings.shadingsea
        || world.landmarbling() != settings.marblingland
        || world.lakemarbling() != settings.marblinglake
        || world.seamarbling() != settings.marblingsea
        || world.minriverflowglobal() != settings.globalriversentry
        || world.minriverflowregional() != settings.regionalriversentry
        || world.shadingdir() != getworldshadingdir(settings.shadingdir)
        || world.snowchange() != settings.snowchange + 1
        || world.seaiceappearance() != settings.seaiceappearance + 1
        || world.showmapoutline() != settings.showmapoutline
        || world.colourcliffs() != settings.colourcliffs
        || world.showmangroves() != settings.mangroves)
    {
        changes.mapAppearanceChanged = true;
    }

    return changes;
}

void applyappearancesettings(planet& world, AppearanceSettings& settings)
{
    for (const ColourBinding& binding : colourbindings)
        applycoloursetting(world, settings, binding);

    applyindexedcolours(world, settings);
    applymapgradients(world, settings);
    world.setlandshading(settings.shadingland);
    world.setlakeshading(settings.shadinglake);
    world.setseashading(settings.shadingsea);
    world.setlandmarbling(settings.marblingland);
    world.setlakemarbling(settings.marblinglake);
    world.setseamarbling(settings.marblingsea);
    world.setminriverflowglobal(settings.globalriversentry);
    world.setminriverflowregional(settings.regionalriversentry);
    world.setshadingdir(getworldshadingdir(settings.shadingdir));
    world.setsnowchange(settings.snowchange + 1);
    world.setseaiceappearance(settings.seaiceappearance + 1);
    world.setshowmapoutline(settings.showmapoutline);
    world.setcolourcliffs(settings.colourcliffs);
    world.setshowmangroves(settings.mangroves);

    syncappearancesettings(world, settings);
}
