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

const std::array<ColourBinding, 24> colourbindings =
{ {
    { &AppearanceSettings::oceancolour,      { &planet::ocean1,      &planet::ocean2,      &planet::ocean3 },      { &planet::setocean1,      &planet::setocean2,      &planet::setocean3 },      true },
    { &AppearanceSettings::deepoceancolour,  { &planet::deepocean1,  &planet::deepocean2,  &planet::deepocean3 },  { &planet::setdeepocean1,  &planet::setdeepocean2,  &planet::setdeepocean3 },  true },
    { &AppearanceSettings::basecolour,       { &planet::base1,       &planet::base2,       &planet::base3 },       { &planet::setbase1,       &planet::setbase2,       &planet::setbase3 },       true },
    { &AppearanceSettings::grasscolour,      { &planet::grass1,      &planet::grass2,      &planet::grass3 },      { &planet::setgrass1,      &planet::setgrass2,      &planet::setgrass3 },      true },
    { &AppearanceSettings::basetempcolour,   { &planet::basetemp1,   &planet::basetemp2,   &planet::basetemp3 },   { &planet::setbasetemp1,   &planet::setbasetemp2,   &planet::setbasetemp3 },   true },
    { &AppearanceSettings::highbasecolour,   { &planet::highbase1,   &planet::highbase2,   &planet::highbase3 },   { &planet::sethighbase1,   &planet::sethighbase2,   &planet::sethighbase3 },   true },
    { &AppearanceSettings::desertcolour,     { &planet::desert1,     &planet::desert2,     &planet::desert3 },     { &planet::setdesert1,     &planet::setdesert2,     &planet::setdesert3 },     true },
    { &AppearanceSettings::highdesertcolour, { &planet::highdesert1, &planet::highdesert2, &planet::highdesert3 }, { &planet::sethighdesert1, &planet::sethighdesert2, &planet::sethighdesert3 }, true },
    { &AppearanceSettings::colddesertcolour, { &planet::colddesert1, &planet::colddesert2, &planet::colddesert3 }, { &planet::setcolddesert1, &planet::setcolddesert2, &planet::setcolddesert3 }, true },
    { &AppearanceSettings::eqtundracolour,   { &planet::eqtundra1,   &planet::eqtundra2,   &planet::eqtundra3 },   { &planet::seteqtundra1,   &planet::seteqtundra2,   &planet::seteqtundra3 },   true },
    { &AppearanceSettings::tundracolour,     { &planet::tundra1,     &planet::tundra2,     &planet::tundra3 },     { &planet::settundra1,     &planet::settundra2,     &planet::settundra3 },     true },
    { &AppearanceSettings::coldcolour,       { &planet::cold1,       &planet::cold2,       &planet::cold3 },       { &planet::setcold1,       &planet::setcold2,       &planet::setcold3 },       true },
    { &AppearanceSettings::seaicecolour,     { &planet::seaice1,     &planet::seaice2,     &planet::seaice3 },     { &planet::setseaice1,     &planet::setseaice2,     &planet::setseaice3 },     true },
    { &AppearanceSettings::glaciercolour,    { &planet::glacier1,    &planet::glacier2,    &planet::glacier3 },    { &planet::setglacier1,    &planet::setglacier2,    &planet::setglacier3 },    true },
    { &AppearanceSettings::saltpancolour,    { &planet::saltpan1,    &planet::saltpan2,    &planet::saltpan3 },    { &planet::setsaltpan1,    &planet::setsaltpan2,    &planet::setsaltpan3 },    true },
    { &AppearanceSettings::ergcolour,        { &planet::erg1,        &planet::erg2,        &planet::erg3 },        { &planet::seterg1,        &planet::seterg2,        &planet::seterg3 },        true },
    { &AppearanceSettings::wetlandscolour,   { &planet::wetlands1,   &planet::wetlands2,   &planet::wetlands3 },   { &planet::setwetlands1,   &planet::setwetlands2,   &planet::setwetlands3 },   true },
    { &AppearanceSettings::lakecolour,       { &planet::lake1,       &planet::lake2,       &planet::lake3 },       { &planet::setlake1,       &planet::setlake2,       &planet::setlake3 },       true },
    { &AppearanceSettings::rivercolour,      { &planet::river1,      &planet::river2,      &planet::river3 },      { &planet::setriver1,      &planet::setriver2,      &planet::setriver3 },      true },
    { &AppearanceSettings::sandcolour,       { &planet::sand1,       &planet::sand2,       &planet::sand3 },       { &planet::setsand1,       &planet::setsand2,       &planet::setsand3 },       true },
    { &AppearanceSettings::mudcolour,        { &planet::mud1,        &planet::mud2,        &planet::mud3 },        { &planet::setmud1,        &planet::setmud2,        &planet::setmud3 },        true },
    { &AppearanceSettings::shinglecolour,    { &planet::shingle1,    &planet::shingle2,    &planet::shingle3 },    { &planet::setshingle1,    &planet::setshingle2,    &planet::setshingle3 },    true },
    { &AppearanceSettings::mangrovecolour,   { &planet::mangrove1,   &planet::mangrove2,   &planet::mangrove3 },   { &planet::setmangrove1,   &planet::setmangrove2,   &planet::setmangrove3 },   true },
    { &AppearanceSettings::highlightcolour,  { &planet::highlight1,  &planet::highlight2,  &planet::highlight3 },  { &planet::sethighlight1,  &planet::sethighlight2,  &planet::sethighlight3 },  false },
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
    world.setcolourcliffs(settings.colourcliffs);
    world.setshowmangroves(settings.mangroves);

    syncappearancesettings(world, settings);
}
