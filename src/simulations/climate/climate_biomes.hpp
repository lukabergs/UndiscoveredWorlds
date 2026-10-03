#pragma once

// Holdridge biome slots shared by classification (climate) and map colours (presentation).
// Codes are stored in saves and bundles; keep the numbering stable.

constexpr int BIOMEMAPCOLOURCOUNT = 44;

enum BiomeMapColourSlot
{
    biomeice = 0,
    biomepolardesert = 1,
    biomepolardrytundra = 2,
    biomepolarmoisttundra = 3,
    biomepolarwettundra = 4,
    biomepolarraindtundra = 5,
    biomesubpolardesert = 6,
    biomesubpolardrytundra = 7,
    biomesubpolarmoisttundra = 8,
    biomesubpolarwettundra = 9,
    biomesubpolarraindtundra = 10,
    biomeborealdesert = 11,
    biomeborealdrybush = 12,
    biomeborealmoistforest = 13,
    biomeborealwetforest = 14,
    biomeborealrainforest = 15,
    biomecooltemperatedesert = 16,
    biomecooltemperatedesertbush = 17,
    biomecooltemperatesteppe = 18,
    biomecooltemperatemoistforest = 19,
    biomecooltemperatewetforest = 20,
    biomecooltemperaterainforest = 21,
    biomewarmtemperatedesert = 22,
    biomewarmtemperatedesertbush = 23,
    biomewarmtemperatethornsteppe = 24,
    biomewarmtemperatedryforest = 25,
    biomewarmtemperatemoistforest = 26,
    biomewarmtemperatewetforest = 27,
    biomewarmtemperaterainforest = 28,
    biomesubtropicaldesert = 29,
    biomesubtropicaldesertbush = 30,
    biomesubtropicalthornsteppe = 31,
    biomesubtropicaldryforest = 32,
    biomesubtropicalmoistforest = 33,
    biomesubtropicalwetforest = 34,
    biomesubtropicalrainforest = 35,
    biometropicaldesert = 36,
    biometropicaldesertbush = 37,
    biometropicalthornsteppe = 38,
    biometropicalverydryforest = 39,
    biometropicaldryforest = 40,
    biometropicalmoistforest = 41,
    biometropicalwetforest = 42,
    biometropicalrainforest = 43
};
