#pragma once

#include <cstdint>
#include <vector>

struct RgbTriplet
{
    int r = 0;
    int g = 0;
    int b = 0;
};

struct GlobalReliefPalette
{
    std::uint64_t seed = 0;
    int minriverflow = 0;
    int shadingdir = 0;
    int seaiceappearance = 0;
    int snowchange = 0;
    bool colourcliffs = false;
    float landshading = 0.0f;
    float lakeshading = 0.0f;
    float seashading = 0.0f;
    float landmarbling = 0.0f;
    float lakemarbling = 0.0f;
    float seamarbling = 0.0f;
    bool showmangroves = false;
    RgbTriplet seaice;
    RgbTriplet ocean;
    RgbTriplet deepocean;
    RgbTriplet river;
    RgbTriplet saltpan;
    RgbTriplet base;
    RgbTriplet basetemp;
    RgbTriplet highbase;
    RgbTriplet desert;
    RgbTriplet highdesert;
    RgbTriplet colddesert;
    RgbTriplet grass;
    RgbTriplet cold;
    RgbTriplet eqtundra;
    RgbTriplet tundra;
    RgbTriplet erg;
    RgbTriplet wetlands;
    RgbTriplet lake;
    RgbTriplet sand;
    RgbTriplet shingle;
    RgbTriplet mud;
    RgbTriplet mangrove;
    RgbTriplet glacier;
};

struct GlobalCudaRendererInputs
{
    int width = 0;
    int height = 0;
    int sealevel = 0;
    int maxelevation = 0;
    int maxriverflow = 0;
    int worldsize = 0;
    int perihelion = 0;
    float eccentricity = 0.0f;
    float tilt = 0.0f;

    const short* mapnom = nullptr;
    const int* lakemap = nullptr;
    const short* oceanridgeheightmap = nullptr;
    const short* mountainheights = nullptr;
    const short* volcanomap = nullptr;
    const short* extraelevmap = nullptr;
    const short* craterrims = nullptr;
    const short* cratercentres = nullptr;
    const short* oceanridgemap = nullptr;
    const short* oceanriftmap = nullptr;
    const short* jantempmap = nullptr;
    const short* jultempmap = nullptr;
    const short* janrainmap = nullptr;
    const short* julrainmap = nullptr;
    const short* seaicemap = nullptr;
    const short* climatemap = nullptr;
    const int* rivermapjan = nullptr;
    const int* rivermapjul = nullptr;
    const short* specials = nullptr;
    const short* deltamapdir = nullptr;
    const int* deltamapjan = nullptr;
    const int* deltamapjul = nullptr;
    const int* riftlakemapsurface = nullptr;
    const bool* noshademap = nullptr;
    GlobalReliefPalette relief;
};

struct GlobalCudaRendererOutputs
{
    std::vector<std::uint8_t> elevation;
    std::vector<std::uint8_t> elevationDisplay;
    std::vector<std::uint8_t> temperature;
    std::vector<std::uint8_t> temperatureDisplay;
    std::vector<std::uint8_t> precipitation;
    std::vector<std::uint8_t> precipitationDisplay;
    std::vector<std::uint8_t> climate;
    std::vector<std::uint8_t> climateDisplay;
    std::vector<std::uint8_t> rivers;
    std::vector<std::uint8_t> riversDisplay;
};

struct RegionalCudaRendererInputs
{
    int sealevel = 0;
    int maxelevation = 0;
    int maxriverflow = 0;
    int regwidthbegin = 0;
    int regwidthend = 0;
    int regheightbegin = 0;
    int regheightend = 0;
    int regionwidth = 0;
    int regionheight = 0;
    int lefty = 0;
    int globalheight = 0;

    const short* map = nullptr;
    const int* lakemap = nullptr;
    const short* jantempmap = nullptr;
    const short* jultempmap = nullptr;
    const short* extrajantempmap = nullptr;
    const short* extrajultempmap = nullptr;
    const short* janrainmap = nullptr;
    const short* julrainmap = nullptr;
    const short* climatemap = nullptr;
    const short* seaicemap = nullptr;
    const short* rivermapdir = nullptr;
    const int* rivermapjan = nullptr;
    const int* rivermapjul = nullptr;
    const short* fakeriversdir = nullptr;
    const int* fakeriversjan = nullptr;
    const int* fakeriversjul = nullptr;
    const short* specials = nullptr;
    const short* deltamapdir = nullptr;
    const int* deltamapjan = nullptr;
    const int* deltamapjul = nullptr;
    const bool* volcanomap = nullptr;
    const bool* mudmap = nullptr;
    const bool* sandmap = nullptr;
    const bool* shinglemap = nullptr;
    const int* testmap = nullptr;
    GlobalReliefPalette relief;
};

struct RegionalCudaRendererOutputs
{
    std::vector<std::uint8_t> elevation;
    std::vector<std::uint8_t> temperature;
    std::vector<std::uint8_t> precipitation;
    std::vector<std::uint8_t> climate;
    std::vector<std::uint8_t> rivers;
};

bool cudaglobalrenderersavailable();
bool renderglobalnonreliefmapscuda(const GlobalCudaRendererInputs& inputs, GlobalCudaRendererOutputs& outputs);
bool renderglobalreliefmapcuda(const GlobalCudaRendererInputs& inputs, std::vector<std::uint8_t>& relief, std::vector<std::uint8_t>& reliefDisplay);
bool renderregionalnonreliefmapscuda(const RegionalCudaRendererInputs& inputs, RegionalCudaRendererOutputs& outputs);
bool renderregionalreliefmapcuda(const RegionalCudaRendererInputs& inputs, std::vector<std::uint8_t>& relief);
