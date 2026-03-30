#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include <cuda_runtime.h>

#include "cuda_renderers.hpp"

namespace
{
    constexpr int cudaplanetstride = 1025;
    constexpr int cudaregionalstride = 1000;
    constexpr int displaywidth = 1024;
    constexpr int displayheight = 512;
    constexpr int rgbachannels = 4;

    struct GlobalCudaDeviceInputs
    {
        int width;
        int height;
        int sealevel;
        int maxelevation;
        int maxriverflow;
        int worldsize;
        int perihelion;
        float eccentricity;
        float tilt;

        const short* mapnom;
        const int* lakemap;
        const short* oceanridgeheightmap;
        const short* mountainheights;
        const short* volcanomap;
        const short* extraelevmap;
        const short* craterrims;
        const short* cratercentres;
        const short* oceanridgemap;
        const short* oceanriftmap;
        const short* jantempmap;
        const short* jultempmap;
        const short* janrainmap;
        const short* julrainmap;
        const short* seaicemap;
        const short* climatemap;
        const int* rivermapjan;
        const int* rivermapjul;
        const short* specials;
        const short* deltamapdir;
        const int* deltamapjan;
        const int* deltamapjul;
        const int* riftlakemapsurface;
        const bool* noshademap;
        GlobalReliefPalette relief;
    };

    struct RegionalCudaDeviceInputs
    {
        int sealevel;
        int maxelevation;
        int maxriverflow;
        int regwidthbegin;
        int regwidthend;
        int regheightbegin;
        int regheightend;
        int regionwidth;
        int regionheight;
        int lefty;
        int globalheight;

        const short* map;
        const int* lakemap;
        const short* jantempmap;
        const short* jultempmap;
        const short* extrajantempmap;
        const short* extrajultempmap;
        const short* janrainmap;
        const short* julrainmap;
        const short* climatemap;
        const short* seaicemap;
        const short* rivermapdir;
        const int* rivermapjan;
        const int* rivermapjul;
        const short* fakeriversdir;
        const int* fakeriversjan;
        const int* fakeriversjul;
        const short* specials;
        const short* deltamapdir;
        const int* deltamapjan;
        const int* deltamapjul;
        const bool* volcanomap;
        const bool* mudmap;
        const bool* sandmap;
        const bool* shinglemap;
        const int* testmap;
        GlobalReliefPalette relief;
    };

    template <typename T>
    class DeviceBuffer
    {
    public:
        DeviceBuffer() = default;

        ~DeviceBuffer()
        {
            if (data_ != nullptr)
                cudaFree(data_);
        }

        bool ensure(size_t count)
        {
            if (count <= count_)
                return true;

            if (data_ != nullptr)
                cudaFree(data_);

            if (cudaMalloc(&data_, count * sizeof(T)) != cudaSuccess)
            {
                data_ = nullptr;
                count_ = 0;
                return false;
            }

            count_ = count;
            return true;
        }

        T* data() { return data_; }
        const T* data() const { return data_; }

    private:
        T* data_ = nullptr;
        size_t count_ = 0;
    };

    bool cudacheck(cudaError_t error)
    {
        return error == cudaSuccess;
    }

    template <typename T>
    bool uploadbuffer(DeviceBuffer<T>& devicebuffer, const T* hostdata, size_t count)
    {
        if (!devicebuffer.ensure(count))
            return false;

        return cudacheck(cudaMemcpy(devicebuffer.data(), hostdata, count * sizeof(T), cudaMemcpyHostToDevice));
    }

    bool copyoutput(std::vector<std::uint8_t>& output, DeviceBuffer<std::uint8_t>& devicebuffer, size_t count)
    {
        output.resize(count);
        return cudacheck(cudaMemcpy(output.data(), devicebuffer.data(), count, cudaMemcpyDeviceToHost));
    }

    __device__ int inputindex(int x, int y)
    {
        return x * cudaplanetstride + y;
    }

    __device__ int regionalinputindex(int x, int y)
    {
        return x * cudaregionalstride + y;
    }

    __device__ size_t outputindex(int width, int x, int y)
    {
        return (static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)) * rgbachannels;
    }

    __device__ int wrapx(int x, int width)
    {
        const int maxvalue = width + 1;
        x = x % maxvalue;

        if (x < 0)
            x += maxvalue;

        return x;
    }

    __device__ int clamp255(int value)
    {
        if (value < 0)
            return 0;

        if (value > 255)
            return 255;

        return value;
    }

    __device__ void writepixel(std::uint8_t* output, int outputwidth, int x, int y, int r, int g, int b, int a = 255)
    {
        const size_t offset = outputindex(outputwidth, x, y);
        output[offset] = static_cast<std::uint8_t>(r);
        output[offset + 1] = static_cast<std::uint8_t>(g);
        output[offset + 2] = static_cast<std::uint8_t>(b);
        output[offset + 3] = static_cast<std::uint8_t>(a);
    }

    __device__ const std::uint8_t* readpixel(const std::uint8_t* pixels, int width, int x, int y)
    {
        return pixels + outputindex(width, x, y);
    }

    __device__ void assignrgb(const RgbTriplet& colour, int& r, int& g, int& b)
    {
        r = colour.r;
        g = colour.g;
        b = colour.b;
    }

    __device__ unsigned int hashcoords(std::uint64_t seed, int x, int y, int salt)
    {
        std::uint64_t value = seed;
        value ^= static_cast<std::uint64_t>(x + 0x9e3779b9u) + (value << 6) + (value >> 2);
        value ^= static_cast<std::uint64_t>(y + 0x85ebca6bu) + (value << 6) + (value >> 2);
        value ^= static_cast<std::uint64_t>(salt + 0xc2b2ae35u) + (value << 6) + (value >> 2);
        value ^= value >> 33;
        value *= 0xff51afd7ed558ccdULL;
        value ^= value >> 33;
        value *= 0xc4ceb9fe1a85ec53ULL;
        value ^= value >> 33;
        return static_cast<unsigned int>(value & 0xffffffffu);
    }

    __device__ int deterministicrandom(std::uint64_t seed, int x, int y, int salt, int minimum, int maximum)
    {
        if (maximum <= minimum)
            return minimum;

        const unsigned int range = static_cast<unsigned int>((maximum - minimum) + 1);
        return minimum + static_cast<int>(hashcoords(seed, x, y, salt) % range);
    }

    __device__ int deterministicrandomsigned(std::uint64_t seed, int x, int y, int salt, int magnitude)
    {
        if (magnitude <= 0)
            return 0;

        int value = deterministicrandom(seed, x, y, salt, 0, magnitude);

        if ((hashcoords(seed, x, y, salt + 97) & 1u) != 0u)
            value = -value;

        return value;
    }

    __device__ int computemap(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);
        const int thisvolcano = abs(static_cast<int>(inputs.volcanomap[index]));

        if (inputs.mapnom[index] <= inputs.sealevel)
            return static_cast<int>(inputs.mapnom[index]) + static_cast<int>(inputs.oceanridgeheightmap[index]) + thisvolcano;

        if (thisvolcano > inputs.mountainheights[index])
            return static_cast<int>(inputs.mapnom[index]) + thisvolcano + static_cast<int>(inputs.extraelevmap[index]) + static_cast<int>(inputs.craterrims[index]) + static_cast<int>(inputs.cratercentres[index]);

        return static_cast<int>(inputs.mapnom[index]) + static_cast<int>(inputs.mountainheights[index]) + static_cast<int>(inputs.extraelevmap[index]) + static_cast<int>(inputs.craterrims[index]) + static_cast<int>(inputs.cratercentres[index]);
    }

    __device__ bool computesea(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);

        if (inputs.mapnom[index] <= inputs.sealevel && inputs.lakemap[index] == 0)
            return true;

        if (inputs.volcanomap[index] == 0)
            return false;

        const int thisvolcano = abs(static_cast<int>(inputs.volcanomap[index]));
        return static_cast<int>(inputs.mapnom[index]) + thisvolcano <= inputs.sealevel;
    }

    __device__ bool computeoutline(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        if (y < 1 || y > inputs.height - 1)
            return false;

        if (computesea(inputs, x, y))
            return false;

        return computesea(inputs, wrapx(x - 1, inputs.width), y)
            || computesea(inputs, x, y - 1)
            || computesea(inputs, wrapx(x + 1, inputs.width), y)
            || computesea(inputs, x, y + 1);
    }

    __device__ int computeaprtemp(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);
        float summertemp = static_cast<float>(inputs.jultempmap[index]);
        float wintertemp = static_cast<float>(inputs.jantempmap[index]);

        if (inputs.perihelion == 1)
        {
            summertemp = static_cast<float>(inputs.jantempmap[index]);
            wintertemp = static_cast<float>(inputs.jultempmap[index]);
        }

        const float winterstrength = 0.5f + inputs.eccentricity * 0.5f;
        const float summerstrength = 1.0f - winterstrength;

        float thistemp = summertemp * summerstrength + wintertemp * winterstrength;
        const float fourseason = inputs.tilt * 0.294592f - 2.45428f;

        float lat = static_cast<float>(y);
        if (y > inputs.height / 2.0f)
            lat = static_cast<float>(inputs.height - y);

        const float fourseasonstrength = lat / (static_cast<float>(inputs.height) / 2.0f);
        thistemp += (fourseason * fourseasonstrength) / 2.0f;

        return static_cast<int>(thistemp);
    }

    __device__ int computewinterrain(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);

        if (inputs.jantempmap[index] < inputs.jultempmap[index])
            return static_cast<int>(inputs.janrainmap[index]);

        return static_cast<int>(inputs.julrainmap[index]);
    }

    __device__ int computesummerrain(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);

        if (inputs.jantempmap[index] >= inputs.jultempmap[index])
            return static_cast<int>(inputs.janrainmap[index]);

        return static_cast<int>(inputs.julrainmap[index]);
    }

    __device__ bool computetruelake(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);
        return inputs.lakemap[index] != 0 && inputs.specials[index] < 110;
    }

    __device__ int computeavetemp(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);
        return (static_cast<int>(inputs.jantempmap[index]) + static_cast<int>(inputs.jultempmap[index])) / 2;
    }

    __device__ int computemintemp(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = inputindex(x, y);
        const int jantemp = static_cast<int>(inputs.jantempmap[index]);
        const int jultemp = static_cast<int>(inputs.jultempmap[index]);
        return jantemp < jultemp ? jantemp : jultemp;
    }

    __device__ bool isseaicevisible(const GlobalCudaDeviceInputs& inputs, int x, int y)
    {
        const int seaice = static_cast<int>(inputs.seaicemap[inputindex(x, y)]);
        return (seaice == 2 && (inputs.relief.seaiceappearance == 1 || inputs.relief.seaiceappearance == 3))
            || (seaice == 1 && inputs.relief.seaiceappearance == 3);
    }

    __device__ int computetundradegrees(const GlobalCudaDeviceInputs& inputs, int y)
    {
        const float fy = static_cast<float>(y);
        const float worldheight = static_cast<float>(inputs.height);
        const float worldhalfheight = worldheight / 2.0f;
        const float pixelsperlat = worldhalfheight / 90.0f;

        float latitude = 0.0f;

        if (fy <= worldhalfheight)
            latitude = (worldhalfheight - fy) / pixelsperlat;
        else
            latitude = (fy - worldhalfheight) / pixelsperlat;

        return static_cast<int>(latitude);
    }

    __device__ RgbTriplet computetundracolour(const GlobalCudaDeviceInputs& inputs, int y)
    {
        const int lat = computetundradegrees(inputs, y);
        const int lat2 = 90 - lat;

        RgbTriplet colour;
        colour.r = (lat * inputs.relief.tundra.r + lat2 * inputs.relief.eqtundra.r) / 90;
        colour.g = (lat * inputs.relief.tundra.g + lat2 * inputs.relief.eqtundra.g) / 90;
        colour.b = (lat * inputs.relief.tundra.b + lat2 * inputs.relief.eqtundra.b) / 90;
        return colour;
    }

    __device__ int getslopecached(const GlobalCudaDeviceInputs& inputs, int x, int y, int currentelevation, bool& valid)
    {
        if (y < 0 || y > inputs.height)
        {
            valid = false;
            return 0;
        }

        valid = true;
        return computemap(inputs, wrapx(x, inputs.width), y) - currentelevation;
    }

    __device__ bool regionalsea(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return inputs.map[index] <= inputs.sealevel && inputs.lakemap[index] == 0;
    }

    __device__ bool regionaloutline(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        if (y < 1 || y > inputs.regionheight - 1)
            return false;

        if (regionalsea(inputs, x, y))
            return false;

        if (x > 0 && regionalsea(inputs, x - 1, y))
            return true;

        if (regionalsea(inputs, x, y - 1))
            return true;

        if (x < inputs.regionwidth && regionalsea(inputs, x + 1, y))
            return true;

        if (regionalsea(inputs, x, y + 1))
            return true;

        return false;
    }

    __device__ int regionalmintemp(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        const int jantemp = static_cast<int>(inputs.jantempmap[index]);
        const int jultemp = static_cast<int>(inputs.jultempmap[index]);
        return jantemp >= jultemp ? jultemp : jantemp;
    }

    __device__ int regionalmaxtemp(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        const int jantemp = static_cast<int>(inputs.jantempmap[index]);
        const int jultemp = static_cast<int>(inputs.jultempmap[index]);
        return jantemp < jultemp ? jultemp : jantemp;
    }

    __device__ int regionalextramintemp(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return inputs.jantempmap[index] >= inputs.jultempmap[index] ? static_cast<int>(inputs.extrajultempmap[index]) : static_cast<int>(inputs.extrajantempmap[index]);
    }

    __device__ int regionalextramaxtemp(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return inputs.jantempmap[index] < inputs.jultempmap[index] ? static_cast<int>(inputs.extrajultempmap[index]) : static_cast<int>(inputs.extrajantempmap[index]);
    }

    __device__ int regionalwinterrain(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return inputs.jantempmap[index] >= inputs.jultempmap[index] ? static_cast<int>(inputs.julrainmap[index]) : static_cast<int>(inputs.janrainmap[index]);
    }

    __device__ int regionalsummerrain(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return inputs.jantempmap[index] < inputs.jultempmap[index] ? static_cast<int>(inputs.julrainmap[index]) : static_cast<int>(inputs.janrainmap[index]);
    }

    __device__ int regionalriveraveflow(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return (inputs.rivermapjan[index] + inputs.rivermapjul[index]) / 2;
    }

    __device__ int regionalfakeaveflow(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return (inputs.fakeriversjan[index] + inputs.fakeriversjul[index]) / 2;
    }

    __device__ bool regionaltruelake(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        const int special = static_cast<int>(inputs.specials[index]);
        return inputs.lakemap[index] != 0 && (special < 110 || special == 140);
    }

    __device__ bool regionalrivervalley(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);
        return inputs.fakeriversdir[index] == -1 && inputs.map[index] > inputs.sealevel && inputs.specials[index] != 140;
    }

    __device__ bool regionalmangrove(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int index = regionalinputindex(x, y);

        if (!inputs.relief.showmangroves)
            return false;

        if (!inputs.mudmap[index] && inputs.specials[index] != 131 && inputs.specials[index] != 132)
            return false;

        if (inputs.sandmap[index] || inputs.shinglemap[index])
            return false;

        if (inputs.jantempmap[index] < 18 || inputs.jultempmap[index] < 18)
            return false;

        return true;
    }

    __device__ bool regionalseaicevisible(const RegionalCudaDeviceInputs& inputs, int x, int y)
    {
        const int seaice = static_cast<int>(inputs.seaicemap[regionalinputindex(x, y)]);
        return (seaice == 2 && (inputs.relief.seaiceappearance == 1 || inputs.relief.seaiceappearance == 3))
            || (seaice == 1 && inputs.relief.seaiceappearance == 3);
    }

    __device__ RgbTriplet computetundracolourforrow(int globalheight, const GlobalReliefPalette& relief, int globalrow)
    {
        const float fy = static_cast<float>(globalrow);
        const float worldheight = static_cast<float>(globalheight);
        const float worldhalfheight = worldheight / 2.0f;
        const float pixelsperlat = worldhalfheight / 90.0f;

        float latitude = 0.0f;

        if (fy <= worldhalfheight)
            latitude = (worldhalfheight - fy) / pixelsperlat;
        else
            latitude = (fy - worldhalfheight) / pixelsperlat;

        const int lat = static_cast<int>(latitude);
        const int lat2 = 90 - lat;

        RgbTriplet colour;
        colour.r = (lat * relief.tundra.r + lat2 * relief.eqtundra.r) / 90;
        colour.g = (lat * relief.tundra.g + lat2 * relief.eqtundra.g) / 90;
        colour.b = (lat * relief.tundra.b + lat2 * relief.eqtundra.b) / 90;
        return colour;
    }

    __device__ void writeshortpixel(short* pixels, int x, int y, short value)
    {
        pixels[regionalinputindex(x, y)] = value;
    }

    __device__ short readshortpixel(const short* pixels, int x, int y)
    {
        return pixels[regionalinputindex(x, y)];
    }

    __global__ void renderelevationkernel(GlobalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x > inputs.width || y > inputs.height)
            return;

        const int div = max(inputs.maxelevation / 255, 1);
        int colour = computemap(inputs, x, y) / div;
        colour = clamp255(colour);

        writepixel(output, outputwidth, x, y, colour, colour, colour);
    }

    __global__ void rendertemperaturekernel(GlobalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x > inputs.width || y > inputs.height)
            return;

        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;

        if (!computeoutline(inputs, x, y))
        {
            const int index = inputindex(x, y);
            int temperature = (static_cast<int>(inputs.jantempmap[index]) + computeaprtemp(inputs, x, y) + static_cast<int>(inputs.jultempmap[index]) + computeaprtemp(inputs, x, y)) / 4;
            temperature += 10;

            if (temperature > 0)
            {
                colour1 = 250;
                colour2 = 250 - (temperature * 3);
                colour3 = 250 - (temperature * 7);
            }
            else
            {
                temperature = abs(temperature);
                colour1 = 250 - (temperature * 7);
                colour2 = 250 - (temperature * 7);
                colour3 = 250;
            }
        }

        writepixel(output, outputwidth, x, y, clamp255(colour1), clamp255(colour2), clamp255(colour3));
    }

    __global__ void renderprecipitationkernel(GlobalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x > inputs.width || y > inputs.height)
            return;

        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;

        if (!computeoutline(inputs, x, y))
        {
            const int rainfall = (computesummerrain(inputs, x, y) + computewinterrain(inputs, x, y)) / 8;
            colour1 = 255 - rainfall;
            colour2 = 255 - rainfall;
            colour3 = 255;
        }

        writepixel(output, outputwidth, x, y, clamp255(colour1), clamp255(colour2), clamp255(colour3));
    }

    __device__ void climatecolour(short climate, int& r, int& g, int& b)
    {
        switch (climate)
        {
        case 1: r = 0; g = 0; b = 254; break;
        case 2: r = 1; g = 119; b = 255; break;
        case 3:
        case 4: r = 70; g = 169; b = 250; break;
        case 5: r = 249; g = 15; b = 0; break;
        case 6: r = 251; g = 150; b = 149; break;
        case 7: r = 245; g = 163; b = 1; break;
        case 8: r = 254; g = 219; b = 99; break;
        case 9: r = 255; g = 255; b = 0; break;
        case 10: r = 198; g = 199; b = 1; break;
        case 11: r = 184; g = 184; b = 114; break;
        case 12: r = 138; g = 255; b = 162; break;
        case 13: r = 86; g = 199; b = 112; break;
        case 14: r = 30; g = 150; b = 66; break;
        case 15: r = 192; g = 254; b = 109; break;
        case 16: r = 76; g = 255; b = 93; break;
        case 17: r = 19; g = 203; b = 74; break;
        case 18: r = 255; g = 8; b = 245; break;
        case 19: r = 204; g = 3; b = 192; break;
        case 20: r = 154; g = 51; b = 144; break;
        case 21: r = 153; g = 100; b = 146; break;
        case 22: r = 172; g = 178; b = 249; break;
        case 23: r = 91; g = 121; b = 213; break;
        case 24: r = 78; g = 83; b = 175; break;
        case 25: r = 54; g = 3; b = 130; break;
        case 26: r = 0; g = 255; b = 245; break;
        case 27: r = 32; g = 200; b = 250; break;
        case 28: r = 0; g = 126; b = 125; break;
        case 29: r = 0; g = 69; b = 92; break;
        case 30: r = 178; g = 178; b = 178; break;
        case 31: r = 104; g = 104; b = 104; break;
        default: r = 0; g = 0; b = 0; break;
        }
    }

    __global__ void renderclimatekernel(GlobalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x > inputs.width || y > inputs.height)
            return;

        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;

        if (computesea(inputs, x, y))
        {
            const int seaice = static_cast<int>(inputs.seaicemap[inputindex(x, y)]);

            if (seaice == 2)
            {
                colour1 = 243;
                colour2 = 243;
                colour3 = 255;
            }
            else if (seaice == 1)
            {
                colour1 = 228;
                colour2 = 228;
                colour3 = 255;
            }
            else
            {
                colour1 = 13;
                colour2 = 49;
                colour3 = 109;
            }
        }
        else
        {
            climatecolour(inputs.climatemap[inputindex(x, y)], colour1, colour2, colour3);
        }

        writepixel(output, outputwidth, x, y, colour1, colour2, colour3);
    }

    __global__ void renderriverskernel(GlobalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x > inputs.width || y > inputs.height)
            return;

        const int index = inputindex(x, y);
        const bool sea = computesea(inputs, x, y);
        int colour1 = 255;
        int colour2 = 255;
        int colour3 = 255;

        if (computeoutline(inputs, x, y))
        {
            colour1 = 0;
            colour2 = 0;
            colour3 = 0;
        }
        else
        {
            int flow = (inputs.rivermapjan[index] + inputs.rivermapjul[index]) / 2;

            if (flow > 0 && !sea)
            {
                const int mult = max(inputs.maxriverflow / 255, 1);
                flow *= 10;
                colour1 = clamp255(255 - (flow / mult));
                colour2 = colour1;
            }
            else if (inputs.deltamapdir[index] != 0 && !sea)
            {
                const int mult = max(inputs.maxriverflow / 255, 1);
                flow = (inputs.deltamapjan[index] + inputs.deltamapjul[index]) / 2;
                flow *= 10;
                colour1 = clamp255(255 - (flow / mult));
                colour2 = colour1;
            }

            if (computetruelake(inputs, x, y))
            {
                colour1 = 150;
                colour2 = 150;
                colour3 = 250;
            }

            if (inputs.specials[index] > 100 && !sea && inputs.rivermapjan[index] + inputs.rivermapjul[index] < 600)
            {
                if (inputs.specials[index] == 110)
                {
                    colour1 = 150;
                    colour2 = 150;
                    colour3 = 150;
                }
                else if (inputs.specials[index] == 120)
                {
                    colour1 = 250;
                    colour2 = 250;
                    colour3 = 50;
                }
                else if (inputs.specials[index] >= 130)
                {
                    colour1 = 50;
                    colour2 = 250;
                    colour3 = 100;
                }
            }
        }

        if (computeoutline(inputs, x, y))
        {
            colour1 = 0;
            colour2 = 0;
            colour3 = 0;
        }

        if (inputs.volcanomap[index] > 0)
        {
            colour1 = 240;
            colour2 = 0;
            colour3 = 0;
        }

        writepixel(output, outputwidth, x, y, colour1, colour2, colour3);
    }

    __global__ void renderglobalreliefkernel(GlobalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x > inputs.width || y > inputs.height)
            return;

        const int index = inputindex(x, y);
        const int mapvalue = computemap(inputs, x, y);
        const bool sea = computesea(inputs, x, y);
        const int special = static_cast<int>(inputs.specials[index]);
        const int riverflow = (inputs.rivermapjan[index] + inputs.rivermapjul[index]) / 2;
        const bool truelake = computetruelake(inputs, x, y);
        const bool riftlake = inputs.riftlakemapsurface[index] != 0;

        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;
        int var = 10;

        if (sea)
        {
            if (isseaicevisible(inputs, x, y))
            {
                assignrgb(inputs.relief.seaice, colour1, colour2, colour3);
                var = 0;
            }
            else
            {
                colour1 = (inputs.relief.ocean.r * mapvalue + inputs.relief.deepocean.r * (inputs.sealevel - mapvalue)) / inputs.sealevel;
                colour2 = (inputs.relief.ocean.g * mapvalue + inputs.relief.deepocean.g * (inputs.sealevel - mapvalue)) / inputs.sealevel;
                colour3 = (inputs.relief.ocean.b * mapvalue + inputs.relief.deepocean.b * (inputs.sealevel - mapvalue)) / inputs.sealevel;
                var = 5;
            }
        }
        else
        {
            if (riverflow >= inputs.relief.minriverflow)
            {
                assignrgb(inputs.relief.river, colour1, colour2, colour3);
            }
            else if (special == 110)
            {
                assignrgb(inputs.relief.saltpan, colour1, colour2, colour3);
                var = 20;
            }
            else
            {
                const int avetemp = computeavetemp(inputs, x, y) + 10;

                int thisbase1 = 0;
                int thisbase2 = 0;
                int thisbase3 = 0;
                int newdesert1 = 0;
                int newdesert2 = 0;
                int newdesert3 = 0;

                if (avetemp > 30)
                {
                    assignrgb(inputs.relief.base, thisbase1, thisbase2, thisbase3);
                    assignrgb(inputs.relief.desert, newdesert1, newdesert2, newdesert3);
                }
                else
                {
                    const int hotno = avetemp / 3;
                    const int coldno = 10 - hotno;

                    thisbase1 = (hotno * inputs.relief.base.r + coldno * inputs.relief.basetemp.r) / 10;
                    thisbase2 = (hotno * inputs.relief.base.g + coldno * inputs.relief.basetemp.g) / 10;
                    thisbase3 = (hotno * inputs.relief.base.b + coldno * inputs.relief.basetemp.b) / 10;
                }

                if (avetemp > 30)
                {
                    assignrgb(inputs.relief.desert, newdesert1, newdesert2, newdesert3);
                }
                else if (avetemp <= 10)
                {
                    assignrgb(inputs.relief.colddesert, newdesert1, newdesert2, newdesert3);
                }
                else
                {
                    const int hotno = avetemp - 10;
                    const int coldno = 20 - hotno;

                    newdesert1 = (hotno * inputs.relief.desert.r + coldno * inputs.relief.colddesert.r) / 20;
                    newdesert2 = (hotno * inputs.relief.desert.g + coldno * inputs.relief.colddesert.g) / 20;
                    newdesert3 = (hotno * inputs.relief.desert.b + coldno * inputs.relief.colddesert.b) / 20;
                }

                float winterrain = static_cast<float>(computewinterrain(inputs, x, y));
                const float summerrain = static_cast<float>(computesummerrain(inputs, x, y));
                float totalrain = winterrain + summerrain;

                if (winterrain < 1.0f)
                    winterrain = 1.0f;

                float rainforestmult = static_cast<float>(computemintemp(inputs, x, y)) / 18.0f;
                rainforestmult = rainforestmult * winterrain / 80.0f;

                if (rainforestmult < 1.0f)
                    rainforestmult = 1.0f;

                totalrain = totalrain * rainforestmult;

                int mapelev = mapvalue - inputs.sealevel;
                const int desertmapelev = mapelev;

                if (inputs.relief.colourcliffs)
                {
                    int biggestslope = 0;

                    for (int k = x - 1; k <= x + 1; k++)
                    {
                        const int kk = wrapx(k, inputs.width);

                        for (int l = y - 1; l <= y + 1; l++)
                        {
                            if (l >= 0 && l <= inputs.height)
                            {
                                const int thisslope = mapvalue - computemap(inputs, kk, l);

                                if (thisslope > biggestslope)
                                    biggestslope = thisslope;
                            }
                        }
                    }

                    biggestslope -= 240;

                    if (biggestslope < 0)
                        biggestslope = 0;

                    const float adjustedelev = static_cast<float>(mapelev) * (static_cast<float>(biggestslope) / 240.0f);
                    mapelev = static_cast<int>(adjustedelev);
                }

                int newbase1 = 0;
                int newbase2 = 0;
                int newbase3 = 0;
                int newgrass1 = 0;
                int newgrass2 = 0;
                int newgrass3 = 0;

                if (desertmapelev > 2000)
                {
                    assignrgb(inputs.relief.highdesert, newdesert1, newdesert2, newdesert3);
                }
                else
                {
                    const int highno = desertmapelev / 50;
                    const int lowno = 40 - highno;

                    newdesert1 = (highno * inputs.relief.highdesert.r + lowno * newdesert1) / 40;
                    newdesert2 = (highno * inputs.relief.highdesert.g + lowno * newdesert2) / 40;
                    newdesert3 = (highno * inputs.relief.highdesert.b + lowno * newdesert3) / 40;
                }

                if (mapelev > 3000)
                {
                    assignrgb(inputs.relief.highbase, newbase1, newbase2, newbase3);
                    assignrgb(inputs.relief.highbase, newgrass1, newgrass2, newgrass3);
                }
                else
                {
                    const int highno = mapelev / 75;
                    const int lowno = 40 - highno;

                    newbase1 = (highno * inputs.relief.highbase.r + lowno * thisbase1) / 40;
                    newbase2 = (highno * inputs.relief.highbase.g + lowno * thisbase2) / 40;
                    newbase3 = (highno * inputs.relief.highbase.b + lowno * thisbase3) / 40;

                    newgrass1 = (highno * inputs.relief.highbase.r + lowno * inputs.relief.grass.r) / 40;
                    newgrass2 = (highno * inputs.relief.highbase.g + lowno * inputs.relief.grass.g) / 40;
                    newgrass3 = (highno * inputs.relief.highbase.b + lowno * inputs.relief.grass.b) / 40;
                }

                if (totalrain > 800.0f)
                {
                    colour1 = newbase1;
                    colour2 = newbase2;
                    colour3 = newbase3;
                }
                else if (totalrain > 200.0f)
                {
                    int wetno = (static_cast<int>(totalrain) - 200) / 40;

                    if (wetno > 20)
                        wetno = 20;

                    const int dryno = 20 - wetno;

                    colour1 = (wetno * newbase1 + dryno * newgrass1) / 20;
                    colour2 = (wetno * newbase2 + dryno * newgrass2) / 20;
                    colour3 = (wetno * newbase3 + dryno * newgrass3) / 20;
                }
                else
                {
                    float ftotalrain = 200.0f - totalrain;
                    ftotalrain = ftotalrain / 200.0f;

                    int powamount = static_cast<int>(totalrain) - 150;
                    if (powamount < 3)
                        powamount = 3;

                    ftotalrain = powf(ftotalrain, static_cast<float>(powamount));
                    ftotalrain = ftotalrain * 200.0f;
                    totalrain = 200.0f - ftotalrain;

                    const int wetno = static_cast<int>(totalrain);
                    const int dryno = 200 - wetno;

                    colour1 = (wetno * newgrass1 + dryno * newdesert1) / 200;
                    colour2 = (wetno * newgrass2 + dryno * newdesert2) / 200;
                    colour3 = (wetno * newgrass3 + dryno * newdesert3) / 200;
                }

                if (avetemp <= 0)
                {
                    assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                }
                else
                {
                    const RgbTriplet tundra = computetundracolour(inputs, y);

                    if (inputs.relief.snowchange == 1 && avetemp < 20)
                    {
                        if (avetemp < 6)
                        {
                            assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                        }
                        else if (avetemp < 10)
                        {
                            assignrgb(tundra, colour1, colour2, colour3);
                        }
                        else
                        {
                            const int hotno = avetemp - 10;
                            const int coldno = 10 - hotno;

                            colour1 = (hotno * colour1 + coldno * tundra.r) / 10;
                            colour2 = (hotno * colour2 + coldno * tundra.g) / 10;
                            colour3 = (hotno * colour3 + coldno * tundra.b) / 10;
                        }
                    }

                    if (inputs.relief.snowchange == 2 && avetemp < 20)
                    {
                        if (avetemp < 6)
                        {
                            assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                        }
                        else if (avetemp < 10)
                        {
                            if (deterministicrandom(inputs.relief.seed, x, y, 41, 6, 10) < avetemp)
                                assignrgb(tundra, colour1, colour2, colour3);
                            else
                                assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                        }
                        else
                        {
                            const int hotno = avetemp - 10;
                            const int coldno = 10 - hotno;

                            colour1 = (hotno * colour1 + coldno * tundra.r) / 10;
                            colour2 = (hotno * colour2 + coldno * tundra.g) / 10;
                            colour3 = (hotno * colour3 + coldno * tundra.b) / 10;
                        }
                    }

                    if (inputs.relief.snowchange == 3 && avetemp < 20)
                    {
                        if (avetemp < 10)
                        {
                            const int hotno = avetemp;
                            const int coldno = 10 - hotno;

                            colour1 = (hotno * tundra.r + coldno * inputs.relief.cold.r) / 10;
                            colour2 = (hotno * tundra.g + coldno * inputs.relief.cold.g) / 10;
                            colour3 = (hotno * tundra.b + coldno * inputs.relief.cold.b) / 10;
                        }
                        else
                        {
                            const int hotno = avetemp - 10;
                            const int coldno = 10 - hotno;

                            colour1 = (hotno * colour1 + coldno * tundra.r) / 10;
                            colour2 = (hotno * colour2 + coldno * tundra.g) / 10;
                            colour3 = (hotno * colour3 + coldno * tundra.b) / 10;
                        }
                    }
                }

                if (special == 120)
                {
                    colour1 = (colour1 * 2 + inputs.relief.erg.r) / 3;
                    colour2 = (colour2 * 2 + inputs.relief.erg.g) / 3;
                    colour3 = (colour3 * 2 + inputs.relief.erg.b) / 3;
                }

                if (special >= 130 && special < 140)
                {
                    colour1 = (colour1 * 2 + inputs.relief.wetlands.r) / 3;
                    colour2 = (colour2 * 2 + inputs.relief.wetlands.g) / 3;
                    colour3 = (colour3 * 2 + inputs.relief.wetlands.b) / 3;
                }
            }
        }

        if (sea)
        {
            const int amount = deterministicrandomsigned(inputs.relief.seed, x, y, 51, var);
            colour1 += amount;
            colour2 += amount;
            colour3 += amount;
        }
        else
        {
            colour1 += deterministicrandomsigned(inputs.relief.seed, x, y, 52, var);
            colour2 += deterministicrandomsigned(inputs.relief.seed, x, y, 53, var);
            colour3 += deterministicrandomsigned(inputs.relief.seed, x, y, 54, var);

            if (truelake || riftlake)
                assignrgb(inputs.relief.lake, colour1, colour2, colour3);
        }

        colour1 = clamp255(colour1);
        colour2 = clamp255(colour2);
        colour3 = clamp255(colour3);

        int r = colour1;
        int g = colour2;
        int b = colour3;

        if (!inputs.noshademap[index])
        {
            bool goahead = true;

            if (isseaicevisible(inputs, x, y))
                goahead = false;

            if (goahead || !sea)
            {
                int slope1 = 0;
                int slope2 = 0;
                bool hasslope1 = false;
                bool hasslope2 = false;

                if (inputs.relief.shadingdir == 2)
                {
                    slope1 = getslopecached(inputs, x - 1, y, mapvalue, hasslope1);
                    slope2 = getslopecached(inputs, x, y + 1, mapvalue, hasslope2);
                }

                if (inputs.relief.shadingdir == 4)
                {
                    slope1 = getslopecached(inputs, x - 1, y, mapvalue, hasslope1);
                    slope2 = getslopecached(inputs, x, y - 1, mapvalue, hasslope2);
                }

                if (inputs.relief.shadingdir == 6)
                {
                    slope1 = getslopecached(inputs, x + 1, y, mapvalue, hasslope1);
                    slope2 = getslopecached(inputs, x, y - 1, mapvalue, hasslope2);
                }

                if (inputs.relief.shadingdir == 8)
                {
                    slope1 = getslopecached(inputs, x + 1, y, mapvalue, hasslope1);
                    slope2 = getslopecached(inputs, x, y + 1, mapvalue, hasslope2);
                }

                if (hasslope1 && hasslope2)
                {
                    int totalslope = (slope1 + slope2) / 10;

                    if (totalslope > 40)
                        totalslope = 40;

                    if (totalslope < -40)
                        totalslope = -40;

                    if (sea)
                        totalslope = static_cast<int>(static_cast<float>(totalslope) * (inputs.relief.seashading * 2.0f));
                    else if (truelake)
                        totalslope = static_cast<int>(static_cast<float>(totalslope) * (inputs.relief.lakeshading * 2.0f));
                    else
                        totalslope = static_cast<int>(static_cast<float>(totalslope) * (inputs.relief.landshading * 2.0f));

                    if (mapvalue <= inputs.sealevel && inputs.oceanriftmap[index] == 0)
                    {
                        bool found = false;
                        bool ignore = false;

                        for (int k = x - 3; k <= x + 3 && !ignore; k++)
                        {
                            const int kk = wrapx(k, inputs.width);

                            for (int l = y - 3; l <= y + 3; l++)
                            {
                                if (l >= 0 && l <= inputs.height && inputs.oceanriftmap[inputindex(kk, l)] != 0)
                                {
                                    ignore = true;
                                    break;
                                }
                            }
                        }

                        if (!ignore)
                        {
                            for (int k = x - 1; k <= x + 1 && !found; k++)
                            {
                                const int kk = wrapx(k, inputs.width);

                                for (int l = y - 1; l <= y + 1; l++)
                                {
                                    if (l >= 0 && l <= inputs.height && inputs.oceanridgemap[inputindex(kk, l)] != 0)
                                    {
                                        found = true;
                                        break;
                                    }
                                }
                            }

                            if (found)
                                totalslope /= 4;
                        }
                    }

                    r += totalslope;
                    g += totalslope;
                    b += totalslope;
                }

                r = clamp255(r);
                g = clamp255(g);
                b = clamp255(b);
            }
        }

        writepixel(output, outputwidth, x, y, r, g, b);
    }

    __global__ void renderregionalelevationkernel(RegionalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        int heightpoint = static_cast<int>(inputs.map[index]);

        if (inputs.specials[index] > 100 && inputs.specials[index] < 130)
            heightpoint = inputs.lakemap[index];

        const int div = max(inputs.maxelevation / 255, 1);
        const int colour = clamp255(heightpoint / div);
        writepixel(output, outputwidth, localx, localy, colour, colour, colour);
    }

    __global__ void renderregionaltemperaturekernel(RegionalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;

        if (!regionaloutline(inputs, x, y))
        {
            int temperature = (regionalmintemp(inputs, x, y) + regionalmaxtemp(inputs, x, y)) / 2;
            temperature += 10;

            if (temperature > 0)
            {
                colour1 = 250;
                colour2 = 250 - (temperature * 3);
                colour3 = 250 - (temperature * 7);
            }
            else
            {
                temperature = abs(temperature);
                colour1 = 250 - (temperature * 7);
                colour2 = 250 - (temperature * 7);
                colour3 = 250;
            }
        }

        writepixel(output, outputwidth, localx, localy, clamp255(colour1), clamp255(colour2), clamp255(colour3));
    }

    __global__ void renderregionalprecipitationkernel(RegionalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;

        if (!regionaloutline(inputs, x, y))
        {
            const int rainfall = (regionalsummerrain(inputs, x, y) + regionalwinterrain(inputs, x, y)) / 8;
            colour1 = 255 - rainfall;
            colour2 = 255 - rainfall;
            colour3 = 255;
        }

        if (inputs.testmap[index] != 0)
        {
            colour1 = 255;
            colour2 = 0;
            colour3 = 255;
        }

        writepixel(output, outputwidth, localx, localy, clamp255(colour1), clamp255(colour2), clamp255(colour3));
    }

    __global__ void renderregionalclimatekernel(RegionalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;

        if (regionalsea(inputs, x, y))
        {
            const int seaice = static_cast<int>(inputs.seaicemap[index]);

            if (seaice == 2)
            {
                colour1 = 243;
                colour2 = 243;
                colour3 = 255;
            }
            else if (seaice == 1)
            {
                colour1 = 228;
                colour2 = 228;
                colour3 = 255;
            }
            else
            {
                colour1 = 13;
                colour2 = 49;
                colour3 = 109;
            }
        }
        else
        {
            climatecolour(inputs.climatemap[index], colour1, colour2, colour3);
        }

        writepixel(output, outputwidth, localx, localy, colour1, colour2, colour3);
    }

    __global__ void renderregionalriverskernel(RegionalCudaDeviceInputs inputs, std::uint8_t* output, int outputwidth)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        int colour1 = 255;
        int colour2 = 255;
        int colour3 = 255;

        if (regionaloutline(inputs, x, y))
        {
            colour1 = 0;
            colour2 = 0;
            colour3 = 0;
        }
        else
        {
            int flow = regionalriveraveflow(inputs, x, y);
            const int mult = max(inputs.maxriverflow / 400, 1);

            if (flow > 0)
            {
                flow *= 100;
                colour1 = clamp255(255 - (flow / mult));
                colour2 = colour1;
            }

            colour3 = 255;

            if (regionaltruelake(inputs, x, y))
            {
                colour1 = 150;
                colour2 = 150;
                colour3 = 250;
            }

            if (inputs.volcanomap[index])
            {
                colour1 = 240;
                colour2 = 0;
                colour3 = 0;
            }

            if (colour1 == 255 && colour2 == 255 && colour3 == 255)
            {
                const int special = static_cast<int>(inputs.specials[index]);

                if (special >= 130 && special < 140)
                {
                    colour1 = 30;
                    colour2 = 250;
                    colour3 = 150;
                }

                if (inputs.mudmap[index])
                {
                    colour1 = 131;
                    colour2 = 98;
                    colour3 = 75;
                }

                if (inputs.sandmap[index] || inputs.shinglemap[index])
                {
                    colour1 = 255;
                    colour2 = 255;
                    colour3 = 50;
                }
            }
        }

        writepixel(output, outputwidth, localx, localy, colour1, colour2, colour3);
    }

    __global__ void renderregionalreliefbasekernel(RegionalCudaDeviceInputs inputs, short* relief1, short* relief2, short* relief3)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputwidth = inputs.regwidthend - inputs.regwidthbegin + 1;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        const int yy = inputs.lefty + (y / 16);
        const int mapvalue = static_cast<int>(inputs.map[index]);
        const int special = static_cast<int>(inputs.specials[index]);
        const int riverflow = regionalriveraveflow(inputs, x, y);
        const int fakeflow = regionalfakeaveflow(inputs, x, y);
        const bool sea = regionalsea(inputs, x, y);
        const bool truelake = regionaltruelake(inputs, x, y);
        const bool mangrove = regionalmangrove(inputs, x, y);

        int avetemp = (regionalextramaxtemp(inputs, x, y) + regionalextramintemp(inputs, x, y)) / 2;
        avetemp += 1000;

        int totalrain = regionalsummerrain(inputs, x, y) + regionalwinterrain(inputs, x, y);

        if (special == 120)
            totalrain /= 4;

        int colour1 = 0;
        int colour2 = 0;
        int colour3 = 0;
        int var = 10;
        bool stripe = true;

        if (sea)
        {
            if (regionalseaicevisible(inputs, x, y))
            {
                assignrgb(inputs.relief.seaice, colour1, colour2, colour3);
                stripe = false;
                var = 0;
            }
            else
            {
                colour1 = (inputs.relief.ocean.r * mapvalue + inputs.relief.deepocean.r * (inputs.sealevel - mapvalue)) / inputs.sealevel;
                colour2 = (inputs.relief.ocean.g * mapvalue + inputs.relief.deepocean.g * (inputs.sealevel - mapvalue)) / inputs.sealevel;
                colour3 = (inputs.relief.ocean.b * mapvalue + inputs.relief.deepocean.b * (inputs.sealevel - mapvalue)) / inputs.sealevel;
                var = 5;
            }
        }
        else if (riverflow >= inputs.relief.minriverflow || fakeflow >= inputs.relief.minriverflow)
        {
            assignrgb(inputs.relief.river, colour1, colour2, colour3);
            stripe = false;
        }
        else if (special == 110)
        {
            assignrgb(inputs.relief.saltpan, colour1, colour2, colour3);
            var = 20;
            stripe = false;
        }
        else if (inputs.sandmap[index] || inputs.shinglemap[index])
        {
            if (inputs.sandmap[index] && !inputs.shinglemap[index])
            {
                assignrgb(inputs.relief.sand, colour1, colour2, colour3);
                var = 5;
            }

            if (inputs.shinglemap[index] && !inputs.sandmap[index])
            {
                assignrgb(inputs.relief.shingle, colour1, colour2, colour3);
                var = 40;
            }

            if (inputs.sandmap[index] && inputs.shinglemap[index])
            {
                colour1 = (inputs.relief.sand.r + inputs.relief.shingle.r) / 2;
                colour2 = (inputs.relief.sand.g + inputs.relief.shingle.g) / 2;
                colour3 = (inputs.relief.sand.b + inputs.relief.shingle.b) / 2;
                var = 20;
            }

            if (inputs.mudmap[index])
            {
                colour1 = (colour1 + inputs.relief.mud.r) / 2;
                colour2 = (colour2 + inputs.relief.mud.g) / 2;
                colour3 = (colour3 + inputs.relief.mud.b) / 2;
            }
        }
        else if (mangrove)
        {
            assignrgb(inputs.relief.mangrove, colour1, colour2, colour3);
            var = 20;
        }
        else
        {
            int thisbase1 = 0;
            int thisbase2 = 0;
            int thisbase3 = 0;
            int newdesert1 = 0;
            int newdesert2 = 0;
            int newdesert3 = 0;

            if (avetemp > 3000)
            {
                assignrgb(inputs.relief.base, thisbase1, thisbase2, thisbase3);
            }
            else
            {
                const int hotno = avetemp / 3;
                const int coldno = 1000 - hotno;
                thisbase1 = (hotno * inputs.relief.base.r + coldno * inputs.relief.basetemp.r) / 1000;
                thisbase2 = (hotno * inputs.relief.base.g + coldno * inputs.relief.basetemp.g) / 1000;
                thisbase3 = (hotno * inputs.relief.base.b + coldno * inputs.relief.basetemp.b) / 1000;
            }

            if (avetemp > 30)
            {
                assignrgb(inputs.relief.desert, newdesert1, newdesert2, newdesert3);
            }
            else if (avetemp <= 10)
            {
                assignrgb(inputs.relief.colddesert, newdesert1, newdesert2, newdesert3);
            }
            else
            {
                const int hotno = avetemp - 10;
                const int coldno = 20 - hotno;
                newdesert1 = (hotno * inputs.relief.desert.r + coldno * inputs.relief.colddesert.r) / 20;
                newdesert2 = (hotno * inputs.relief.desert.g + coldno * inputs.relief.colddesert.g) / 20;
                newdesert3 = (hotno * inputs.relief.desert.b + coldno * inputs.relief.colddesert.b) / 20;
            }

            float rainforestmult = static_cast<float>(regionalmintemp(inputs, x, y)) / 18.0f;
            rainforestmult = rainforestmult * static_cast<float>(regionalwinterrain(inputs, x, y)) / 80.0f;

            if (rainforestmult < 1.0f)
                rainforestmult = 1.0f;

            totalrain = totalrain * static_cast<int>(rainforestmult);

            int mapelev = mapvalue - inputs.sealevel;

            if (special == 110 || special == 120)
                mapelev = inputs.lakemap[index] - inputs.sealevel;

            const int desertmapelev = mapelev;

            if (inputs.relief.colourcliffs)
            {
                int biggestslope = 0;

                for (int k = x - 1; k <= x + 1; k++)
                {
                    if (k >= 0 && k <= inputs.regionwidth)
                    {
                        for (int l = y - 1; l <= y + 1; l++)
                        {
                            if (l >= 0 && l <= inputs.regionheight)
                            {
                                const int thisslope = mapelev + inputs.sealevel - static_cast<int>(inputs.map[regionalinputindex(k, l)]);

                                if (thisslope > biggestslope)
                                    biggestslope = thisslope;
                            }
                        }
                    }
                }

                biggestslope -= 240;

                if (biggestslope < 0)
                    biggestslope = 0;

                float adjustedelev = static_cast<float>(mapelev);
                adjustedelev = adjustedelev * (static_cast<float>(biggestslope) / 240.0f);

                if (adjustedelev > static_cast<float>(mapelev))
                    adjustedelev = static_cast<float>(mapelev);

                mapelev = static_cast<int>(adjustedelev);
            }

            int newbase1 = 0;
            int newbase2 = 0;
            int newbase3 = 0;
            int newgrass1 = 0;
            int newgrass2 = 0;
            int newgrass3 = 0;

            if (desertmapelev > 2000)
            {
                assignrgb(inputs.relief.highdesert, newdesert1, newdesert2, newdesert3);
            }
            else
            {
                const int highno = desertmapelev / 50;
                const int lowno = 40 - highno;
                newdesert1 = (highno * inputs.relief.highdesert.r + lowno * newdesert1) / 40;
                newdesert2 = (highno * inputs.relief.highdesert.g + lowno * newdesert2) / 40;
                newdesert3 = (highno * inputs.relief.highdesert.b + lowno * newdesert3) / 40;
            }

            if (mapelev > 3000)
            {
                assignrgb(inputs.relief.highbase, newbase1, newbase2, newbase3);
                assignrgb(inputs.relief.highbase, newgrass1, newgrass2, newgrass3);
            }
            else
            {
                const int highno = mapelev / 75;
                const int lowno = 40 - highno;
                newbase1 = (highno * inputs.relief.highbase.r + lowno * thisbase1) / 40;
                newbase2 = (highno * inputs.relief.highbase.g + lowno * thisbase2) / 40;
                newbase3 = (highno * inputs.relief.highbase.b + lowno * thisbase3) / 40;
                newgrass1 = (highno * inputs.relief.highbase.r + lowno * inputs.relief.grass.r) / 40;
                newgrass2 = (highno * inputs.relief.highbase.g + lowno * inputs.relief.grass.g) / 40;
                newgrass3 = (highno * inputs.relief.highbase.b + lowno * inputs.relief.grass.b) / 40;
            }

            if (regionalrivervalley(inputs, x, y) || inputs.rivermapdir[index] != 0 || inputs.fakeriversdir[index] != 0)
            {
                float biggestflow = 0.0f;

                for (int k = x - 20; k <= x + 20; k++)
                {
                    for (int l = y - 20; l <= y + 20; l++)
                    {
                        if (k >= 0 && k <= inputs.regionwidth && l >= 0 && l <= inputs.regionheight)
                        {
                            const float thisflow = static_cast<float>(inputs.rivermapjan[regionalinputindex(k, l)] + inputs.rivermapjul[regionalinputindex(k, l)]);
                            if (thisflow > biggestflow)
                                biggestflow = thisflow;
                        }
                    }
                }

                if (biggestflow > 12000.0f)
                    biggestflow = 1200.0f;

                float mult = static_cast<float>(totalrain);

                if (mult < 1.0f)
                    mult = 1.0f;

                if (mult > 1000.0f)
                    mult = 1000.0f;

                biggestflow = biggestflow / mult;
                totalrain += static_cast<int>(biggestflow);
            }

            if (totalrain > 800)
            {
                colour1 = newbase1;
                colour2 = newbase2;
                colour3 = newbase3;
            }
            else if (totalrain > 200)
            {
                int wetno = (totalrain - 200) / 40;

                if (wetno > 20)
                    wetno = 20;

                const int dryno = 20 - wetno;
                colour1 = (wetno * newbase1 + dryno * newgrass1) / 20;
                colour2 = (wetno * newbase2 + dryno * newgrass2) / 20;
                colour3 = (wetno * newbase3 + dryno * newgrass3) / 20;
            }
            else
            {
                float ftotalrain = 200.0f - static_cast<float>(totalrain);
                ftotalrain = ftotalrain / 200.0f;

                int powamount = totalrain - 150;
                if (powamount < 3)
                    powamount = 3;

                ftotalrain = powf(ftotalrain, static_cast<float>(powamount));
                ftotalrain = ftotalrain * 200.0f;
                totalrain = 200 - static_cast<int>(ftotalrain);

                const int wetno = totalrain;
                const int dryno = 200 - wetno;
                colour1 = (wetno * newgrass1 + dryno * newdesert1) / 200;
                colour2 = (wetno * newgrass2 + dryno * newdesert2) / 200;
                colour3 = (wetno * newgrass3 + dryno * newdesert3) / 200;
            }

            if (avetemp <= 0 || yy > inputs.globalheight - 3)
            {
                assignrgb(inputs.relief.cold, colour1, colour2, colour3);
            }
            else
            {
                const RgbTriplet tundra = computetundracolourforrow(inputs.globalheight, inputs.relief, yy);

                if (inputs.relief.snowchange == 1 && avetemp < 2000)
                {
                    if (avetemp < 600)
                    {
                        assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                    }
                    else if (avetemp < 1000)
                    {
                        assignrgb(tundra, colour1, colour2, colour3);
                    }
                    else
                    {
                        const int hotno = avetemp - 1000;
                        const int coldno = 1000 - hotno;
                        colour1 = (hotno * colour1 + coldno * tundra.r) / 1000;
                        colour2 = (hotno * colour2 + coldno * tundra.g) / 1000;
                        colour3 = (hotno * colour3 + coldno * tundra.b) / 1000;
                    }
                }

                if (inputs.relief.snowchange == 2 && avetemp < 2000)
                {
                    if (avetemp < 600)
                    {
                        assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                    }
                    else if (avetemp < 1000)
                    {
                        if (deterministicrandom(inputs.relief.seed, x, y, 201, 600, 1000) < avetemp)
                            assignrgb(tundra, colour1, colour2, colour3);
                        else
                            assignrgb(inputs.relief.cold, colour1, colour2, colour3);
                    }
                    else
                    {
                        const int hotno = avetemp - 1000;
                        const int coldno = 1000 - hotno;
                        colour1 = (hotno * colour1 + coldno * tundra.r) / 1000;
                        colour2 = (hotno * colour2 + coldno * tundra.g) / 1000;
                        colour3 = (hotno * colour3 + coldno * tundra.b) / 1000;
                    }
                }

                if (inputs.relief.snowchange == 3 && avetemp < 2000)
                {
                    if (avetemp < 1000)
                    {
                        const int hotno = avetemp;
                        const int coldno = 1000 - hotno;
                        colour1 = (hotno * tundra.r + coldno * inputs.relief.cold.r) / 1000;
                        colour2 = (hotno * tundra.g + coldno * inputs.relief.cold.g) / 1000;
                        colour3 = (hotno * tundra.b + coldno * inputs.relief.cold.b) / 1000;
                    }
                    else
                    {
                        const int hotno = avetemp - 1000;
                        const int coldno = 1000 - hotno;
                        colour1 = (hotno * colour1 + coldno * tundra.r) / 1000;
                        colour2 = (hotno * colour2 + coldno * tundra.g) / 1000;
                        colour3 = (hotno * colour3 + coldno * tundra.b) / 1000;
                    }
                }
            }

            if (special == 120)
            {
                colour1 = (colour1 * 6 + inputs.relief.erg.r) / 7;
                colour2 = (colour2 * 6 + inputs.relief.erg.g) / 7;
                colour3 = (colour3 * 6 + inputs.relief.erg.b) / 7;
                var = 10;
            }

            if (inputs.mudmap[index])
            {
                colour1 = (colour1 + inputs.relief.mud.r * 2) / 3;
                colour2 = (colour2 + inputs.relief.mud.g * 2) / 3;
                colour3 = (colour3 + inputs.relief.mud.b * 2) / 3;
                var = 10;
            }

            if (special >= 130 && special < 140)
            {
                colour1 = (colour1 * 2 + inputs.relief.wetlands.r) / 3;
                colour2 = (colour2 * 2 + inputs.relief.wetlands.g) / 3;
                colour3 = (colour3 * 2 + inputs.relief.wetlands.b) / 3;
            }
        }

        if (sea)
        {
            const int amount = deterministicrandomsigned(inputs.relief.seed, x, y, 210, var);
            colour1 += amount;
            colour2 += amount;
            colour3 += amount;
        }
        else
        {
            colour1 += deterministicrandomsigned(inputs.relief.seed, x, y, 211, var);
            colour2 += deterministicrandomsigned(inputs.relief.seed, x, y, 212, var);
            colour3 += deterministicrandomsigned(inputs.relief.seed, x, y, 213, var);

            if (truelake)
            {
                assignrgb(inputs.relief.lake, colour1, colour2, colour3);

                for (int k = x - 1; k <= x + 1; k++)
                {
                    for (int l = y - 1; l <= y + 1; l++)
                    {
                        if (k >= 0 && k <= inputs.regionwidth && l >= 0 && l <= inputs.regionheight)
                        {
                            const int neighbour = regionalinputindex(k, l);
                            if ((inputs.rivermapdir[neighbour] != 0 || inputs.fakeriversdir[neighbour] > 0) && inputs.lakemap[neighbour] == 0)
                            {
                                colour1 = (inputs.relief.lake.r + inputs.relief.river.r) / 2;
                                colour2 = (inputs.relief.lake.g + inputs.relief.river.g) / 2;
                                colour3 = (inputs.relief.lake.b + inputs.relief.river.b) / 2;
                                k = x + 1;
                                l = y + 1;
                            }
                        }
                    }
                }
            }
        }

        if (special == 140)
        {
            assignrgb(inputs.relief.glacier, colour1, colour2, colour3);
            stripe = false;
        }

        if (stripe)
        {
            float stripevar = static_cast<float>(avetemp);

            if (stripevar > 5.0f)
                stripevar = 5.0f;

            if (stripevar < 0.0f)
                stripevar = 0.0f;

            if (special > 100)
                stripevar = 1.0f;

            if (sea)
                stripevar = stripevar * inputs.relief.seamarbling;
            else if (truelake)
                stripevar = stripevar * inputs.relief.lakemarbling;
            else
                stripevar = stripevar * inputs.relief.landmarbling;

            const std::uint64_t stripeseed = inputs.relief.seed ^ static_cast<std::uint64_t>(static_cast<long long>(regionalmaxtemp(inputs, x, y) - 5000));
            const int stripeamount = static_cast<int>(stripevar);
            colour1 += deterministicrandomsigned(stripeseed, x, y, 214, stripeamount);
            colour2 += deterministicrandomsigned(stripeseed, x, y, 215, stripeamount);
            colour2 += deterministicrandomsigned(stripeseed, x, y, 216, stripeamount);
        }

        colour1 = clamp255(colour1);
        colour2 = clamp255(colour2);
        colour3 = clamp255(colour3);

        writeshortpixel(relief1, x, y, static_cast<short>(colour1));
        writeshortpixel(relief2, x, y, static_cast<short>(colour2));
        writeshortpixel(relief3, x, y, static_cast<short>(colour3));
    }

    __global__ void blurregionalwetlandskernel(RegionalCudaDeviceInputs inputs, const short* input1, const short* input2, const short* input3, short* output1, short* output2, short* output3)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputwidth = inputs.regwidthend - inputs.regwidthbegin + 1;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        short r = input1[index];
        short g = input2[index];
        short b = input3[index];
        const int special = static_cast<int>(inputs.specials[index]);

        if (x > inputs.regwidthbegin && x < inputs.regwidthend && y > inputs.regheightbegin && y < inputs.regheightend)
        {
            const bool mangrove = regionalmangrove(inputs, x, y);
            if ((special == 130 || special == 131 || special == 132 || inputs.mudmap[index]) && !mangrove)
            {
                if (regionalriveraveflow(inputs, x, y) < inputs.relief.minriverflow || regionalfakeaveflow(inputs, x, y) < inputs.relief.minriverflow)
                {
                    float colred = 0.0f;
                    float colgreen = 0.0f;
                    float colblue = 0.0f;

                    for (int k = x - 1; k <= x + 1; k++)
                    {
                        for (int l = y - 1; l <= y + 1; l++)
                        {
                            colred += static_cast<float>(readshortpixel(input1, k, l));
                            colgreen += static_cast<float>(readshortpixel(input2, k, l));
                            colblue += static_cast<float>(readshortpixel(input3, k, l));
                        }
                    }

                    r = static_cast<short>(colred / 9.0f);
                    g = static_cast<short>(colgreen / 9.0f);
                    b = static_cast<short>(colblue / 9.0f);
                }
            }
        }

        writeshortpixel(output1, x, y, r);
        writeshortpixel(output2, x, y, g);
        writeshortpixel(output3, x, y, b);
    }

    __global__ void blurregionalvalleyskernel(RegionalCudaDeviceInputs inputs, const short* input1, const short* input2, const short* input3, short* output1, short* output2, short* output3)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputwidth = inputs.regwidthend - inputs.regwidthbegin + 1;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        short r = input1[index];
        short g = input2[index];
        short b = input3[index];

        if (x > inputs.regwidthbegin && x < inputs.regwidthend && y > inputs.regheightbegin && y < inputs.regheightend)
        {
            if (regionalrivervalley(inputs, x, y) && inputs.specials[index] < 130)
            {
                if (!(regionalriveraveflow(inputs, x, y) >= inputs.relief.minriverflow || regionalfakeaveflow(inputs, x, y) >= inputs.relief.minriverflow))
                {
                    float colred = 0.0f;
                    float colgreen = 0.0f;
                    float colblue = 0.0f;
                    float count = 0.0f;

                    for (int k = x - 1; k <= x + 1; k++)
                    {
                        for (int l = y - 1; l <= y + 1; l++)
                        {
                            const int neighbour = regionalinputindex(k, l);
                            if (inputs.rivermapjan[neighbour] == 0 && inputs.rivermapjul[neighbour] == 0
                                && inputs.fakeriversjan[neighbour] == 0 && inputs.fakeriversjul[neighbour] == 0
                                && inputs.deltamapjan[neighbour] == 0 && inputs.deltamapjul[neighbour] == 0)
                            {
                                colred += static_cast<float>(readshortpixel(input1, k, l));
                                colgreen += static_cast<float>(readshortpixel(input2, k, l));
                                colblue += static_cast<float>(readshortpixel(input3, k, l));
                                count += 1.0f;
                            }
                        }
                    }

                    if (count > 0.0f)
                    {
                        r = static_cast<short>(colred / count);
                        g = static_cast<short>(colgreen / count);
                        b = static_cast<short>(colblue / count);
                    }
                }
            }
        }

        writeshortpixel(output1, x, y, r);
        writeshortpixel(output2, x, y, g);
        writeshortpixel(output3, x, y, b);
    }

    __global__ void redrawregionalreliefriverskernel(RegionalCudaDeviceInputs inputs, short* relief1, short* relief2, short* relief3)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputwidth = inputs.regwidthend - inputs.regwidthbegin + 1;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);

        if (!regionalsea(inputs, x, y) && inputs.specials[index] < 130 && !regionaltruelake(inputs, x, y))
        {
            if (regionalriveraveflow(inputs, x, y) >= inputs.relief.minriverflow || regionalfakeaveflow(inputs, x, y) >= inputs.relief.minriverflow)
            {
                int colour1 = inputs.relief.river.r;
                int colour2 = inputs.relief.river.g;
                int colour3 = inputs.relief.river.b;
                const std::uint64_t stripeseed = inputs.relief.seed ^ static_cast<std::uint64_t>(static_cast<long long>(regionalmaxtemp(inputs, x, y) + regionalwinterrain(inputs, x, y)));

                colour1 += deterministicrandomsigned(stripeseed, x, y, 220, 5);
                colour2 += deterministicrandomsigned(stripeseed, x, y, 221, 5);
                colour2 += deterministicrandomsigned(stripeseed, x, y, 222, 5);

                writeshortpixel(relief1, x, y, static_cast<short>(clamp255(colour1)));
                writeshortpixel(relief2, x, y, static_cast<short>(clamp255(colour2)));
                writeshortpixel(relief3, x, y, static_cast<short>(clamp255(colour3)));
            }
        }
    }

    __global__ void shaderegionalreliefkernel(RegionalCudaDeviceInputs inputs, const short* relief1, const short* relief2, const short* relief3, std::uint8_t* output, int outputwidth)
    {
        const int localx = blockIdx.x * blockDim.x + threadIdx.x;
        const int localy = blockIdx.y * blockDim.y + threadIdx.y;
        const int outputheight = inputs.regheightend - inputs.regheightbegin + 1;

        if (localx >= outputwidth || localy >= outputheight)
            return;

        const int x = inputs.regwidthbegin + localx;
        const int y = inputs.regheightbegin + localy;
        const int index = regionalinputindex(x, y);
        int r = static_cast<int>(relief1[index]);
        int g = static_cast<int>(relief2[index]);
        int b = static_cast<int>(relief3[index]);
        const int mapvalue = static_cast<int>(inputs.map[index]);
        const bool sea = regionalsea(inputs, x, y);
        const bool truelake = regionaltruelake(inputs, x, y);

        if (inputs.specials[index] == 0 && inputs.rivermapdir[index] == 0 && inputs.fakeriversdir[index] == 0)
        {
            bool goahead = true;

            if (regionalseaicevisible(inputs, x, y))
                goahead = false;

            if (goahead || !sea)
            {
                int slope1 = 0;
                int slope2 = 0;
                bool hasslope1 = false;
                bool hasslope2 = false;

                if (inputs.relief.shadingdir == 2)
                {
                    if (x - 1 >= 0 && x - 1 <= inputs.regionwidth)
                    {
                        slope1 = static_cast<int>(inputs.map[regionalinputindex(x - 1, y)]) - mapvalue;
                        hasslope1 = true;
                    }

                    if (y + 1 >= 0 && y + 1 <= inputs.regionheight)
                    {
                        slope2 = static_cast<int>(inputs.map[regionalinputindex(x, y + 1)]) - mapvalue;
                        hasslope2 = true;
                    }
                }

                if (inputs.relief.shadingdir == 4)
                {
                    if (x - 1 >= 0 && x - 1 <= inputs.regionwidth)
                    {
                        slope1 = static_cast<int>(inputs.map[regionalinputindex(x - 1, y)]) - mapvalue;
                        hasslope1 = true;
                    }

                    if (y - 1 >= 0 && y - 1 <= inputs.regionheight)
                    {
                        slope2 = static_cast<int>(inputs.map[regionalinputindex(x, y - 1)]) - mapvalue;
                        hasslope2 = true;
                    }
                }

                if (inputs.relief.shadingdir == 6)
                {
                    if (x + 1 >= 0 && x + 1 <= inputs.regionwidth)
                    {
                        slope1 = static_cast<int>(inputs.map[regionalinputindex(x + 1, y)]) - mapvalue;
                        hasslope1 = true;
                    }

                    if (y - 1 >= 0 && y - 1 <= inputs.regionheight)
                    {
                        slope2 = static_cast<int>(inputs.map[regionalinputindex(x, y - 1)]) - mapvalue;
                        hasslope2 = true;
                    }
                }

                if (inputs.relief.shadingdir == 8)
                {
                    if (x + 1 >= 0 && x + 1 <= inputs.regionwidth)
                    {
                        slope1 = static_cast<int>(inputs.map[regionalinputindex(x + 1, y)]) - mapvalue;
                        hasslope1 = true;
                    }

                    if (y + 1 >= 0 && y + 1 <= inputs.regionheight)
                    {
                        slope2 = static_cast<int>(inputs.map[regionalinputindex(x, y + 1)]) - mapvalue;
                        hasslope2 = true;
                    }
                }

                if (hasslope1 && hasslope2)
                {
                    int totalslope = (slope1 + slope2) / 10;

                    if (totalslope > 40)
                        totalslope = 40;

                    if (totalslope < -40)
                        totalslope = -40;

                    float thisshading = inputs.relief.landshading;

                    if (truelake)
                        thisshading = inputs.relief.lakeshading;

                    if (sea)
                        thisshading = inputs.relief.seashading;

                    totalslope = static_cast<int>(static_cast<float>(totalslope) * (thisshading * 2.0f));
                    r += totalslope;
                    g += totalslope;
                    b += totalslope;
                }
            }
        }

        writepixel(output, outputwidth, localx, localy, clamp255(r), clamp255(g), clamp255(b));
    }

    __global__ void nearestdisplaykernel(const std::uint8_t* sourcepixels, int sourcewidth, int worldwidth, std::uint8_t* displaypixels)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x >= displaywidth || y >= displayheight)
            return;

        const float mapdiv = (static_cast<float>(worldwidth) + 1.0f) / static_cast<float>(displaywidth);
        const int sourcex = static_cast<int>(static_cast<float>(x) * mapdiv);
        const int sourcey = static_cast<int>(static_cast<float>(y) * mapdiv);
        const std::uint8_t* sourcepixel = readpixel(sourcepixels, sourcewidth, sourcex, sourcey);
        writepixel(displaypixels, displaywidth, x, y, sourcepixel[0], sourcepixel[1], sourcepixel[2], sourcepixel[3]);
    }

    __global__ void outlineddisplaykernel(const std::uint8_t* sourcepixels, int sourcewidth, int worldwidth, std::uint8_t* displaypixels)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x >= displaywidth || y >= displayheight)
            return;

        const float mapdiv = (static_cast<float>(worldwidth) + 1.0f) / static_cast<float>(displaywidth);
        const int sourcex = static_cast<int>(static_cast<float>(x) * mapdiv);
        const int sourcey = static_cast<int>(static_cast<float>(y) * mapdiv);
        const std::uint8_t* sourcepixel = readpixel(sourcepixels, sourcewidth, sourcex, sourcey);

        int r = sourcepixel[0];
        int g = sourcepixel[1];
        int b = sourcepixel[2];
        int a = sourcepixel[3];

        if (worldwidth > displaywidth && sourcex > 0 && sourcey > 0)
        {
            for (int k = sourcex - 1; k <= sourcex; k++)
            {
                for (int l = sourcey - 1; l <= sourcey; l++)
                {
                    const std::uint8_t* nearby = readpixel(sourcepixels, sourcewidth, k, l);

                    if (nearby[0] == 0 && nearby[1] == 0 && nearby[2] == 0)
                    {
                        r = nearby[0];
                        g = nearby[1];
                        b = nearby[2];
                        a = nearby[3];
                    }
                }
            }
        }

        writepixel(displaypixels, displaywidth, x, y, r, g, b, a);
    }

    __global__ void riversdisplaykernel(const std::uint8_t* sourcepixels, int sourcewidth, int worldwidth, int worldsize, std::uint8_t* displaypixels)
    {
        const int x = blockIdx.x * blockDim.x + threadIdx.x;
        const int y = blockIdx.y * blockDim.y + threadIdx.y;

        if (x >= displaywidth || y >= displayheight)
            return;

        int sourcex = x;
        int sourcey = y;

        if (worldsize == 0)
        {
            sourcex = x / 2;
            sourcey = y / 2;
        }
        else if (worldsize == 2)
        {
            const float mapdiv = (static_cast<float>(worldwidth) + 1.0f) / static_cast<float>(displaywidth);
            sourcex = static_cast<int>(static_cast<float>(x) * mapdiv);
            sourcey = static_cast<int>(static_cast<float>(y) * mapdiv);
        }

        const std::uint8_t* sourcepixel = readpixel(sourcepixels, sourcewidth, sourcex, sourcey);
        int r = sourcepixel[0];
        int g = sourcepixel[1];
        int b = sourcepixel[2];
        int a = sourcepixel[3];

        if (worldsize == 2 && sourcex > 0 && sourcey > 0)
        {
            int lowestred = 255;
            int lowestgreen = 255;
            int lowestblue = 255;

            for (int k = sourcex - 1; k <= sourcex; k++)
            {
                for (int l = sourcey - 1; l <= sourcey; l++)
                {
                    const std::uint8_t* nearby = readpixel(sourcepixels, sourcewidth, k, l);

                    if (nearby[0] < lowestred && nearby[0] == nearby[1])
                    {
                        lowestred = nearby[0];
                        lowestgreen = nearby[1];
                        lowestblue = nearby[2];
                    }
                }
            }

            if (lowestred < 255)
            {
                r = lowestred;
                g = lowestgreen;
                b = lowestblue;
            }

            if (worldwidth > displaywidth)
            {
                for (int k = sourcex - 1; k <= sourcex; k++)
                {
                    for (int l = sourcey - 1; l <= sourcey; l++)
                    {
                        const std::uint8_t* nearby = readpixel(sourcepixels, sourcewidth, k, l);

                        if (nearby[0] == 0 && nearby[1] == 0 && nearby[2] == 0)
                        {
                            r = nearby[0];
                            g = nearby[1];
                            b = nearby[2];
                            a = nearby[3];
                        }
                    }
                }
            }
        }

        writepixel(displaypixels, displaywidth, x, y, r, g, b, a);
    }

    template <typename Kernel, typename DisplayKernel>
    bool renderandcopy(Kernel kernel, DisplayKernel displaykernel, const GlobalCudaDeviceInputs& inputs, DeviceBuffer<std::uint8_t>& fullpixels, DeviceBuffer<std::uint8_t>& displaypixels, std::vector<std::uint8_t>& fulloutput, std::vector<std::uint8_t>& displayoutput)
    {
        const int imagewidth = inputs.width + 1;
        const int imageheight = inputs.height + 2;
        const size_t fullbytes = static_cast<size_t>(imagewidth) * static_cast<size_t>(imageheight) * rgbachannels;
        const size_t displaybytes = static_cast<size_t>(displaywidth) * static_cast<size_t>(displayheight) * rgbachannels;

        if (!fullpixels.ensure(fullbytes) || !displaypixels.ensure(displaybytes))
            return false;

        if (!cudacheck(cudaMemset(fullpixels.data(), 0, fullbytes)))
            return false;

        const dim3 block(16, 16);
        const dim3 grid((inputs.width + 16) / 16, (inputs.height + 16) / 16);
        kernel<<<grid, block>>>(inputs, fullpixels.data(), imagewidth);

        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        const dim3 displaygrid((displaywidth + 15) / 16, (displayheight + 15) / 16);
        displaykernel<<<displaygrid, block>>>(fullpixels.data(), imagewidth, inputs.width, displaypixels.data());

        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        return copyoutput(fulloutput, fullpixels, fullbytes) && copyoutput(displayoutput, displaypixels, displaybytes);
    }

    bool renderrivers(const GlobalCudaDeviceInputs& inputs, DeviceBuffer<std::uint8_t>& fullpixels, DeviceBuffer<std::uint8_t>& displaypixels, std::vector<std::uint8_t>& fulloutput, std::vector<std::uint8_t>& displayoutput)
    {
        const int imagewidth = inputs.width + 1;
        const int imageheight = inputs.height + 2;
        const size_t fullbytes = static_cast<size_t>(imagewidth) * static_cast<size_t>(imageheight) * rgbachannels;
        const size_t displaybytes = static_cast<size_t>(displaywidth) * static_cast<size_t>(displayheight) * rgbachannels;

        if (!fullpixels.ensure(fullbytes) || !displaypixels.ensure(displaybytes))
            return false;

        if (!cudacheck(cudaMemset(fullpixels.data(), 0, fullbytes)))
            return false;

        const dim3 block(16, 16);
        const dim3 grid((inputs.width + 16) / 16, (inputs.height + 16) / 16);
        renderriverskernel<<<grid, block>>>(inputs, fullpixels.data(), imagewidth);

        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        const dim3 displaygrid((displaywidth + 15) / 16, (displayheight + 15) / 16);
        riversdisplaykernel<<<displaygrid, block>>>(fullpixels.data(), imagewidth, inputs.width, inputs.worldsize, displaypixels.data());

        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        return copyoutput(fulloutput, fullpixels, fullbytes) && copyoutput(displayoutput, displaypixels, displaybytes);
    }

    struct GlobalCudaDeviceBuffers
    {
        DeviceBuffer<short> mapnom;
        DeviceBuffer<int> lakemap;
        DeviceBuffer<short> oceanridgeheightmap;
        DeviceBuffer<short> mountainheights;
        DeviceBuffer<short> volcanomap;
        DeviceBuffer<short> extraelevmap;
        DeviceBuffer<short> craterrims;
        DeviceBuffer<short> cratercentres;
        DeviceBuffer<short> oceanridgemap;
        DeviceBuffer<short> oceanriftmap;
        DeviceBuffer<short> jantempmap;
        DeviceBuffer<short> jultempmap;
        DeviceBuffer<short> janrainmap;
        DeviceBuffer<short> julrainmap;
        DeviceBuffer<short> seaicemap;
        DeviceBuffer<short> climatemap;
        DeviceBuffer<int> rivermapjan;
        DeviceBuffer<int> rivermapjul;
        DeviceBuffer<short> specials;
        DeviceBuffer<short> deltamapdir;
        DeviceBuffer<int> deltamapjan;
        DeviceBuffer<int> deltamapjul;
        DeviceBuffer<int> riftlakemapsurface;
        DeviceBuffer<bool> noshademap;
        DeviceBuffer<std::uint8_t> fullpixels;
        DeviceBuffer<std::uint8_t> displaypixels;
    };

    struct RegionalCudaDeviceBuffers
    {
        DeviceBuffer<short> map;
        DeviceBuffer<int> lakemap;
        DeviceBuffer<short> jantempmap;
        DeviceBuffer<short> jultempmap;
        DeviceBuffer<short> extrajantempmap;
        DeviceBuffer<short> extrajultempmap;
        DeviceBuffer<short> janrainmap;
        DeviceBuffer<short> julrainmap;
        DeviceBuffer<short> climatemap;
        DeviceBuffer<short> seaicemap;
        DeviceBuffer<short> rivermapdir;
        DeviceBuffer<int> rivermapjan;
        DeviceBuffer<int> rivermapjul;
        DeviceBuffer<short> fakeriversdir;
        DeviceBuffer<int> fakeriversjan;
        DeviceBuffer<int> fakeriversjul;
        DeviceBuffer<short> specials;
        DeviceBuffer<short> deltamapdir;
        DeviceBuffer<int> deltamapjan;
        DeviceBuffer<int> deltamapjul;
        DeviceBuffer<bool> volcanomap;
        DeviceBuffer<bool> mudmap;
        DeviceBuffer<bool> sandmap;
        DeviceBuffer<bool> shinglemap;
        DeviceBuffer<int> testmap;
        DeviceBuffer<short> relief1;
        DeviceBuffer<short> relief2;
        DeviceBuffer<short> relief3;
        DeviceBuffer<short> scratch1;
        DeviceBuffer<short> scratch2;
        DeviceBuffer<short> scratch3;
        DeviceBuffer<std::uint8_t> fullpixels;
    };

    bool uploadglobaldevicebuffers(const GlobalCudaRendererInputs& inputs, size_t cellcount, GlobalCudaDeviceBuffers& buffers)
    {
        return uploadbuffer(buffers.mapnom, inputs.mapnom, cellcount)
            && uploadbuffer(buffers.lakemap, inputs.lakemap, cellcount)
            && uploadbuffer(buffers.oceanridgeheightmap, inputs.oceanridgeheightmap, cellcount)
            && uploadbuffer(buffers.mountainheights, inputs.mountainheights, cellcount)
            && uploadbuffer(buffers.volcanomap, inputs.volcanomap, cellcount)
            && uploadbuffer(buffers.extraelevmap, inputs.extraelevmap, cellcount)
            && uploadbuffer(buffers.craterrims, inputs.craterrims, cellcount)
            && uploadbuffer(buffers.cratercentres, inputs.cratercentres, cellcount)
            && uploadbuffer(buffers.oceanridgemap, inputs.oceanridgemap, cellcount)
            && uploadbuffer(buffers.oceanriftmap, inputs.oceanriftmap, cellcount)
            && uploadbuffer(buffers.jantempmap, inputs.jantempmap, cellcount)
            && uploadbuffer(buffers.jultempmap, inputs.jultempmap, cellcount)
            && uploadbuffer(buffers.janrainmap, inputs.janrainmap, cellcount)
            && uploadbuffer(buffers.julrainmap, inputs.julrainmap, cellcount)
            && uploadbuffer(buffers.seaicemap, inputs.seaicemap, cellcount)
            && uploadbuffer(buffers.climatemap, inputs.climatemap, cellcount)
            && uploadbuffer(buffers.rivermapjan, inputs.rivermapjan, cellcount)
            && uploadbuffer(buffers.rivermapjul, inputs.rivermapjul, cellcount)
            && uploadbuffer(buffers.specials, inputs.specials, cellcount)
            && uploadbuffer(buffers.deltamapdir, inputs.deltamapdir, cellcount)
            && uploadbuffer(buffers.deltamapjan, inputs.deltamapjan, cellcount)
            && uploadbuffer(buffers.deltamapjul, inputs.deltamapjul, cellcount)
            && uploadbuffer(buffers.riftlakemapsurface, inputs.riftlakemapsurface, cellcount)
            && uploadbuffer(buffers.noshademap, inputs.noshademap, cellcount);
    }

    GlobalCudaDeviceInputs builddeviceinputs(const GlobalCudaRendererInputs& inputs, const GlobalCudaDeviceBuffers& buffers)
    {
        GlobalCudaDeviceInputs deviceinputs = {};
        deviceinputs.width = inputs.width;
        deviceinputs.height = inputs.height;
        deviceinputs.sealevel = inputs.sealevel;
        deviceinputs.maxelevation = inputs.maxelevation;
        deviceinputs.maxriverflow = inputs.maxriverflow;
        deviceinputs.worldsize = inputs.worldsize;
        deviceinputs.perihelion = inputs.perihelion;
        deviceinputs.eccentricity = inputs.eccentricity;
        deviceinputs.tilt = inputs.tilt;
        deviceinputs.mapnom = buffers.mapnom.data();
        deviceinputs.lakemap = buffers.lakemap.data();
        deviceinputs.oceanridgeheightmap = buffers.oceanridgeheightmap.data();
        deviceinputs.mountainheights = buffers.mountainheights.data();
        deviceinputs.volcanomap = buffers.volcanomap.data();
        deviceinputs.extraelevmap = buffers.extraelevmap.data();
        deviceinputs.craterrims = buffers.craterrims.data();
        deviceinputs.cratercentres = buffers.cratercentres.data();
        deviceinputs.oceanridgemap = buffers.oceanridgemap.data();
        deviceinputs.oceanriftmap = buffers.oceanriftmap.data();
        deviceinputs.jantempmap = buffers.jantempmap.data();
        deviceinputs.jultempmap = buffers.jultempmap.data();
        deviceinputs.janrainmap = buffers.janrainmap.data();
        deviceinputs.julrainmap = buffers.julrainmap.data();
        deviceinputs.seaicemap = buffers.seaicemap.data();
        deviceinputs.climatemap = buffers.climatemap.data();
        deviceinputs.rivermapjan = buffers.rivermapjan.data();
        deviceinputs.rivermapjul = buffers.rivermapjul.data();
        deviceinputs.specials = buffers.specials.data();
        deviceinputs.deltamapdir = buffers.deltamapdir.data();
        deviceinputs.deltamapjan = buffers.deltamapjan.data();
        deviceinputs.deltamapjul = buffers.deltamapjul.data();
        deviceinputs.riftlakemapsurface = buffers.riftlakemapsurface.data();
        deviceinputs.noshademap = buffers.noshademap.data();
        deviceinputs.relief = inputs.relief;
        return deviceinputs;
    }

    bool uploadregionaldevicebuffers(const RegionalCudaRendererInputs& inputs, RegionalCudaDeviceBuffers& buffers)
    {
        constexpr size_t cellcount = static_cast<size_t>(cudaregionalstride) * static_cast<size_t>(cudaregionalstride);

        return uploadbuffer(buffers.map, inputs.map, cellcount)
            && uploadbuffer(buffers.lakemap, inputs.lakemap, cellcount)
            && uploadbuffer(buffers.jantempmap, inputs.jantempmap, cellcount)
            && uploadbuffer(buffers.jultempmap, inputs.jultempmap, cellcount)
            && uploadbuffer(buffers.extrajantempmap, inputs.extrajantempmap, cellcount)
            && uploadbuffer(buffers.extrajultempmap, inputs.extrajultempmap, cellcount)
            && uploadbuffer(buffers.janrainmap, inputs.janrainmap, cellcount)
            && uploadbuffer(buffers.julrainmap, inputs.julrainmap, cellcount)
            && uploadbuffer(buffers.climatemap, inputs.climatemap, cellcount)
            && uploadbuffer(buffers.seaicemap, inputs.seaicemap, cellcount)
            && uploadbuffer(buffers.rivermapdir, inputs.rivermapdir, cellcount)
            && uploadbuffer(buffers.rivermapjan, inputs.rivermapjan, cellcount)
            && uploadbuffer(buffers.rivermapjul, inputs.rivermapjul, cellcount)
            && uploadbuffer(buffers.fakeriversdir, inputs.fakeriversdir, cellcount)
            && uploadbuffer(buffers.fakeriversjan, inputs.fakeriversjan, cellcount)
            && uploadbuffer(buffers.fakeriversjul, inputs.fakeriversjul, cellcount)
            && uploadbuffer(buffers.specials, inputs.specials, cellcount)
            && uploadbuffer(buffers.deltamapdir, inputs.deltamapdir, cellcount)
            && uploadbuffer(buffers.deltamapjan, inputs.deltamapjan, cellcount)
            && uploadbuffer(buffers.deltamapjul, inputs.deltamapjul, cellcount)
            && uploadbuffer(buffers.volcanomap, inputs.volcanomap, cellcount)
            && uploadbuffer(buffers.mudmap, inputs.mudmap, cellcount)
            && uploadbuffer(buffers.sandmap, inputs.sandmap, cellcount)
            && uploadbuffer(buffers.shinglemap, inputs.shinglemap, cellcount)
            && uploadbuffer(buffers.testmap, inputs.testmap, cellcount);
    }

    RegionalCudaDeviceInputs buildregionaldeviceinputs(const RegionalCudaRendererInputs& inputs, const RegionalCudaDeviceBuffers& buffers)
    {
        RegionalCudaDeviceInputs deviceinputs = {};
        deviceinputs.sealevel = inputs.sealevel;
        deviceinputs.maxelevation = inputs.maxelevation;
        deviceinputs.maxriverflow = inputs.maxriverflow;
        deviceinputs.regwidthbegin = inputs.regwidthbegin;
        deviceinputs.regwidthend = inputs.regwidthend;
        deviceinputs.regheightbegin = inputs.regheightbegin;
        deviceinputs.regheightend = inputs.regheightend;
        deviceinputs.regionwidth = inputs.regionwidth;
        deviceinputs.regionheight = inputs.regionheight;
        deviceinputs.lefty = inputs.lefty;
        deviceinputs.globalheight = inputs.globalheight;
        deviceinputs.map = buffers.map.data();
        deviceinputs.lakemap = buffers.lakemap.data();
        deviceinputs.jantempmap = buffers.jantempmap.data();
        deviceinputs.jultempmap = buffers.jultempmap.data();
        deviceinputs.extrajantempmap = buffers.extrajantempmap.data();
        deviceinputs.extrajultempmap = buffers.extrajultempmap.data();
        deviceinputs.janrainmap = buffers.janrainmap.data();
        deviceinputs.julrainmap = buffers.julrainmap.data();
        deviceinputs.climatemap = buffers.climatemap.data();
        deviceinputs.seaicemap = buffers.seaicemap.data();
        deviceinputs.rivermapdir = buffers.rivermapdir.data();
        deviceinputs.rivermapjan = buffers.rivermapjan.data();
        deviceinputs.rivermapjul = buffers.rivermapjul.data();
        deviceinputs.fakeriversdir = buffers.fakeriversdir.data();
        deviceinputs.fakeriversjan = buffers.fakeriversjan.data();
        deviceinputs.fakeriversjul = buffers.fakeriversjul.data();
        deviceinputs.specials = buffers.specials.data();
        deviceinputs.deltamapdir = buffers.deltamapdir.data();
        deviceinputs.deltamapjan = buffers.deltamapjan.data();
        deviceinputs.deltamapjul = buffers.deltamapjul.data();
        deviceinputs.volcanomap = buffers.volcanomap.data();
        deviceinputs.mudmap = buffers.mudmap.data();
        deviceinputs.sandmap = buffers.sandmap.data();
        deviceinputs.shinglemap = buffers.shinglemap.data();
        deviceinputs.testmap = buffers.testmap.data();
        deviceinputs.relief = inputs.relief;
        return deviceinputs;
    }

    template <typename Kernel>
    bool renderregionalmap(Kernel kernel, const RegionalCudaDeviceInputs& inputs, DeviceBuffer<std::uint8_t>& fullpixels, std::vector<std::uint8_t>& output)
    {
        const int imagewidth = inputs.regwidthend - inputs.regwidthbegin + 1;
        const int imageheight = inputs.regheightend - inputs.regheightbegin + 1;
        const size_t fullbytes = static_cast<size_t>(imagewidth) * static_cast<size_t>(imageheight) * rgbachannels;

        if (!fullpixels.ensure(fullbytes))
            return false;

        if (!cudacheck(cudaMemset(fullpixels.data(), 0, fullbytes)))
            return false;

        const dim3 block(16, 16);
        const dim3 grid((imagewidth + 15) / 16, (imageheight + 15) / 16);
        kernel<<<grid, block>>>(inputs, fullpixels.data(), imagewidth);

        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        return copyoutput(output, fullpixels, fullbytes);
    }

    bool renderregionalrelief(const RegionalCudaDeviceInputs& inputs, RegionalCudaDeviceBuffers& buffers, std::vector<std::uint8_t>& output)
    {
        const int imagewidth = inputs.regwidthend - inputs.regwidthbegin + 1;
        const int imageheight = inputs.regheightend - inputs.regheightbegin + 1;
        const size_t pixelcount = static_cast<size_t>(cudaregionalstride) * static_cast<size_t>(cudaregionalstride);
        const size_t fullbytes = static_cast<size_t>(imagewidth) * static_cast<size_t>(imageheight) * rgbachannels;

        if (!buffers.relief1.ensure(pixelcount) || !buffers.relief2.ensure(pixelcount) || !buffers.relief3.ensure(pixelcount)
            || !buffers.scratch1.ensure(pixelcount) || !buffers.scratch2.ensure(pixelcount) || !buffers.scratch3.ensure(pixelcount)
            || !buffers.fullpixels.ensure(fullbytes))
        {
            return false;
        }

        const dim3 block(16, 16);
        const dim3 grid((imagewidth + 15) / 16, (imageheight + 15) / 16);

        renderregionalreliefbasekernel<<<grid, block>>>(inputs, buffers.relief1.data(), buffers.relief2.data(), buffers.relief3.data());
        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        blurregionalwetlandskernel<<<grid, block>>>(inputs, buffers.relief1.data(), buffers.relief2.data(), buffers.relief3.data(), buffers.scratch1.data(), buffers.scratch2.data(), buffers.scratch3.data());
        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        blurregionalvalleyskernel<<<grid, block>>>(inputs, buffers.scratch1.data(), buffers.scratch2.data(), buffers.scratch3.data(), buffers.relief1.data(), buffers.relief2.data(), buffers.relief3.data());
        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        redrawregionalreliefriverskernel<<<grid, block>>>(inputs, buffers.relief1.data(), buffers.relief2.data(), buffers.relief3.data());
        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        if (!cudacheck(cudaMemset(buffers.fullpixels.data(), 0, fullbytes)))
            return false;

        shaderegionalreliefkernel<<<grid, block>>>(inputs, buffers.relief1.data(), buffers.relief2.data(), buffers.relief3.data(), buffers.fullpixels.data(), imagewidth);
        if (!cudacheck(cudaGetLastError()) || !cudacheck(cudaDeviceSynchronize()))
            return false;

        return copyoutput(output, buffers.fullpixels, fullbytes);
    }
}

bool cudaglobalrenderersavailable()
{
    static int availability = -1;

    if (availability != -1)
        return availability == 1;

    int devicecount = 0;

    if (!cudacheck(cudaGetDeviceCount(&devicecount)) || devicecount == 0)
    {
        availability = 0;
        return false;
    }

    availability = 1;
    return true;
}

bool renderglobalnonreliefmapscuda(const GlobalCudaRendererInputs& inputs, GlobalCudaRendererOutputs& outputs)
{
    if (!cudaglobalrenderersavailable())
        return false;

    const size_t cellcount = static_cast<size_t>(inputs.width + 1) * cudaplanetstride;

    GlobalCudaDeviceBuffers buffers;

    if (!uploadglobaldevicebuffers(inputs, cellcount, buffers))
        return false;

    const GlobalCudaDeviceInputs deviceinputs = builddeviceinputs(inputs, buffers);

    if (!renderandcopy(renderelevationkernel, nearestdisplaykernel, deviceinputs, buffers.fullpixels, buffers.displaypixels, outputs.elevation, outputs.elevationDisplay))
        return false;

    if (!renderandcopy(rendertemperaturekernel, outlineddisplaykernel, deviceinputs, buffers.fullpixels, buffers.displaypixels, outputs.temperature, outputs.temperatureDisplay))
        return false;

    if (!renderandcopy(renderprecipitationkernel, outlineddisplaykernel, deviceinputs, buffers.fullpixels, buffers.displaypixels, outputs.precipitation, outputs.precipitationDisplay))
        return false;

    if (!renderandcopy(renderclimatekernel, outlineddisplaykernel, deviceinputs, buffers.fullpixels, buffers.displaypixels, outputs.climate, outputs.climateDisplay))
        return false;

    if (!renderrivers(deviceinputs, buffers.fullpixels, buffers.displaypixels, outputs.rivers, outputs.riversDisplay))
        return false;

    return true;
}

bool renderglobalreliefmapcuda(const GlobalCudaRendererInputs& inputs, std::vector<std::uint8_t>& relief, std::vector<std::uint8_t>& reliefDisplay)
{
    if (!cudaglobalrenderersavailable())
        return false;

    const size_t cellcount = static_cast<size_t>(inputs.width + 1) * cudaplanetstride;
    GlobalCudaDeviceBuffers buffers;

    if (!uploadglobaldevicebuffers(inputs, cellcount, buffers))
        return false;

    const GlobalCudaDeviceInputs deviceinputs = builddeviceinputs(inputs, buffers);
    return renderandcopy(renderglobalreliefkernel, nearestdisplaykernel, deviceinputs, buffers.fullpixels, buffers.displaypixels, relief, reliefDisplay);
}

bool renderregionalnonreliefmapscuda(const RegionalCudaRendererInputs& inputs, RegionalCudaRendererOutputs& outputs)
{
    if (!cudaglobalrenderersavailable())
        return false;

    RegionalCudaDeviceBuffers buffers;

    if (!uploadregionaldevicebuffers(inputs, buffers))
        return false;

    const RegionalCudaDeviceInputs deviceinputs = buildregionaldeviceinputs(inputs, buffers);

    if (!renderregionalmap(renderregionalelevationkernel, deviceinputs, buffers.fullpixels, outputs.elevation))
        return false;

    if (!renderregionalmap(renderregionaltemperaturekernel, deviceinputs, buffers.fullpixels, outputs.temperature))
        return false;

    if (!renderregionalmap(renderregionalprecipitationkernel, deviceinputs, buffers.fullpixels, outputs.precipitation))
        return false;

    if (!renderregionalmap(renderregionalclimatekernel, deviceinputs, buffers.fullpixels, outputs.climate))
        return false;

    if (!renderregionalmap(renderregionalriverskernel, deviceinputs, buffers.fullpixels, outputs.rivers))
        return false;

    return true;
}

bool renderregionalreliefmapcuda(const RegionalCudaRendererInputs& inputs, std::vector<std::uint8_t>& relief)
{
    if (!cudaglobalrenderersavailable())
        return false;

    RegionalCudaDeviceBuffers buffers;

    if (!uploadregionaldevicebuffers(inputs, buffers))
        return false;

    const RegionalCudaDeviceInputs deviceinputs = buildregionaldeviceinputs(inputs, buffers);
    return renderregionalrelief(deviceinputs, buffers, relief);
}
