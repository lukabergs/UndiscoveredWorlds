#include "cuda_renderers.hpp"

#ifndef UW_ENABLE_CUDA_RENDERERS

bool cudaglobalrenderersavailable()
{
    return false;
}

bool renderglobalnonreliefmapscuda(const GlobalCudaRendererInputs&, GlobalCudaRendererOutputs&)
{
    return false;
}

bool renderglobalreliefmapcuda(const GlobalCudaRendererInputs&, std::vector<std::uint8_t>&, std::vector<std::uint8_t>&)
{
    return false;
}

bool renderregionalnonreliefmapscuda(const RegionalCudaRendererInputs&, RegionalCudaRendererOutputs&)
{
    return false;
}

bool renderregionalreliefmapcuda(const RegionalCudaRendererInputs&, std::vector<std::uint8_t>&)
{
    return false;
}

#endif
