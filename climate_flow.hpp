#pragma once

#include <vector>

namespace climateflow
{
struct ParticlePoint
{
    float x = 0.0f;
    float y = 0.0f;
    float elapsedSeconds = 0.0f;
};

struct ParticleTraceConfig
{
    float planetRadiusMetres = 6371000.0f;
    float maximumSubstepCells = 0.65f;
    float maximumDirectionChangeRadians = 0.20f;
    int maximumSubsteps = 128;
    bool poleInclusiveRows = true;
};

struct ParticleStep
{
    std::vector<ParticlePoint> points;
    bool remainedInDomain = true;
};

ParticleStep advectParticleRk2(
    int columns,
    int rows,
    const std::vector<float>& eastWindMps,
    const std::vector<float>& southWindMps,
    float startX,
    float startY,
    float elapsedSeconds,
    const ParticleTraceConfig& config = {});
}
