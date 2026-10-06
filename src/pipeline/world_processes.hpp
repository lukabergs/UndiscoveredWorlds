#pragma once

// Process registry: each generation stage is one process with a kind and a
// declared access set over the world catalogue (docs/architecture/world-model.md).
// Access sets were measured with the read-tracking analysis build and the
// scheduler's write audit; keep them current when a stage changes.

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "pipeline/generation_stage_registry.hpp"
#include "core/world/world_fields.hpp"

enum class ProcessKind : std::uint8_t
{
    initializer,   // creates fields from a seed and parameters
    evolver,       // advances state through time or solves for an equilibrium
    diagnoser,     // derives fields from existing state
    detailer,      // adds plausible small-scale detail from noise
    editor         // user-directed, non-physical adjustment
};

// The time of an evolver of the time model (docs/architecture/world-model.md):
// its duration and step in Myr. A module parameter holds them where one is named
// (the compiled default, set by --set, a recipe step's override or
// --duration <evolver>=<Myr>); `source` says where they come from otherwise.
// A process without a source has no time of its own.
struct EvolverTime
{
    const char* durationParameter = nullptr;   // <module>.<group>.<key>
    const char* stepParameter = nullptr;
    const char* source = nullptr;
};

struct WorldProcessDefinition
{
    GenerationStageId stage;
    const char* id;                  // module-qualified, e.g. "climate.winds"
    ProcessKind kind;
    bool stochastic;                 // draws from a seeded random stream
    std::vector<worldmodel::FieldId> reads;
    std::vector<worldmodel::FieldId> writes;
    std::vector<worldmodel::TableId> tablereads;
    std::vector<worldmodel::TableId> tablewrites;
    const char* note;                // known responsibilities outside the label
    EvolverTime time;                // evolvers of the time model only
};

const char* processkindname(ProcessKind kind);
const std::vector<WorldProcessDefinition>& getworldprocesses();
const WorldProcessDefinition& getworldprocess(GenerationStageId stage);
const WorldProcessDefinition* findworldprocess(std::string_view id);
// The ids of the evolvers of the time model, comma-separated, for messages and help.
std::string timedevolverids();
// Seeds the shared random generator for one step from the world seed and the
// process id: each process has its own stream, so skipping or repeating a step
// leaves the streams of the others unchanged (docs/architecture/world-model.md).
void seedprocessstream(long worldseed, GenerationStageId stage);
bool declareswrite(const WorldProcessDefinition& process, worldmodel::FieldId field);
bool declareswrite(const WorldProcessDefinition& process, worldmodel::TableId table);
