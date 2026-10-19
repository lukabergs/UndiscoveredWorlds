#pragma once

// The climate adapter (docs/architecture/sphere-mesh.md section 10). Climate keeps
// its lat-lon process grids; on a cube world it runs on a persistent lat-lon
// shadow world of the export width (4 N), kept across climate processes in the
// generation scratch and reset whenever the world is cleared.
// - Before each climate process every field a climate process reads or writes
//   that changed on the cube since the last sync is exported to the shadow by the
//   catalogue's remap rules.
// - After it the declared writes the process changed on the shadow are imported
//   to the cube (importfields, once per process); a cube cell that no changed
//   shadow cell reaches keeps its value, so a process that rewrites part of a
//   field does not send the rest round the remap. A field the process changed
//   but did not declare is an error naming the process and the fields.
// - Climate-owned state stays authoritative on the shadow: the solver's own grids,
//   the atmosphere levels, the solve diagnostics, the zonal profiles and the
//   direction-coded fields (the orographic rain directions, D8 codes in the
//   lat-lon frame) are never exported or imported; the cube has none of them, and
//   whatever needs them on a cube world reads the shadow. Cell counts are
//   recounted on the shadow's own terrain.
// On a lat-lon world the adapter is the identity: the process runs on the world
// itself and no shadow exists.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

#include "core/world/world_fields.hpp"
#include "pipeline/generation_stage_registry.hpp"

class planet;

namespace mesh
{
struct RemapWeights;
}

// What the last climate process moved between the cube and its shadow.
struct ClimateShadowSync
{
    bool created = false;            // the shadow was built by a full export before this process
    std::size_t exported = 0;        // changed fields sent to the shadow (0 when created)
    std::size_t imported = 0;        // written fields received from the shadow
    std::size_t importcalls = 0;     // importfields calls: 1 when a field was imported, else 0
    std::size_t keptcells = 0;       // cube cell values of imported fields left as they were
    double exportseconds = 0.0;      // change detection, export, recount and the shadow snapshot
    double importseconds = 0.0;      // write check and import
};

class ClimateShadow
{
public:
    ClimateShadow();
    ~ClimateShadow();
    ClimateShadow(const ClimateShadow& other);
    ClimateShadow& operator=(const ClimateShadow& other);
    ClimateShadow(ClimateShadow&&) noexcept;
    ClimateShadow& operator=(ClimateShadow&&) noexcept;

    // The lat-lon shadow world; null before the first climate process on a cube world.
    const planet* world() const { return shadow_.get(); }
    const ClimateShadowSync& lastsync() const { return last_; }
    // Drops the shadow: the next climate process on a cube world rebuilds it. Called
    // where the world is cleared (a new generation).
    void reset();

private:
    friend void runclimateprocess(GenerationStageId, planet&, ClimateShadow&, const std::function<void(planet&)>&);

    std::unique_ptr<planet> shadow_;
    int facecells_ = 0;                                    // the cube's N when the shadow was built
    std::vector<std::uint64_t> hashes_;                    // per field id: the cube's value at the last sync
    std::vector<std::uint8_t> tracked_;                    // per field id: read or written by a climate process
    std::shared_ptr<const mesh::RemapWeights> importweights_;   // shadow -> cube, for the reach of changed cells
    ClimateShadowSync last_;
};

// Climate-owned state that stays authoritative on the shadow (above; of the direction
// fields only the atmosphere's): what a cube world's exports and bundles take from
// the shadow (export_routing.hpp).
bool climateshadowowned(worldmodel::FieldId id);

// Runs one climate process on the world. On a lat-lon world `body` runs on the
// world itself. On a cube world it runs on the shadow, between the export of the
// changed fields and the import of its writes; the shadow's climate settings
// (resolution, coupling iterations, solver options) follow the world's.
void runclimateprocess(GenerationStageId stage, planet& world, ClimateShadow& shadow, const std::function<void(planet&)>& body);
