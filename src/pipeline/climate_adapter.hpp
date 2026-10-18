#pragma once

// The climate adapter (docs/architecture/sphere-mesh.md section 10). Climate keeps
// its lat-lon process grids; on a cube world it runs on a persistent lat-lon
// shadow world of the export width (4 N), kept across climate processes in the
// generation scratch. Before each climate process the fields it reads (its
// declared access) that changed on the cube since the last sync are exported to
// the shadow by the catalogue's remap rules; after it the fields it writes are
// imported to the cube (importfields). Climate-owned state stays authoritative on
// the shadow: a field the cube received from the shadow is not sent back unless
// another process changes it on the cube. On a lat-lon world the adapter is the
// identity: the process runs on the world itself and no shadow exists.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

#include "pipeline/generation_stage_registry.hpp"

class planet;
class FieldRemap;

namespace mesh
{
class CellGrid;
}

// What the last climate process moved between the cube and its shadow.
struct ClimateShadowSync
{
    bool created = false;            // the shadow was (re)built by a full export before this process
    std::size_t exported = 0;        // changed read fields sent to the shadow (0 when created)
    std::size_t imported = 0;        // written fields received from the shadow
    std::size_t importcalls = 0;     // importfields calls (one per process with fields to import)
    double exportseconds = 0.0;      // change detection, export and recount
    double importseconds = 0.0;
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
    // Drops the shadow: the next climate process on a cube world rebuilds it.
    void reset();

private:
    friend void runclimateprocess(GenerationStageId, planet&, ClimateShadow&, const std::function<void(planet&)>&);

    std::unique_ptr<planet> shadow_;
    const mesh::CellGrid* source_ = nullptr;     // the cube grid the shadow mirrors
    std::vector<std::uint64_t> hashes_;          // per field id: the cube's value at the last sync
    std::vector<std::uint8_t> tracked_;          // per field id: read by some climate process
    mutable std::unique_ptr<FieldRemap> remap_;  // cube -> shadow, weights built on first use
    ClimateShadowSync last_;
};

// Runs one climate process on the world. On a lat-lon world `body` runs on the
// world itself. On a cube world it runs on the shadow, between the export of the
// process's changed reads and the import of its writes; the shadow's climate
// settings (resolution, coupling iterations, solver options) follow the world's.
void runclimateprocess(GenerationStageId stage, planet& world, ClimateShadow& shadow, const std::function<void(planet&)>& body);
