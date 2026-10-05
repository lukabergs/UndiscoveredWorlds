#ifndef TECTONIC_SNAPSHOTS_HPP
#define TECTONIC_SNAPSHOTS_HPP

#include "tectonic_contract.hpp"

#include <cstdint>
#include <filesystem>
#include <string>

namespace tectonics::tooling {

struct SnapshotBundleOptions {
    std::string schema = "tectonic-snapshot-bundle/v1";
    std::string label = "snapshot";
    std::string source = "runtime";
    bool include_update_count = false;
    uint32_t update_count = 0;
};

struct LoadedSnapshotBundle {
    tectonics::contract::Snapshot snapshot;
    std::string schema;
    std::string label;
    std::string source;
    uint32_t declared_boundary_segment_count = 0;
    uint32_t declared_deforming_region_count = 0;
    uint32_t declared_tectonic_junction_count = 0;
};

void write_snapshot_bundle(const std::filesystem::path& output_dir,
                           const tectonics::contract::Snapshot& snapshot,
                           const SnapshotBundleOptions& options);
LoadedSnapshotBundle load_snapshot_bundle(const std::filesystem::path& input_dir);
// The record of a snapshot as JSON: label, time, scenario, engine diagnostics,
// statistics and plates. It is the manifest of a snapshot bundle without the
// file and hash tables, and the geosphere.tectonic_snapshot table of a run bundle.
std::string snapshot_record_json(const tectonics::contract::Snapshot& snapshot,
                                 const SnapshotBundleOptions& options);

} // namespace tectonics::tooling

#endif
