#pragma once

// Run bundles: one NumPy .npy file per field, tables as text and a
// manifest.json (schema uw-run-bundle/v2; v1 is read as grid kind latlon). The
// writer and the reader are independent of the world class; see
// docs/architecture/world-model.md and docs/architecture/sphere-mesh.md section 7.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace runbundle
{
inline constexpr const char* kSchema = "uw-run-bundle/v2";
inline constexpr const char* kSchemaV1 = "uw-run-bundle/v1";   // read, never written

struct FieldEntry
{
    std::string id;           // "component.name"
    std::string component;
    std::string units;
    std::string description;
    std::string layout;
    std::string dtype;        // NumPy descriptor, e.g. "<f4"
    std::vector<std::size_t> shape;
    const void* data = nullptr;
    std::size_t bytes = 0;
    std::string writtenby;    // step label from provenance; empty when unknown
};

// A plain array of the cell table (grid/<id>.npy).
struct ArrayEntry
{
    std::string id;           // "centre", "latitude_deg", ...
    std::string dtype;
    std::vector<std::size_t> shape;
    const void* data = nullptr;
    std::size_t bytes = 0;
};

struct TableEntry
{
    std::string id;
    std::string component;
    std::string description;
    std::string extension;    // ".csv", ".json", ".hex"
    std::string text;
    std::string writtenby;
};

struct StepRecord
{
    int index = 0;
    std::string process;
    std::string label;
    // Settings the recipe changed for this step only (name, JSON value); written when not empty.
    std::vector<std::pair<std::string, std::string>> overrides;
    double seconds = 0.0;
    std::vector<std::string> writes;     // fields whose contents changed (audit)
};

// The manifest's grid block beyond Metadata::columns and rows (sphere-mesh.md
// section 7). latlon: the columns x rows grid. cube: facecells N per face edge,
// orientation uw-cube-1, fields [layer,] 6, N, N; columns x rows is its lat-lon
// export (export_columns W = 4 N). The cube-window kind comes with the windows (W2-W).
struct GridSpec
{
    std::string kind = "latlon";
    int facecells = 0;
    double radiuskm = 0.0;    // written when positive
};

struct Metadata
{
    long seed = 0;
    int columns = 0;          // the lat-lon grid, or the cube's lat-lon export
    int rows = 0;
    GridSpec grid;
    std::string recipe;
    std::vector<std::string> command;
    std::vector<StepRecord> steps;
    std::vector<std::pair<std::string, std::string>> properties;   // extra string facts
    std::string parameters;   // a JSON object of the modules' tunable values, or empty
    // Extra top-level members (name, JSON value) written after "parameters" in this order;
    // a regional bundle's window, source and diagnostics (io/regional_bundle.hpp). World
    // bundles have none.
    std::vector<std::pair<std::string, std::string>> blocks;
};

class Writer
{
public:
    // Refuses an existing non-empty directory so earlier evidence is kept.
    explicit Writer(std::filesystem::path directory);
    void addfield(const FieldEntry& field);
    // A field's lat-lon export (latlon/<component>/<name>.npy) by its remap rule; a
    // direction field records whether it was recomputed on the export.
    void addexport(const FieldEntry& field, const std::string& rule, bool recomputed);
    // One array of the cell table (grid/<id>.npy); cube kinds only.
    void addcells(const ArrayEntry& cells);
    void addtable(const TableEntry& table);
    void finish(const Metadata& metadata);
    const std::filesystem::path& directory() const { return directory_; }

private:
    struct WrittenField { FieldEntry entry; std::string file; std::uint64_t hash; };
    struct WrittenExport { FieldEntry entry; std::string file; std::uint64_t hash; std::string rule; bool recomputed; };
    struct WrittenArray { ArrayEntry entry; std::string file; std::uint64_t hash; };
    struct WrittenTable { TableEntry entry; std::string file; std::uint64_t hash; };
    std::filesystem::path directory_;
    std::vector<WrittenField> fields_;
    std::vector<WrittenExport> exports_;
    std::vector<WrittenArray> cells_;
    std::vector<WrittenTable> tables_;
};

// A field, export, cell array or table as the manifest lists it.
struct ManifestEntry
{
    std::string id;
    std::string component;
    std::string layout;       // fields only
    std::string dtype;        // fields, exports and cell arrays
    std::vector<std::size_t> shape;
    std::string units;
    std::string description;
    std::string file;         // relative to the bundle directory
    std::string writtenby;
    std::string rule;         // exports only
    bool recomputed = false;  // exports of direction fields
};

// Reads a bundle's manifest; field data and table text are read on request.
class Reader
{
public:
    // Throws when the directory holds no manifest of schema v1 or v2.
    explicit Reader(std::filesystem::path directory);
    const std::filesystem::path& directory() const { return directory_; }
    const std::string& schema() const { return schema_; }
    long seed() const { return seed_; }
    // The lat-lon grid, or the cube's lat-lon export.
    int columns() const { return columns_; }
    int rows() const { return rows_; }
    const GridSpec& grid() const { return grid_; }   // v1: kind latlon, no radius
    const std::string& recipe() const { return recipe_; }
    const std::vector<ManifestEntry>& fields() const { return fields_; }
    const std::vector<ManifestEntry>& exports() const { return exports_; }
    const std::vector<ManifestEntry>& cells() const { return cells_; }
    const std::vector<ManifestEntry>& tables() const { return tables_; }
    // The array's elements (a field, export or cell array), checked against the manifest's dtype and shape.
    std::vector<unsigned char> fielddata(const ManifestEntry& field) const;
    std::string tabletext(const ManifestEntry& table) const;

private:
    std::filesystem::path directory_;
    std::string schema_;
    long seed_ = 0;
    int columns_ = 0;
    int rows_ = 0;
    GridSpec grid_;
    std::string recipe_;
    std::vector<ManifestEntry> fields_;
    std::vector<ManifestEntry> exports_;
    std::vector<ManifestEntry> cells_;
    std::vector<ManifestEntry> tables_;
};

void writenpy(const std::filesystem::path& path, const std::string& dtype,
    const std::vector<std::size_t>& shape, const void* data, std::size_t bytes);
std::string jsonstring(std::string_view text);
std::uint64_t hashbytes(const void* data, std::size_t size);
}
