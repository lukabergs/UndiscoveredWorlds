#pragma once

// Run bundles: one NumPy .npy file per field, tables as text and a
// manifest.json (schema uw-run-bundle/v1). The writer and the reader are
// independent of the world class; see docs/architecture/world-model.md.

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace runbundle
{
inline constexpr const char* kSchema = "uw-run-bundle/v1";

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
    double seconds = 0.0;
    std::vector<std::string> writes;     // fields whose contents changed (audit)
};

struct Metadata
{
    long seed = 0;
    int columns = 0;
    int rows = 0;
    std::string recipe;
    std::vector<std::string> command;
    std::vector<StepRecord> steps;
    std::vector<std::pair<std::string, std::string>> properties;   // extra string facts
};

class Writer
{
public:
    // Refuses an existing non-empty directory so earlier evidence is kept.
    explicit Writer(std::filesystem::path directory);
    void addfield(const FieldEntry& field);
    void addtable(const TableEntry& table);
    void finish(const Metadata& metadata);
    const std::filesystem::path& directory() const { return directory_; }

private:
    struct WrittenField { FieldEntry entry; std::string file; std::uint64_t hash; };
    struct WrittenTable { TableEntry entry; std::string file; std::uint64_t hash; };
    std::filesystem::path directory_;
    std::vector<WrittenField> fields_;
    std::vector<WrittenTable> tables_;
};

// A field or table as the manifest lists it.
struct ManifestEntry
{
    std::string id;
    std::string component;
    std::string layout;       // fields only
    std::string dtype;        // fields only
    std::vector<std::size_t> shape;
    std::string units;
    std::string description;
    std::string file;         // relative to the bundle directory
    std::string writtenby;
};

// Reads a bundle's manifest; field data and table text are read on request.
class Reader
{
public:
    // Throws when the directory holds no manifest of this schema.
    explicit Reader(std::filesystem::path directory);
    const std::filesystem::path& directory() const { return directory_; }
    long seed() const { return seed_; }
    int columns() const { return columns_; }
    int rows() const { return rows_; }
    const std::string& recipe() const { return recipe_; }
    const std::vector<ManifestEntry>& fields() const { return fields_; }
    const std::vector<ManifestEntry>& tables() const { return tables_; }
    // The field's elements, checked against the manifest's dtype and shape.
    std::vector<unsigned char> fielddata(const ManifestEntry& field) const;
    std::string tabletext(const ManifestEntry& table) const;

private:
    std::filesystem::path directory_;
    long seed_ = 0;
    int columns_ = 0;
    int rows_ = 0;
    std::string recipe_;
    std::vector<ManifestEntry> fields_;
    std::vector<ManifestEntry> tables_;
};

void writenpy(const std::filesystem::path& path, const std::string& dtype,
    const std::vector<std::size_t>& shape, const void* data, std::size_t bytes);
std::string jsonstring(std::string_view text);
std::uint64_t hashbytes(const void* data, std::size_t size);
}
