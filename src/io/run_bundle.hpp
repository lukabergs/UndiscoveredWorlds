#pragma once

// Run bundle writer: one NumPy .npy file per field, tables as text and a
// manifest.json (schema uw-run-bundle/v1). Independent of the world class;
// see docs/architecture/world-model.md.

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

void writenpy(const std::filesystem::path& path, const std::string& dtype,
    const std::vector<std::size_t>& shape, const void* data, std::size_t bytes);
std::string jsonstring(std::string_view text);
std::uint64_t hashbytes(const void* data, std::size_t size);
}
