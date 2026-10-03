#pragma once

#include <filesystem>
#include <map>
#include <string>

#include "run_bundle.hpp"
#include "world/world_fields.hpp"

class planet;

// Last writer of each field or table, by step label; supplied by the scheduler.
struct BundleProvenance
{
    std::map<worldmodel::FieldId, std::string> fields;
    std::map<worldmodel::TableId, std::string> tables;
};

// Writes every catalogue field and entity table of the world as a run bundle.
void writeplanetbundle(const planet& world, const std::filesystem::path& directory,
    const runbundle::Metadata& metadata, const BundleProvenance& provenance = {});
