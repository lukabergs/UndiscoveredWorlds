#include "planet_bundle.hpp"

#include "planet.hpp"
#include "planet_fields.hpp"

void writeplanetbundle(const planet& world, const std::filesystem::path& directory,
    const runbundle::Metadata& metadata, const BundleProvenance& provenance)
{
    using namespace worldmodel;
    runbundle::Writer writer(directory);

    for (std::size_t index = 0; index < kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        const FieldSpec& spec = fieldspec(id);
        const PlanetFieldData data = readplanetfield(world, id);
        if (data.shape.empty() && data.bytes.empty())
            continue;   // absent (a diagnostic the world has not recorded)
        runbundle::FieldEntry entry;
        entry.id = qualifiedname(id);
        entry.component = componentname(spec.component);
        entry.units = spec.units;
        entry.description = spec.description;
        entry.layout = layoutname(spec.layout);
        entry.dtype = elementname(spec.element);
        entry.shape = data.shape;
        entry.data = data.bytes.data();
        entry.bytes = data.bytes.size();
        if (const auto found = provenance.fields.find(id); found != provenance.fields.end())
            entry.writtenby = found->second;
        writer.addfield(entry);
    }

    for (std::size_t index = 0; index < kTableCount; index++)
    {
        const TableId id = static_cast<TableId>(index);
        const TableSpec& spec = tablespec(id);
        runbundle::TableEntry entry;
        entry.id = qualifiedname(id);
        entry.component = componentname(spec.component);
        entry.description = spec.description;
        entry.text = planettabletext(world, id, &entry.extension);
        if (const auto found = provenance.tables.find(id); found != provenance.tables.end())
            entry.writtenby = found->second;
        writer.addtable(entry);
    }

    writer.finish(metadata);
}
