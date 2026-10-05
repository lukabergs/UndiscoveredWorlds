#include "world/planet_bundle.hpp"

#include <algorithm>
#include <stdexcept>

#include "world/planet.hpp"
#include "world/planet_fields.hpp"

namespace
{
template<typename Id>
bool contains(const std::vector<Id>& ids, Id id)
{
    return std::find(ids.begin(), ids.end(), id) != ids.end();
}

const runbundle::ManifestEntry* findentry(const std::vector<runbundle::ManifestEntry>& entries, const std::string& id)
{
    for (const runbundle::ManifestEntry& entry : entries)
        if (entry.id == id)
            return &entry;
    return nullptr;
}
}

PlanetBundleSource loadplanetbundle(planet& world, const std::filesystem::path& directory)
{
    using namespace worldmodel;
    PlanetBundleSource source;
    source.reader = std::make_shared<runbundle::Reader>(directory);
    const runbundle::Reader& reader = *source.reader;
    if (reader.columns() < 2 || reader.columns() != reader.rows() * 2)
        throw std::runtime_error("Run bundle grid is not W x W/2: " + directory.string());
    world.setwidth(reader.columns() - 1);
    if (world.width() + 1 != reader.columns())
        throw std::runtime_error("Run bundle grid exceeds the world's capacity: " + directory.string());

    // Tables are compared before fields load: the crater table follows a scalar field.
    for (const runbundle::ManifestEntry& entry : reader.tables())
    {
        TableId id{};
        if (!findtable(entry.id, id))
            throw std::runtime_error("Run bundle has a table outside the catalogue: " + entry.id);
        if (reader.tabletext(entry) == planettabletext(world, id))
            source.loadedtables.push_back(id);
        else
            source.carriedtables.push_back(id);
        if (!entry.writtenby.empty())
            source.provenance.tables[id] = entry.writtenby;
    }

    for (const runbundle::ManifestEntry& entry : reader.fields())
    {
        PlanetFieldData data;
        if (!findfield(entry.id, data.id))
            throw std::runtime_error("Run bundle has a field outside the catalogue: " + entry.id);
        if (entry.dtype != elementname(fieldspec(data.id).element))
            throw std::runtime_error("Run bundle field has another element type than the catalogue: " + entry.id);
        data.shape = entry.shape;
        data.bytes = reader.fielddata(entry);
        if (loadplanetfield(world, data))
            source.loadedfields.push_back(data.id);
        else
            source.carriedfields.push_back(data.id);
        if (!entry.writtenby.empty())
            source.provenance.fields[data.id] = entry.writtenby;
    }
    return source;
}

void writeplanetbundle(const planet& world, const std::filesystem::path& directory,
    const runbundle::Metadata& metadata, const BundleProvenance& provenance,
    const PlanetBundleSource* source, const std::vector<worldmodel::FieldId>& writtenfields,
    const std::vector<worldmodel::TableId>& writtentables)
{
    using namespace worldmodel;
    runbundle::Writer writer(directory);

    for (std::size_t index = 0; index < kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        const FieldSpec& spec = fieldspec(id);
        runbundle::FieldEntry entry;
        entry.id = qualifiedname(id);
        entry.component = componentname(spec.component);
        entry.units = spec.units;
        entry.description = spec.description;
        entry.layout = layoutname(spec.layout);
        entry.dtype = elementname(spec.element);
        if (const auto found = provenance.fields.find(id); found != provenance.fields.end())
            entry.writtenby = found->second;

        PlanetFieldData data;
        if (source != nullptr && contains(source->carriedfields, id) && !contains(writtenfields, id))
        {
            const runbundle::ManifestEntry& carried = *findentry(source->reader->fields(), entry.id);
            data.shape = carried.shape;
            data.bytes = source->reader->fielddata(carried);
        }
        else
        {
            data = readplanetfield(world, id);
            if (data.shape.empty() && data.bytes.empty())
                continue;   // absent (a diagnostic the world has not recorded)
        }
        entry.shape = data.shape;
        entry.data = data.bytes.data();
        entry.bytes = data.bytes.size();
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
        if (source != nullptr && contains(source->carriedtables, id) && !contains(writtentables, id))
            entry.text = source->reader->tabletext(*findentry(source->reader->tables(), entry.id));
        if (const auto found = provenance.tables.find(id); found != provenance.tables.end())
            entry.writtenby = found->second;
        writer.addtable(entry);
    }

    writer.finish(metadata);
}
