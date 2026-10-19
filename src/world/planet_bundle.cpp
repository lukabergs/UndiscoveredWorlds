#include "world/planet_bundle.hpp"

#include <algorithm>
#include <cstdint>
#include <memory>
#include <stdexcept>

#include "core/mesh/cell_grid.hpp"
#include "world/planet.hpp"
#include "world/planet_export.hpp"
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

template<typename T>
void addcells(runbundle::Writer& writer, const char* id, const char* dtype, std::vector<std::size_t> shape, const std::vector<T>& values)
{
    runbundle::ArrayEntry entry;
    entry.id = id;
    entry.dtype = dtype;
    entry.shape = std::move(shape);
    entry.data = values.data();
    entry.bytes = values.size() * sizeof(T);
    writer.addcells(entry);
}

// The cell table of a cube grid (sphere-mesh.md section 7), in id order.
void writecelltable(runbundle::Writer& writer, const mesh::CellGrid& grid)
{
    const std::size_t count = grid.cellcount();
    constexpr double degrees = 180.0 / 3.14159265358979323846;
    std::vector<double> centre(3 * count), latitude(count), longitude(count), area(count);
    std::vector<std::int32_t> neighbours(mesh::slotcount * count);
    std::vector<std::int8_t> reverse(mesh::slotcount * count);
    for (mesh::CellId c = 0; c < count; c++)
    {
        const mesh::Vec3 p = grid.centre(c);
        centre[3 * c] = p.x;
        centre[3 * c + 1] = p.y;
        centre[3 * c + 2] = p.z;
        latitude[c] = grid.latitude_rad(c) * degrees;
        longitude[c] = grid.longitude_rad(c) * degrees;
        area[c] = grid.area_m2(c);
        for (int k = 1; k <= mesh::slotcount; k++)
        {
            neighbours[mesh::slotcount * c + (k - 1)] = grid.neighbour(c, k);
            reverse[mesh::slotcount * c + (k - 1)] = static_cast<std::int8_t>(grid.reverse(c, k));   // 0 where no neighbour
        }
    }
    const std::size_t slots = static_cast<std::size_t>(mesh::slotcount);
    addcells(writer, "centre", "<f8", { count, 3 }, centre);
    addcells(writer, "latitude_deg", "<f8", { count }, latitude);
    addcells(writer, "longitude_deg", "<f8", { count }, longitude);
    addcells(writer, "area_m2", "<f8", { count }, area);
    addcells(writer, "neighbours", "<i4", { count, slots }, neighbours);
    addcells(writer, "reverse", "|i1", { count, slots }, reverse);
}
}

PlanetBundleSource loadplanetbundle(planet& world, const std::filesystem::path& directory)
{
    using namespace worldmodel;
    PlanetBundleSource source;
    source.reader = std::make_shared<runbundle::Reader>(directory);
    const runbundle::Reader& reader = *source.reader;
    const bool cube = reader.grid().kind == "cube";
    const std::size_t facecells = cube ? static_cast<std::size_t>(reader.grid().facecells) : 0;
    if (cube)
    {
        // A cube bundle gives a cube world (sphere-mesh.md section 7).
        if (reader.grid().facecells < 1 || reader.columns() != 4 * reader.grid().facecells)
            throw std::runtime_error("Run bundle cube grid has no N >= 1 with export columns 4 N: " + directory.string());
        world.setgrid(mesh::GridKind::cube, reader.columns());
    }
    else if (reader.grid().kind == "latlon")
    {
        if (reader.columns() < 2 || reader.columns() != reader.rows() * 2)
            throw std::runtime_error("Run bundle grid is not W x W/2: " + directory.string());
        world.setwidth(reader.columns() - 1);
        if (world.width() + 1 != reader.columns())
            throw std::runtime_error("Run bundle grid exceeds the world's capacity: " + directory.string());
    }
    else
        throw std::runtime_error("Run bundle grid kind " + reader.grid().kind + " does not load as a world: " + directory.string());

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
        if (cube)
            encoding::cubeworldshape(data.shape, facecells);   // [layer,] 6, N, N -> [layer,] cells
        data.bytes = reader.fielddata(entry);
        // A cube bundle's zonal profiles are climate's lat-lon ones (its shadow's): carried.
        if (cube && fieldspec(data.id).layout == Layout::zonalprofile)
            source.carriedfields.push_back(data.id);
        else if (loadplanetfield(world, data))
        {
            source.loadedfields.push_back(data.id);
            // The radius is the grid's too: its areas and distances follow it.
            if (data.id == FieldId::radius)
            {
                float radius = 0.0f;
                encoding::loadscalar(data, radius);
                world.setradius(radius);
            }
        }
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
    const std::vector<worldmodel::TableId>& writtentables, const PlanetBundleExport& exports)
{
    using namespace worldmodel;
    runbundle::Writer writer(directory);

    // The grid block is the world's (sphere-mesh.md section 7).
    const mesh::CellGrid& grid = world.grid();
    const bool cube = grid.kind() == mesh::GridKind::cube;
    const std::size_t facecells = cube ? static_cast<std::size_t>(grid.facecells()) : 0;
    runbundle::Metadata manifest = metadata;
    manifest.grid.kind = cube ? "cube" : "latlon";
    manifest.grid.facecells = cube ? grid.facecells() : 0;
    manifest.grid.radiuskm = grid.radius_m() / 1000.0;
    manifest.columns = cube ? 4 * grid.facecells() : grid.columns();
    manifest.rows = manifest.columns / 2;

    // A cube bundle carries the lat-lon export of every field on the world grid (section 6):
    // the recomputed fields from the world's export, the climate shadow's own from the shadow.
    std::shared_ptr<const mesh::CellGrid> exportgrid;
    std::unique_ptr<FieldRemap> exporter;
    std::unique_ptr<planet> exported;
    const planet* shadow = nullptr;
    std::unique_ptr<FieldRemap> shadowshape;   // an identity on the shadow's grid: checks a field is on it
    if (cube)
    {
        exportgrid = mesh::makelatlongrid(manifest.columns, grid.radius_m());
        exporter = std::make_unique<FieldRemap>(grid, *exportgrid);
        writecelltable(writer, grid);
        if (exports.recompute)
            exported = exportlatlon(world, manifest.columns, exports.recompute);
        if (exports.climateshadow != nullptr && exports.shadowowned)
        {
            shadow = exports.climateshadow;
            if (shadow->grid().kind() != mesh::GridKind::latlon || shadow->grid().columns() != manifest.columns)
                throw std::logic_error("writeplanetbundle: the climate shadow is not on the export's grid");
            shadowshape = std::make_unique<FieldRemap>(shadow->grid(), shadow->grid());
        }
    }

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
        const bool carried = source != nullptr && contains(source->carriedfields, id) && !contains(writtenfields, id);

        // A field the climate shadow owns (and the source does not carry) is the shadow's:
        // exported from it on its grid, else stored from it in its lat-lon shape.
        PlanetFieldData shadowdata;
        bool fromshadow = false, shadowexport = false;
        if (shadow != nullptr && !carried && exports.shadowowned(id))
        {
            shadowdata = readplanetfield(*shadow, id);
            fromshadow = !(shadowdata.shape.empty() && shadowdata.bytes.empty());
            shadowexport = fromshadow && shadowshape->ongrid(shadowdata);
        }

        // data.shape is the bundle's, nativeshape the planet's encoding ([layer,] cells on the cube).
        PlanetFieldData data;
        std::vector<std::size_t> nativeshape;
        bool native = true;
        if (carried)
        {
            const runbundle::ManifestEntry& carriedentry = *findentry(source->reader->fields(), entry.id);
            data.id = id;
            data.shape = carriedentry.shape;
            data.bytes = source->reader->fielddata(carriedentry);
            nativeshape = data.shape;
            if (cube)
                encoding::cubeworldshape(nativeshape, facecells);
        }
        else if (fromshadow && !shadowexport)
        {
            data = std::move(shadowdata);
            nativeshape = data.shape;
        }
        else
        {
            data = readplanetfield(world, id);
            native = !(data.shape.empty() && data.bytes.empty());   // absent: a diagnostic the world has not recorded
            if (!native && !shadowexport)
                continue;
            nativeshape = data.shape;
            if (cube && native)
                encoding::cubebundleshape(data.shape, facecells);
        }

        if (exporter && !(fromshadow && !shadowexport))
        {
            runbundle::FieldEntry exportentry = entry;
            if (shadowexport)
            {
                // Climate's own values on the export's grid; its direction codes are lat-lon codes.
                exportentry.shape = shadowdata.shape;
                exportentry.data = shadowdata.bytes.data();
                exportentry.bytes = shadowdata.bytes.size();
                writer.addexport(exportentry, remaprulename(spec.rule), spec.rule == RemapRule::direction);
            }
            else
            {
                std::vector<std::size_t> bundleshape = data.shape;
                data.shape = nativeshape;
                if (exporter->ongrid(data))
                {
                    PlanetFieldData out;
                    bool recomputed = false;
                    if (exported != nullptr && !carried && contains(exports.recomputed, id))
                    {
                        out = readplanetfield(*exported, id);
                        recomputed = true;
                    }
                    else
                        out = spec.rule == RemapRule::direction ? exporter->zeros(data) : exporter->remap(data);
                    exportentry.shape = out.shape;
                    exportentry.data = out.bytes.data();
                    exportentry.bytes = out.bytes.size();
                    // Without the routing's callback a direction field is exported as 0.
                    writer.addexport(exportentry, remaprulename(spec.rule), recomputed);
                }
                data.shape = std::move(bundleshape);
            }
        }

        if (!native)
            continue;
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
        // The climate tables describe the shadow's levels and weather on a cube world.
        const bool climatetable = id == TableId::atmosphere_levels || id == TableId::weather_state;
        entry.text = planettabletext(shadow != nullptr && climatetable ? *shadow : world, id, &entry.extension);
        if (source != nullptr && contains(source->carriedtables, id) && !contains(writtentables, id))
            entry.text = source->reader->tabletext(*findentry(source->reader->tables(), entry.id));
        if (const auto found = provenance.tables.find(id); found != provenance.tables.end())
            entry.writtenby = found->second;
        writer.addtable(entry);
    }

    writer.finish(manifest);
}
