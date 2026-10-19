#include "world/planet_export.hpp"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>

#include "world/planet.hpp"
#include "world/planet_fields.hpp"

using worldmodel::ElementType;
using worldmodel::FieldData;
using worldmodel::FieldId;
using worldmodel::Layout;
using worldmodel::RemapRule;

namespace
{
// The number of trailing axes a field on the grid has: rows and columns on lat-lon, cells on the cube.
std::size_t gridaxes(const mesh::CellGrid& grid)
{
    return grid.kind() == mesh::GridKind::latlon ? 2 : 1;
}
}

FieldRemap::FieldRemap(const mesh::CellGrid& source, const mesh::CellGrid& destination, mesh::RemapOptions options)
    : source_(source), destination_(destination), options_(options)
{
    if (source.issubgrid() || destination.issubgrid())
        throw std::invalid_argument("FieldRemap moves fields between full grids");
    identity_ = source.kind() == destination.kind() && source.columns() == destination.columns() &&
        source.rows() == destination.rows() && source.facecells() == destination.facecells() &&
        source.cellcount() == destination.cellcount();
}

FieldRemap::~FieldRemap() = default;

std::size_t FieldRemap::layers(const FieldData& data) const
{
    std::size_t count = 1;
    for (std::size_t axis = 0; axis + gridaxes(source_) < data.shape.size(); axis++)
        count *= data.shape[axis];
    return count;
}

bool FieldRemap::ongrid(const FieldData& data) const
{
    const worldmodel::FieldSpec& spec = worldmodel::fieldspec(data.id);
    if (spec.layout == Layout::scalar || spec.layout == Layout::zonalprofile || data.bytes.empty())
        return false;
    const std::vector<std::size_t>& shape = data.shape;
    if (source_.kind() == mesh::GridKind::latlon)
    {
        if (shape.size() < 2 || shape[shape.size() - 2] != static_cast<std::size_t>(source_.rows()) ||
            shape.back() != static_cast<std::size_t>(source_.columns()))
            return false;
    }
    else if (shape.empty() || shape.back() != source_.cellcount())
        return false;
    return data.bytes.size() == layers(data) * source_.cellcount() * worldmodel::elementsize(spec.element);
}

std::vector<std::size_t> FieldRemap::destinationshape(const FieldData& data) const
{
    std::vector<std::size_t> shape(data.shape.begin(), data.shape.end() - static_cast<std::ptrdiff_t>(gridaxes(source_)));
    if (destination_.kind() == mesh::GridKind::latlon)
    {
        shape.push_back(static_cast<std::size_t>(destination_.rows()));
        shape.push_back(static_cast<std::size_t>(destination_.columns()));
    }
    else
        shape.push_back(destination_.cellcount());
    return shape;
}

const mesh::RemapWeights& FieldRemap::weights() const
{
    if (!weights_)
        weights_ = std::make_unique<mesh::RemapWeights>(mesh::makeremap(source_, destination_, options_));
    return *weights_;
}

FieldData FieldRemap::zeros(const FieldData& data) const
{
    if (!ongrid(data))
        throw std::logic_error("Field is not on the remap's source grid: " + worldmodel::qualifiedname(data.id));
    FieldData out;
    out.id = data.id;
    out.shape = destinationshape(data);
    out.bytes.assign(layers(data) * destination_.cellcount() * worldmodel::elementsize(worldmodel::fieldspec(data.id).element), 0);
    return out;
}

FieldData FieldRemap::remap(const FieldData& data) const
{
    if (!ongrid(data))
        throw std::logic_error("Field is not on the remap's source grid: " + worldmodel::qualifiedname(data.id));
    const worldmodel::FieldSpec& spec = worldmodel::fieldspec(data.id);
    FieldData out;
    out.id = data.id;
    out.shape = destinationshape(data);
    if (identity_)
    {
        out.bytes = data.bytes;
        return out;
    }
    if (spec.rule == RemapRule::direction)
        throw std::logic_error("A direction field is recomputed on the destination grid, never remapped: " +
            worldmodel::qualifiedname(data.id));

    const std::size_t layercount = layers(data);
    const std::size_t sourcecells = source_.cellcount(), destinationcells = destination_.cellcount();
    const auto apply = [&](auto tag)
    {
        using T = decltype(tag);
        out.bytes.resize(layercount * destinationcells * sizeof(T));
        std::vector<T> values(sourcecells);
        for (std::size_t layer = 0; layer < layercount; layer++)
        {
            std::memcpy(values.data(), data.bytes.data() + layer * sourcecells * sizeof(T), sourcecells * sizeof(T));
            std::vector<T> result;
            switch (spec.rule)
            {
            case RemapRule::mean: result = mesh::remapmean(weights(), source_, destination_, values); break;
            case RemapRule::total: result = mesh::remaptotal(weights(), values); break;
            default: result = mesh::remapnearest(source_, destination_, values); break;
            }
            std::memcpy(out.bytes.data() + layer * destinationcells * sizeof(T), result.data(), destinationcells * sizeof(T));
        }
    };
    switch (spec.element)
    {
    case ElementType::b8:   // one byte, 0 or 1: nearest keeps it so
    case ElementType::u8: apply(std::uint8_t{}); break;
    case ElementType::i16: apply(std::int16_t{}); break;
    case ElementType::u16: apply(std::uint16_t{}); break;
    case ElementType::i32: apply(std::int32_t{}); break;
    case ElementType::f32: apply(float{}); break;
    }
    return out;
}

std::unique_ptr<planet> exportlatlon(const planet& world, int columns, const LatlonDirectionRecompute& recompute)
{
    auto latlon = std::make_unique<planet>(world);
    const mesh::CellGrid& source = world.grid();
    if (source.kind() == mesh::GridKind::latlon && source.columns() == columns)
        return latlon;
    latlon->setgrid(mesh::GridKind::latlon, columns);
    if (latlon->grid().columns() != columns)
        throw std::invalid_argument("exportlatlon: " + std::to_string(columns) + " columns is not an even width within the world's capacity");

    const FieldRemap remap(source, latlon->grid());
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        const FieldData data = readplanetfield(world, id);
        if (!remap.ongrid(data))
            continue;
        // False for what a planet does not load (the climate solver's levels and diagnostics).
        loadplanetfield(*latlon, worldmodel::fieldspec(id).rule == RemapRule::direction ? remap.zeros(data) : remap.remap(data));
    }
    // Records on cells move to the export's cell of their cell's centre; x, y are its column and row.
    const mesh::CellGrid& destination = latlon->grid();
    for (Settlement& settlement : latlon->settlements())
    {
        if (settlement.cell == mesh::nocell || settlement.cell >= source.cellcount())
            continue;
        settlement.cell = destination.locate(source.centre(settlement.cell));
        settlement.x = destination.column(settlement.cell);
        settlement.y = destination.row(settlement.cell);
    }
    if (recompute)
        recompute(*latlon);
    return latlon;
}

void importfields(planet& world, const planet& latlon, const std::vector<FieldId>& fields)
{
    const FieldRemap remap(latlon.grid(), world.grid());
    for (const FieldId id : fields)
    {
        const FieldData data = readplanetfield(latlon, id);
        const Layout layout = worldmodel::fieldspec(id).layout;
        if (layout == Layout::scalar || layout == Layout::zonalprofile)
        {
            loadplanetfield(world, data);
            continue;
        }
        if (!remap.ongrid(data))
            throw std::runtime_error("importfields: field is not on the lat-lon world's grid: " + worldmodel::qualifiedname(id));
        if (!loadplanetfield(world, remap.remap(data)))
            throw std::runtime_error("importfields: the world does not load " + worldmodel::qualifiedname(id));
    }
}
