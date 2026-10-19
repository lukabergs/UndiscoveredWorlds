#include "world/planet_export.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <utility>
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

// Calls apply with a value of the element's C++ type.
template <typename Apply>
void byelement(ElementType element, const Apply& apply)
{
    switch (element)
    {
    case ElementType::b8:   // one byte, 0 or 1: nearest and majority keep it so
    case ElementType::u8: apply(std::uint8_t{}); break;
    case ElementType::i16: apply(std::int16_t{}); break;
    case ElementType::u16: apply(std::uint16_t{}); break;
    case ElementType::i32: apply(std::int32_t{}); break;
    case ElementType::f32: apply(float{}); break;
    }
}

// The source cell at p; if it is of the other class, the closest cell of p's class
// among it and its neighbours two rings out (by the neighbour table); the cell at p
// when none is.
template <typename K>
mesh::CellId nearestofclass(const mesh::CellGrid& source, const mesh::Vec3& p, const std::vector<K>& sourcesea, K sea)
{
    const mesh::CellId located = source.locate(p);
    if (located == mesh::nocell || sourcesea[located] == sea)
        return located;
    std::vector<mesh::CellId> ring = { located }, cells = { located };
    for (int step = 0; step < 2; step++)
    {
        std::vector<mesh::CellId> next;
        for (const mesh::CellId c : ring)
            for (int k = 1; k <= 8; k++)
            {
                const std::int32_t n = source.neighbour(c, k);
                if (n < 0 || std::find(cells.begin(), cells.end(), static_cast<mesh::CellId>(n)) != cells.end())
                    continue;
                cells.push_back(static_cast<mesh::CellId>(n));
                next.push_back(static_cast<mesh::CellId>(n));
            }
        ring = std::move(next);
    }
    mesh::CellId best = located;
    double bestdistance = 0.0;
    for (const mesh::CellId c : cells)
    {
        if (sourcesea[c] != sea)
            continue;
        const mesh::Vec3 q = source.centre(c);
        const double distance = (q.x - p.x) * (q.x - p.x) + (q.y - p.y) * (q.y - p.y) + (q.z - p.z) * (q.z - p.z);
        if (best == located || distance < bestdistance || (distance == bestdistance && c < best))
        {
            best = c;
            bestdistance = distance;
        }
    }
    return best;
}

// remapmajority restricted to the sources of the destination cell's own class;
// where none of its class overlaps, the nearest source of its class (nearestofclass).
template <typename T>
std::vector<T> remapclass(const mesh::RemapWeights& weights, const mesh::CellGrid& source, const mesh::CellGrid& destination,
    const std::vector<T>& values, const std::vector<std::uint8_t>& sourcesea, const std::vector<std::uint8_t>& destinationsea)
{
    std::vector<T> result(weights.destinationcount);
    std::vector<std::pair<T, double>> tally;
    for (std::size_t d = 0; d < weights.destinationcount; d++)
    {
        tally.clear();
        for (std::size_t e = weights.rowstart[d]; e < weights.rowstart[d + 1]; e++)
        {
            const mesh::CellId s = weights.source[e];
            if (sourcesea[s] != destinationsea[d])
                continue;
            const T code = values[s];
            auto found = std::find_if(tally.begin(), tally.end(), [&](const auto& entry) { return entry.first == code; });
            if (found == tally.end())
                tally.emplace_back(code, weights.area_m2[e]);
            else
                found->second += weights.area_m2[e];
        }
        if (tally.empty())
        {
            const mesh::CellId s = nearestofclass(source, destination.centre(static_cast<mesh::CellId>(d)), sourcesea, destinationsea[d]);
            result[d] = s == mesh::nocell ? T{} : values[s];
            continue;
        }
        std::size_t best = 0;
        for (std::size_t n = 1; n < tally.size(); n++)
            if (tally[n].second > tally[best].second ||
                (tally[n].second == tally[best].second && tally[n].first < tally[best].first))
                best = n;
        result[d] = tally[best].first;
    }
    return result;
}

// remapmean restricted to the overlapping sources of the destination cell's own class (the same
// sums in the same order as remapmean where every overlapping source is of it); a row without
// overlap is remapmean's bilinear sample; where none of its class overlaps, the nearest source of
// its class (nearestofclass).
template <typename T>
std::vector<T> remapmeanclass(const mesh::RemapWeights& weights, const mesh::CellGrid& source, const mesh::CellGrid& destination,
    const std::vector<T>& values, const std::vector<float>& sourceclass, const std::vector<float>& destinationclass)
{
    std::vector<T> result(weights.destinationcount);
    for (std::size_t d = 0; d < weights.destinationcount; d++)
    {
        const std::size_t begin = weights.rowstart[d], end = weights.rowstart[d + 1];
        const mesh::Vec3 centre = destination.centre(static_cast<mesh::CellId>(d));
        if (begin == end)
        {
            result[d] = mesh::samplebilinear(source, values, centre);
            continue;
        }
        bool found = false;
        double first = 0.0, weighted = 0.0, total = 0.0;
        for (std::size_t e = begin; e < end; e++)
        {
            const mesh::CellId s = weights.source[e];
            if (sourceclass[s] != destinationclass[d])
                continue;
            const double value = static_cast<double>(values[s]);
            if (!found)
            {
                first = value;
                found = true;
            }
            weighted += weights.area_m2[e] * (value - first);
            total += weights.area_m2[e];
        }
        if (found)
        {
            result[d] = static_cast<T>(first + weighted / total);
            continue;
        }
        const mesh::CellId s = nearestofclass(source, centre, sourceclass, destinationclass[d]);
        result[d] = s == mesh::nocell ? T{} : values[s];
    }
    return result;
}

// Mean fields that hold on land or at sea, moved by class against the sea masks: the glacial fields.
bool seaclassmeanfield(FieldId id)
{
    return id == FieldId::glacial_denudation || id == FieldId::glaciation || id == FieldId::ice_thickness;
}

std::vector<float> asclass(const std::vector<std::uint8_t>& mask)
{
    return std::vector<float>(mask.begin(), mask.end());
}

std::vector<float> lakeclass(const planet& world)
{
    std::vector<float> lake(world.cellcount());
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
        lake[c] = world.lakesurfacefloat(c);
    return lake;
}
}

bool exportedbyclass(FieldId id)
{
    return landawarefield(id) || seaclassmeanfield(id) || id == FieldId::base_elevation;
}

bool landawarefield(FieldId id)
{
    switch (id)
    {
    case FieldId::climate_class:
    case FieldId::biome:
    case FieldId::special_features:
    case FieldId::island:
    case FieldId::mountain_island:
    case FieldId::basin_class:
    case FieldId::prograded:
    case FieldId::owner_settlement:
    case FieldId::owner_polity:
    case FieldId::sea_ice:
    case FieldId::continental_shelf:
        return true;
    default:
        return false;
    }
}

std::vector<std::uint8_t> seamask(const planet& world)
{
    std::vector<std::uint8_t> sea(world.cellcount());
    for (mesh::CellId c = 0; c < world.cellcount(); c++)
        sea[c] = world.surfacebelowsea(c) ? 1 : 0;
    return sea;
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
    byelement(spec.element, apply);
    return out;
}

FieldData FieldRemap::remapbyclass(const FieldData& data, const std::vector<std::uint8_t>& sourcesea,
    const std::vector<std::uint8_t>& destinationsea) const
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
    const std::size_t sourcecells = source_.cellcount(), destinationcells = destination_.cellcount();
    if (sourcesea.size() != sourcecells || destinationsea.size() != destinationcells)
        throw std::invalid_argument("remapbyclass: a sea mask does not match its grid");

    const std::size_t layercount = layers(data);
    const auto apply = [&](auto tag)
    {
        using T = decltype(tag);
        out.bytes.resize(layercount * destinationcells * sizeof(T));
        std::vector<T> values(sourcecells);
        for (std::size_t layer = 0; layer < layercount; layer++)
        {
            std::memcpy(values.data(), data.bytes.data() + layer * sourcecells * sizeof(T), sourcecells * sizeof(T));
            const std::vector<T> result = remapclass(weights(), source_, destination_, values, sourcesea, destinationsea);
            std::memcpy(out.bytes.data() + layer * destinationcells * sizeof(T), result.data(), destinationcells * sizeof(T));
        }
    };
    byelement(spec.element, apply);
    return out;
}

FieldData FieldRemap::remapmeanbyclass(const FieldData& data, const std::vector<float>& sourceclass,
    const std::vector<float>& destinationclass) const
{
    if (!ongrid(data))
        throw std::logic_error("Field is not on the remap's source grid: " + worldmodel::qualifiedname(data.id));
    const worldmodel::FieldSpec& spec = worldmodel::fieldspec(data.id);
    if (spec.element != ElementType::f32 || spec.rule != RemapRule::mean)
        throw std::logic_error("remapmeanbyclass moves float fields of the mean rule: " + worldmodel::qualifiedname(data.id));
    FieldData out;
    out.id = data.id;
    out.shape = destinationshape(data);
    if (identity_)
    {
        out.bytes = data.bytes;
        return out;
    }
    const std::size_t sourcecells = source_.cellcount(), destinationcells = destination_.cellcount();
    if (sourceclass.size() != sourcecells || destinationclass.size() != destinationcells)
        throw std::invalid_argument("remapmeanbyclass: a class does not match its grid");
    const std::size_t layercount = layers(data);
    out.bytes.resize(layercount * destinationcells * sizeof(float));
    std::vector<float> values(sourcecells);
    for (std::size_t layer = 0; layer < layercount; layer++)
    {
        std::memcpy(values.data(), data.bytes.data() + layer * sourcecells * sizeof(float), sourcecells * sizeof(float));
        const std::vector<float> result = remapmeanclass(weights(), source_, destination_, values, sourceclass, destinationclass);
        std::memcpy(out.bytes.data() + layer * destinationcells * sizeof(float), result.data(), destinationcells * sizeof(float));
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
        if (landawarefield(id))
            continue;
        const FieldData data = readplanetfield(world, id);
        if (!remap.ongrid(data))
            continue;
        // False for what a planet does not load (the climate solver's levels and diagnostics).
        loadplanetfield(*latlon, worldmodel::fieldspec(id).rule == RemapRule::direction ? remap.zeros(data) : remap.remap(data));
    }
    // The base elevation by lake (the lake surface, moved by nearest): a lake cell's terrain is
    // the mean over the sources of its own lake, never over its banks; a cell outside lakes the mean
    // over the sources outside lakes. Without lakes this is the mean rule.
    if (const FieldData base = readplanetfield(world, FieldId::base_elevation); remap.ongrid(base))
        loadplanetfield(*latlon, remap.remapmeanbyclass(base, lakeclass(world), lakeclass(*latlon)));
    // The land-aware codes and the glacial fields against the export's own sea mask, from its exported terrain.
    const std::vector<std::uint8_t> sourcesea = seamask(world), destinationsea = seamask(*latlon);
    const std::vector<float> sourceclass = asclass(sourcesea), destinationclass = asclass(destinationsea);
    for (std::size_t index = 0; index < worldmodel::kFieldCount; index++)
    {
        const FieldId id = static_cast<FieldId>(index);
        if (!landawarefield(id) && !seaclassmeanfield(id))
            continue;
        const FieldData data = readplanetfield(world, id);
        if (!remap.ongrid(data))
            continue;
        loadplanetfield(*latlon, landawarefield(id) ? remap.remapbyclass(data, sourcesea, destinationsea)
                                                    : remap.remapmeanbyclass(data, sourceclass, destinationclass));
    }
    // Records on cells (settlements, deposits) move to the export's cell of their cell's centre;
    // their x, y are its column and row. Polities and routes have no positions.
    const mesh::CellGrid& destination = latlon->grid();
    const auto relocate = [&](mesh::CellId& cell, int& x, int& y)
    {
        if (cell == mesh::nocell || cell >= source.cellcount())
            return;
        cell = destination.locate(source.centre(cell));
        x = destination.column(cell);
        y = destination.row(cell);
    };
    for (Settlement& settlement : latlon->settlements())
        relocate(settlement.cell, settlement.x, settlement.y);
    for (mineralsystems::Deposit& deposit : latlon->mineraldeposits())
        relocate(deposit.surface.cell, deposit.surface.worldX, deposit.surface.worldY);
    if (recompute)
        recompute(*latlon);
    return latlon;
}

void importfields(planet& world, const planet& latlon, const std::vector<FieldId>& fields)
{
    const FieldRemap remap(latlon.grid(), world.grid());
    std::vector<std::uint8_t> sourcesea, destinationsea;   // for land-aware codes, built on first use
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
        if (landawarefield(id) && !remap.identity() && sourcesea.empty())
        {
            sourcesea = seamask(latlon);
            destinationsea = seamask(world);
        }
        if (!loadplanetfield(world, landawarefield(id) ? remap.remapbyclass(data, sourcesea, destinationsea) : remap.remap(data)))
            throw std::runtime_error("importfields: the world does not load " + worldmodel::qualifiedname(id));
    }
}
