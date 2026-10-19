#pragma once

// Encodes stored world fields as row-major arrays with their catalogue element
// type, for run bundles and access audits, and loads them back from a bundle. Header-only; shared by the planet
// facade (planet_fields.cpp) and the climate module (climate_world_fields.cpp).

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "core/world/world_fields.hpp"
#include "core/world/world_grid.hpp"

namespace worldmodel
{
struct FieldData
{
    FieldId id{};
    std::vector<std::size_t> shape;      // row-major: [layer,] row, column; empty for scalars
    std::vector<unsigned char> bytes;    // little-endian elements of the field's type
};

namespace encoding
{
template<typename T>
constexpr bool matcheselement(ElementType element)
{
    switch (element)
    {
    case ElementType::b8: return std::is_same_v<T, bool>;
    case ElementType::u8: return std::is_same_v<T, std::uint8_t>;
    case ElementType::i16: return std::is_same_v<T, short>;
    case ElementType::u16: return std::is_same_v<T, unsigned short>;
    case ElementType::i32: return std::is_same_v<T, int>;
    case ElementType::f32: return std::is_same_v<T, float>;
    }
    return false;
}

template<typename T>
void checkelement(FieldId id)
{
    if (!matcheselement<T>(fieldspec(id).element))
        throw std::logic_error("World field element type mismatch: " + qualifiedname(id));
}

// A world grid is stored in its cell grid's id order (row-major on lat-lon): its
// bytes are the field's. Shape [row, column] on lat-lon, [cell] otherwise (bundle
// v2 gives the cube its own shape, docs/architecture/sphere-mesh.md section 7:
// cubebundleshape below). Every grid field's storage is a worldgrid (the former
// grid[x][y] path had no field left, T8 W2-F); on a cube columns and rows are 0.
template<typename T>
FieldData grid(FieldId id, std::size_t columns, std::size_t rows, const worldgrid<T>& source)
{
    checkelement<T>(id);
    FieldData data;
    data.id = id;
    const bool latlon = source.columns() > 0;
    if (latlon && (static_cast<std::size_t>(source.columns()) != columns || static_cast<std::size_t>(source.rows()) != rows))
        throw std::logic_error("World field does not match the world grid: " + qualifiedname(id));
    data.shape = latlon ? std::vector<std::size_t>{ rows, columns } : std::vector<std::size_t>{ source.size() };
    data.bytes.resize(source.size() * sizeof(T));
    if (source.size())
        std::memcpy(data.bytes.data(), source.data(), source.size() * sizeof(T));
    return data;
}

// Bundle v2 stores a cube field as [layer,] 6, N, N (id order: face, row, column);
// the encoding above gives [layer,] cells. cubebundleshape turns a shape ending in
// the 6 N^2 cells into the bundle's, cubeworldshape the reverse; both return false
// and leave the shape unchanged when it does not end in the cube's cells (scalars,
// zonal profiles, the climate solver's own grids).
inline bool cubebundleshape(std::vector<std::size_t>& shape, std::size_t facecells)
{
    if (shape.empty() || facecells == 0 || shape.back() != 6 * facecells * facecells)
        return false;
    shape.back() = 6;
    shape.push_back(facecells);
    shape.push_back(facecells);
    return true;
}

inline bool cubeworldshape(std::vector<std::size_t>& shape, std::size_t facecells)
{
    const std::size_t count = shape.size();
    if (count < 3 || facecells == 0 || shape[count - 3] != 6 || shape[count - 2] != facecells || shape[count - 1] != facecells)
        return false;
    shape.resize(count - 2);
    shape.back() = 6 * facecells * facecells;
    return true;
}

// Seasonal layers stored row-major per layer. Unsized storage keeps its size.
template<typename Layers>
FieldData seasonal(FieldId id, std::size_t columns, std::size_t rows, const Layers& layers)
{
    using T = typename std::remove_reference_t<decltype(layers[0])>::value_type;
    checkelement<T>(id);
    FieldData data;
    data.id = id;
    const std::size_t cells = layers[0].size();
    for (const auto& layer : layers)
        if (layer.size() != cells)
            throw std::logic_error("Seasonal layers differ in size: " + qualifiedname(id));
    if (cells == rows * columns)
        data.shape = { layers.size(), rows, columns };
    else
        data.shape = { layers.size(), cells };
    data.bytes.resize(layers.size() * cells * sizeof(T));
    auto* out = reinterpret_cast<T*>(data.bytes.data());
    for (std::size_t layer = 0; layer < layers.size(); layer++)
        for (std::size_t cell = 0; cell < cells; cell++)
            out[layer * cells + cell] = layers[layer][cell];
    return data;
}

// Solver grids: x-major vector<vector<float>> per layer; shape follows the data.
template<typename Layers>
FieldData solverseasonal(FieldId id, const Layers& layers)
{
    checkelement<float>(id);
    FieldData data;
    data.id = id;
    const std::size_t layercolumns = layers[0].size();
    const std::size_t layerrows = layercolumns ? layers[0][0].size() : 0;
    data.shape = { layers.size(), layerrows, layercolumns };
    data.bytes.resize(layers.size() * layerrows * layercolumns * sizeof(float));
    auto* out = reinterpret_cast<float*>(data.bytes.data());
    for (std::size_t layer = 0; layer < layers.size(); layer++)
    {
        if (layers[layer].size() != layercolumns)
            throw std::logic_error("Solver layers differ in size: " + qualifiedname(id));
        for (std::size_t x = 0; x < layercolumns; x++)
        {
            if (layers[layer][x].size() != layerrows)
                throw std::logic_error("Solver columns differ in size: " + qualifiedname(id));
            for (std::size_t y = 0; y < layerrows; y++)
                out[(layer * layerrows + y) * layercolumns + x] = layers[layer][x][y];
        }
    }
    return data;
}

// Monthly layers stored row-major on a columns x rows solver grid.
template<typename Layers>
FieldData monthly(FieldId id, std::size_t columns, std::size_t rows, const Layers& layers)
{
    checkelement<float>(id);
    FieldData data;
    data.id = id;
    const std::size_t cells = columns * rows;
    const bool complete = cells > 0 && layers[0].size() == cells;
    data.shape = { layers.size(), complete ? rows : 0, complete ? columns : 0 };
    if (!complete)
        return data;
    data.bytes.resize(layers.size() * cells * sizeof(float));
    auto* out = reinterpret_cast<float*>(data.bytes.data());
    for (std::size_t layer = 0; layer < layers.size(); layer++)
    {
        if (layers[layer].size() != cells)
            throw std::logic_error("Monthly layers differ in size: " + qualifiedname(id));
        for (std::size_t cell = 0; cell < cells; cell++)
            out[layer * cells + cell] = layers[layer][cell];
    }
    return data;
}

// Six bands per column, stored grid[x][band].
template<typename Profile>
FieldData zonalprofile(FieldId id, std::size_t columns, const Profile& source)
{
    using T = std::remove_cv_t<std::remove_reference_t<decltype(source[0][0])>>;
    checkelement<T>(id);
    FieldData data;
    data.id = id;
    data.shape = { 6, columns };
    data.bytes.resize(6 * columns * sizeof(T));
    auto* out = reinterpret_cast<T*>(data.bytes.data());
    for (std::size_t band = 0; band < 6; band++)
        for (std::size_t x = 0; x < columns; x++)
            out[band * columns + x] = source[x][band];
    return data;
}

template<typename Out, typename In>
void appendvalue(std::vector<unsigned char>& bytes, In value)
{
    const Out converted = static_cast<Out>(value);
    const auto* raw = reinterpret_cast<const unsigned char*>(&converted);
    bytes.insert(bytes.end(), raw, raw + sizeof(Out));
}

// Scalars are converted to the catalogue element type.
template<typename In>
FieldData scalar(FieldId id, In value)
{
    FieldData data;
    data.id = id;
    switch (fieldspec(id).element)
    {
    case ElementType::b8: appendvalue<bool>(data.bytes, value); break;
    case ElementType::i32: appendvalue<int>(data.bytes, value); break;
    case ElementType::f32: appendvalue<float>(data.bytes, value); break;
    default: throw std::logic_error("Unsupported scalar type: " + qualifiedname(id));
    }
    return data;
}

// Loading is the inverse of encoding: row-major data back into the storage. Data
// whose type or shape does not match the storage is an error, never a partial load.
inline void checkloadsize(const FieldData& data, std::size_t elements, std::size_t elementbytes)
{
    if (data.bytes.size() != elements * elementbytes)
        throw std::runtime_error("Field data does not match its storage: " + qualifiedname(data.id));
}

template<typename T>
void loadgrid(const FieldData& data, std::size_t columns, std::size_t rows, worldgrid<T>& target)
{
    checkelement<T>(data.id);
    const bool latlon = target.columns() > 0;
    if (latlon && (static_cast<std::size_t>(target.columns()) != columns || static_cast<std::size_t>(target.rows()) != rows))
        throw std::runtime_error("Field shape does not match the world grid: " + qualifiedname(data.id));
    if (data.shape != (latlon ? std::vector<std::size_t>{ rows, columns } : std::vector<std::size_t>{ target.size() }))
        throw std::runtime_error("Field shape does not match the world grid: " + qualifiedname(data.id));
    checkloadsize(data, target.size(), sizeof(T));
    if (target.size())
        std::memcpy(target.data(), data.bytes.data(), target.size() * sizeof(T));
}

template<typename Layers>
void loadseasonal(const FieldData& data, Layers& layers)
{
    using T = typename std::remove_reference_t<decltype(layers[0])>::value_type;
    checkelement<T>(data.id);
    if (data.shape.size() < 2 || data.shape[0] != layers.size())
        throw std::runtime_error("Field shape does not match its seasons: " + qualifiedname(data.id));
    std::size_t cells = 1;
    for (std::size_t axis = 1; axis < data.shape.size(); axis++)
        cells *= data.shape[axis];
    checkloadsize(data, layers.size() * cells, sizeof(T));
    const auto* in = reinterpret_cast<const T*>(data.bytes.data());
    for (std::size_t layer = 0; layer < layers.size(); layer++)
        layers[layer].assign(in + layer * cells, in + (layer + 1) * cells);
}

template<typename Profile>
void loadzonalprofile(const FieldData& data, std::size_t columns, Profile& target)
{
    using T = std::remove_reference_t<decltype(target[0][0])>;
    checkelement<T>(data.id);
    if (data.shape != std::vector<std::size_t>{ 6, columns })
        throw std::runtime_error("Field shape does not match the world grid: " + qualifiedname(data.id));
    checkloadsize(data, 6 * columns, sizeof(T));
    const auto* in = reinterpret_cast<const T*>(data.bytes.data());
    for (std::size_t band = 0; band < 6; band++)
        for (std::size_t x = 0; x < columns; x++)
            target[x][band] = in[band * columns + x];
}

template<typename In, typename Out>
void loadvalue(const FieldData& data, Out& value)
{
    checkloadsize(data, 1, sizeof(In));
    In stored;
    std::memcpy(&stored, data.bytes.data(), sizeof(In));
    value = static_cast<Out>(stored);
}

template<typename Out>
void loadscalar(const FieldData& data, Out& value)
{
    switch (fieldspec(data.id).element)
    {
    case ElementType::b8: loadvalue<bool>(data, value); break;
    case ElementType::i32: loadvalue<int>(data, value); break;
    case ElementType::f32: loadvalue<float>(data, value); break;
    default: throw std::logic_error("Unsupported scalar type: " + qualifiedname(data.id));
    }
}

// One mapping from field ids to storage serves both directions: the mapping
// calls scalar, grid, seasonal or zonalprofile on the operation with the storage.
struct ReadOperation
{
    std::size_t columns = 0, rows = 0;
    FieldData data;
    template<typename T> bool scalar(FieldId id, const T& value) { data = encoding::scalar(id, value); return true; }
    template<typename Grid> bool grid(FieldId id, const Grid& source) { data = encoding::grid(id, columns, rows, source); return true; }
    template<typename Layers> bool seasonal(FieldId id, const Layers& layers) { data = encoding::seasonal(id, columns, rows, layers); return true; }
    template<typename Profile> bool zonalprofile(FieldId id, const Profile& source) { data = encoding::zonalprofile(id, columns, source); return true; }
};

struct LoadOperation
{
    std::size_t columns = 0, rows = 0;
    const FieldData* data = nullptr;
    template<typename T> bool scalar(FieldId, T& value) { loadscalar(*data, value); return true; }
    template<typename Grid> bool grid(FieldId, Grid& target) { loadgrid(*data, columns, rows, target); return true; }
    template<typename Layers> bool seasonal(FieldId, Layers& layers) { loadseasonal(*data, layers); return true; }
    template<typename Profile> bool zonalprofile(FieldId, Profile& target) { loadzonalprofile(*data, columns, target); return true; }
};
}
}
