#pragma once

// Encodes stored world fields as row-major arrays with their catalogue element
// type, for run bundles and access audits. Header-only; shared by the planet
// facade (planet_fields.cpp) and the climate module (climate_world_fields.cpp).

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "world_fields.hpp"

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

// x-major storage indexed grid[x][y] over 0..columns-1, 0..rows-1.
template<typename Grid>
FieldData grid(FieldId id, std::size_t columns, std::size_t rows, const Grid& source)
{
    using T = std::remove_cv_t<std::remove_reference_t<decltype(source[0][0])>>;
    checkelement<T>(id);
    FieldData data;
    data.id = id;
    data.shape = { rows, columns };
    data.bytes.resize(rows * columns * sizeof(T));
    auto* out = reinterpret_cast<T*>(data.bytes.data());
    for (std::size_t y = 0; y < rows; y++)
        for (std::size_t x = 0; x < columns; x++)
            out[y * columns + x] = source[x][y];
    return data;
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
}
}
