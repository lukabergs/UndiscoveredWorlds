#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <memory>
#include <type_traits>

#include "core/mesh/cell_grid.hpp"

#ifdef UW_CHECK_GRID_BOUNDS
// Analysis builds: prints the index and the call stack, then aborts (world_grid_check.cpp).
[[noreturn]] void worldgridboundsfailure(int x, int y, int columns, int rows);
#endif

// A surface grid of one world, stored in the cell grid's id order
// (docs/architecture/sphere-mesh.md section 4): on the lat-lon grid row-major,
// row 0 north (cell y * columns + x); on the cube face, row, column. grid[c]
// and grid.at(c) read a cell by its mesh::CellId; the legacy grid[x][y] reads
// a lat-lon grid through a strided column proxy. It holds exactly one value per
// cell; resize() reallocates and zeroes only when the size changes.
template<typename T>
class worldgrid
{
    static_assert(std::is_trivially_copyable_v<T>, "worldgrid holds plain cell values");

public:
    worldgrid() = default;
    worldgrid(int columns, int rows) { resize(columns, rows); }
    explicit worldgrid(const mesh::CellGrid& grid) { resize(grid); }
    worldgrid(const worldgrid& other) { *this = other; }
    worldgrid(worldgrid&&) noexcept = default;
    worldgrid& operator=(worldgrid&&) noexcept = default;

    worldgrid& operator=(const worldgrid& other)
    {
        if (this == &other)
            return *this;
        if (size() != other.size())
            values_.reset(other.size() ? new T[other.size()] : nullptr);
        columns_ = other.columns_;
        rows_ = other.rows_;
        count_ = other.count_;
        if (other.size())
            std::memcpy(values_.get(), other.values_.get(), other.size() * sizeof(T));
        return *this;
    }

    // A lat-lon grid of columns x rows cells.
    void resize(int columns, int rows)
    {
        resizecells(columns, rows, static_cast<std::size_t>(columns) * static_cast<std::size_t>(rows));
    }

    // One value per cell of the grid; columns() and rows() are 0 on the cube.
    void resize(const mesh::CellGrid& grid)
    {
        if (grid.kind() == mesh::GridKind::latlon && !grid.issubgrid())
            resize(grid.columns(), grid.rows());
        else
            resizecells(0, 0, grid.cellcount());
    }

    // A cell by its id.
    T& operator[](mesh::CellId c) { return at(c); }
    const T& operator[](mesh::CellId c) const { return at(c); }
#ifdef UW_CHECK_GRID_BOUNDS
    T& at(mesh::CellId c)
    {
        if (c >= count_)
            worldgridboundsfailure(static_cast<int>(c), -1, columns_, rows_);
        return values_[c];
    }
    const T& at(mesh::CellId c) const
    {
        if (c >= count_)
            worldgridboundsfailure(static_cast<int>(c), -1, columns_, rows_);
        return values_[c];
    }
#else
    T& at(mesh::CellId c) { return values_[c]; }
    const T& at(mesh::CellId c) const { return values_[c]; }
#endif

    // The legacy grid[x][y] of a lat-lon grid: column x, read with a stride of a row.
    template<typename V>
    struct column
    {
        V* values;
        std::size_t stride;
#ifdef UW_CHECK_GRID_BOUNDS
        int x;
        const worldgrid* grid;
        V& operator[](int y) const
        {
            if (y < 0 || y >= grid->rows_)
                worldgridboundsfailure(x, y, grid->columns_, grid->rows_);
            return values[static_cast<std::size_t>(y) * stride];
        }
#else
        V& operator[](int y) const { return values[static_cast<std::size_t>(y) * stride]; }
#endif
    };

    template<typename I, std::enable_if_t<std::is_integral_v<I> && !std::is_same_v<I, mesh::CellId>, int> = 0>
    column<T> operator[](I x) { return columnat(static_cast<int>(x)); }
    template<typename I, std::enable_if_t<std::is_integral_v<I> && !std::is_same_v<I, mesh::CellId>, int> = 0>
    column<const T> operator[](I x) const { return columnat(static_cast<int>(x)); }

    void fill(T value) { std::fill(values_.get(), values_.get() + size(), value); }
    int columns() const { return columns_; }   // lat-lon only (0 on the cube)
    int rows() const { return rows_; }         // lat-lon only (0 on the cube)
    std::size_t size() const { return count_; }
    T* data() { return values_.get(); }
    const T* data() const { return values_.get(); }

    // Bitwise equality, so identical floats (including NaN payloads) compare equal.
    bool operator==(const worldgrid& other) const
    {
        return columns_ == other.columns_ && rows_ == other.rows_ && count_ == other.count_ &&
            (size() == 0 || std::memcmp(values_.get(), other.values_.get(), size() * sizeof(T)) == 0);
    }
    bool operator!=(const worldgrid& other) const { return !(*this == other); }

private:
    void resizecells(int columns, int rows, std::size_t count)
    {
        if (columns == columns_ && rows == rows_ && count == count_ && values_)
            return;
        columns_ = columns;
        rows_ = rows;
        count_ = count;
        values_.reset(count_ ? new T[count_]() : nullptr);
    }

#ifdef UW_CHECK_GRID_BOUNDS
    column<T> columnat(int x)
    {
        if (columns_ == 0 || x < 0 || x >= columns_)
            worldgridboundsfailure(x, 0, columns_, rows_);
        return { values_.get() + x, static_cast<std::size_t>(columns_), x, this };
    }
    column<const T> columnat(int x) const
    {
        if (columns_ == 0 || x < 0 || x >= columns_)
            worldgridboundsfailure(x, 0, columns_, rows_);
        return { values_.get() + x, static_cast<std::size_t>(columns_), x, this };
    }
#else
    column<T> columnat(int x) { return { values_.get() + x, static_cast<std::size_t>(columns_) }; }
    column<const T> columnat(int x) const { return { values_.get() + x, static_cast<std::size_t>(columns_) }; }
#endif

    int columns_ = 0;
    int rows_ = 0;
    std::size_t count_ = 0;
    std::unique_ptr<T[]> values_;
};

// Elevation is stored in metres as float. Code written for whole metres, the
// .uww save format, the renderers and the land/sea mask read it rounded.
inline int wholemetres(float value)
{
    return static_cast<int>(std::lround(value));
}
