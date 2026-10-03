#pragma once

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <memory>
#include <type_traits>

#ifdef UW_CHECK_GRID_BOUNDS
// Analysis builds: prints the index and the call stack, then aborts (world_grid_check.cpp).
[[noreturn]] void worldgridboundsfailure(int x, int y, int columns, int rows);
#endif

// A surface grid of one world, stored x-major (column after column) so that
// grid[x][y] reads like the former fixed [2048][1024] arrays. It holds exactly
// columns x rows cells; resize() reallocates and zeroes only when the size changes.
template<typename T>
class worldgrid
{
    static_assert(std::is_trivially_copyable_v<T>, "worldgrid holds plain cell values");

public:
    worldgrid() = default;
    worldgrid(int columns, int rows) { resize(columns, rows); }
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
        if (other.size())
            std::memcpy(values_.get(), other.values_.get(), other.size() * sizeof(T));
        return *this;
    }

    void resize(int columns, int rows)
    {
        if (columns == columns_ && rows == rows_ && values_)
            return;
        columns_ = columns;
        rows_ = rows;
        values_.reset(size() ? new T[size()]() : nullptr);
    }

#ifdef UW_CHECK_GRID_BOUNDS
    // Analysis builds (-DUW_CHECK_GRID_BOUNDS) report any index outside the grid.
    template<typename V>
    struct checkedcolumn
    {
        V* values;
        int x;
        const worldgrid* grid;
        V& operator[](int y) const
        {
            if (y < 0 || y >= grid->rows_)
                worldgridboundsfailure(x, y, grid->columns_, grid->rows_);
            return values[y];
        }
    };
    checkedcolumn<T> operator[](int x)
    {
        if (x < 0 || x >= columns_)
            worldgridboundsfailure(x, 0, columns_, rows_);
        return { values_.get() + static_cast<std::size_t>(x) * rows_, x, this };
    }
    checkedcolumn<const T> operator[](int x) const
    {
        if (x < 0 || x >= columns_)
            worldgridboundsfailure(x, 0, columns_, rows_);
        return { values_.get() + static_cast<std::size_t>(x) * rows_, x, this };
    }
#else
    T* operator[](int x) { return values_.get() + static_cast<std::size_t>(x) * rows_; }
    const T* operator[](int x) const { return values_.get() + static_cast<std::size_t>(x) * rows_; }
#endif

    void fill(T value) { std::fill(values_.get(), values_.get() + size(), value); }
    int columns() const { return columns_; }
    int rows() const { return rows_; }
    std::size_t size() const { return static_cast<std::size_t>(columns_) * static_cast<std::size_t>(rows_); }
    T* data() { return values_.get(); }
    const T* data() const { return values_.get(); }

    // Bitwise equality, so identical floats (including NaN payloads) compare equal.
    bool operator==(const worldgrid& other) const
    {
        return columns_ == other.columns_ && rows_ == other.rows_ &&
            (size() == 0 || std::memcmp(values_.get(), other.values_.get(), size() * sizeof(T)) == 0);
    }
    bool operator!=(const worldgrid& other) const { return !(*this == other); }

private:
    int columns_ = 0;
    int rows_ = 0;
    std::unique_ptr<T[]> values_;
};
