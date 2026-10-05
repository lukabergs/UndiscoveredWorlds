#pragma once

// Records which catalogue fields and tables are read. Accessors mark reads with
// UW_FIELD_READ/UW_TABLE_READ; they compile to nothing unless the build defines
// UW_TRACK_FIELD_READS (an analysis build used to derive declared process
// access). Writes are found by hashing fields around each step instead.

#include <array>
#include <atomic>
#include <cstddef>
#include <vector>

#include "core/world/world_fields.hpp"

namespace worldmodel::access
{
inline std::array<std::atomic<unsigned char>, kFieldCount> fieldreads{};
inline std::array<std::atomic<unsigned char>, kTableCount> tablereads{};
inline std::atomic<unsigned char> randomdraws{0};   // the shared legacy generator was used

inline void noterandomdraw()
{
    if (randomdraws.load(std::memory_order_relaxed) == 0)
        randomdraws.store(1, std::memory_order_relaxed);
}

inline void notefieldread(FieldId id)
{
    // Check first so that parallel readers do not keep writing a shared line.
    std::atomic<unsigned char>& flag = fieldreads[static_cast<std::size_t>(id)];
    if (flag.load(std::memory_order_relaxed) == 0)
        flag.store(1, std::memory_order_relaxed);
}

inline void notetableread(TableId id)
{
    std::atomic<unsigned char>& flag = tablereads[static_cast<std::size_t>(id)];
    if (flag.load(std::memory_order_relaxed) == 0)
        flag.store(1, std::memory_order_relaxed);
}

inline void resetreads()
{
    for (auto& flag : fieldreads)
        flag.store(0, std::memory_order_relaxed);
    for (auto& flag : tablereads)
        flag.store(0, std::memory_order_relaxed);
    randomdraws.store(0, std::memory_order_relaxed);
}

inline bool drewrandom()
{
    return randomdraws.load(std::memory_order_relaxed) != 0;
}

inline std::vector<FieldId> fieldsread()
{
    std::vector<FieldId> result;
    for (std::size_t index = 0; index < kFieldCount; index++)
        if (fieldreads[index].load(std::memory_order_relaxed))
            result.push_back(static_cast<FieldId>(index));
    return result;
}

inline std::vector<TableId> tablesread()
{
    std::vector<TableId> result;
    for (std::size_t index = 0; index < kTableCount; index++)
        if (tablereads[index].load(std::memory_order_relaxed))
            result.push_back(static_cast<TableId>(index));
    return result;
}

#ifdef UW_TRACK_FIELD_READS
inline constexpr bool tracking = true;
#else
inline constexpr bool tracking = false;
#endif
}

#ifdef UW_TRACK_FIELD_READS
#define UW_FIELD_READ(id) ::worldmodel::access::notefieldread(::worldmodel::FieldId::id)
#define UW_TABLE_READ(id) ::worldmodel::access::notetableread(::worldmodel::TableId::id)
#define UW_RANDOM_DRAW() ::worldmodel::access::noterandomdraw()
#else
#define UW_FIELD_READ(id) ((void)0)
#define UW_TABLE_READ(id) ((void)0)
#define UW_RANDOM_DRAW() ((void)0)
#endif
