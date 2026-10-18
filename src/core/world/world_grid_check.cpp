// Bounds failures of worldgrid in analysis builds (-DUW_CHECK_GRID_BOUNDS):
// an (x, y) outside the lat-lon grid, (x, y) on a cube grid, or a cell id beyond the grid.
#ifdef UW_CHECK_GRID_BOUNDS
#include "core/world/world_grid.hpp"

#include <cstdio>
#include <cstdlib>

#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <dbghelp.h>
#pragma comment(lib, "dbghelp.lib")

void worldgridboundsfailure(int x, int y, int columns, int rows)
{
    // y is -1 for a cell id (x) beyond the cell count; columns 0 is a cube grid read by [x][y].
    if (y == -1)
        std::fprintf(stderr, "worldgrid cell %d outside the grid (%d x %d)\n", x, columns, rows);
    else
        std::fprintf(stderr, "worldgrid index (%d, %d) outside %d x %d\n", x, y, columns, rows);
    void* frames[48];
    const USHORT count = CaptureStackBackTrace(1, 48, frames, nullptr);
    HANDLE process = GetCurrentProcess();
    SymSetOptions(SYMOPT_LOAD_LINES | SYMOPT_UNDNAME);
    SymInitialize(process, nullptr, TRUE);
    alignas(SYMBOL_INFO) char buffer[sizeof(SYMBOL_INFO) + 512];
    auto* symbol = reinterpret_cast<SYMBOL_INFO*>(buffer);
    for (USHORT index = 0; index < count; index++)
    {
        symbol->SizeOfStruct = sizeof(SYMBOL_INFO);
        symbol->MaxNameLen = 511;
        const DWORD64 address = reinterpret_cast<DWORD64>(frames[index]);
        IMAGEHLP_LINE64 line{};
        line.SizeOfStruct = sizeof(line);
        DWORD offset = 0;
        const bool named = SymFromAddr(process, address, nullptr, symbol) != FALSE;
        const bool located = SymGetLineFromAddr64(process, address, &offset, &line) != FALSE;
        std::fprintf(stderr, "  %s  %s:%lu\n", named ? symbol->Name : "?", located ? line.FileName : "?",
                     located ? static_cast<unsigned long>(line.LineNumber) : 0ul);
    }
    std::fflush(stderr);
    std::abort();
}
#endif
