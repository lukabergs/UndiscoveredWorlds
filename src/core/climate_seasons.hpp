#pragma once

// Shared ordering of seasonal world fields and solver diagnostics. Values are
// also used by persistence; keep them stable when changing simulation modules.
constexpr int CLIMATESEASONCOUNT = 4;

enum ClimateSeason
{
    seasonjanuary = 0,
    seasonapril = 1,
    seasonjuly = 2,
    seasonoctober = 3
};
