#pragma once

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

// Appends one row of stage timings to the LOG sheet of an existing Excel
// profiling workbook (Windows, through PowerShell and Excel automation).
// Does nothing when the workbook does not exist.
void appendprofilingworkbookrow(const std::filesystem::path& workbookpath, long seed,
    const std::vector<std::pair<std::string, double>>& timings);
