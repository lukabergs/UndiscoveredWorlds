#pragma once

#include "climate_reference.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>

namespace climateio::experiment
{
// Explicit benchmark-only inputs. Reuses UWCLIM1; never changes saved worlds.
class Inputs
{
public:
    explicit Inputs(const std::filesystem::path& directory = {}) : directory_(directory)
    {
        if (directory.empty()) return;
        if (!std::filesystem::is_directory(directory))
            throw std::runtime_error("Reference experiment directory does not exist");
        for (const auto& size : std::filesystem::directory_iterator(directory))
        {
            if (!size.is_directory()) continue;
            const auto label = size.path().filename().string();
            if (label.empty() || label.find_first_not_of("0123456789") != std::string::npos) continue;
            const int width = std::stoi(label);
            for (const auto& entry : std::filesystem::directory_iterator(size.path()))
            {
                if (entry.path().extension() != ".uwclim") continue;
                const auto role = entry.path().stem().string();
                climatereference::MonthlyGrid field;
                std::string error;
                if (!climatereference::loadMonthlyGrid(entry.path(), role, field, &error) ||
                    field.width != width || field.monthCount != 12 ||
                    std::filesystem::file_size(entry.path()) != 40 + field.values.size() * sizeof(float))
                    throw std::runtime_error("Invalid reference experiment input: " + entry.path().string() + ": " + error);
                for (const float value : field.values)
                    if (std::isinf(value)) throw std::runtime_error("Infinite reference experiment input");
                fields_[width].emplace(role, std::move(field));
            }
        }
        if (fields_.empty()) throw std::runtime_error("Reference experiment has no fields");
    }

    const climatereference::MonthlyGrid* field(int width, const std::string& role) const
    {
        const auto size = fields_.find(width);
        if (size == fields_.end()) return nullptr;
        const auto value = size->second.find(role);
        return value == size->second.end() ? nullptr : &value->second;
    }
    const std::filesystem::path& directory() const { return directory_; }

private:
    std::filesystem::path directory_;
    std::map<int, std::map<std::string, climatereference::MonthlyGrid>> fields_;
};

inline float sample(const climatereference::MonthlyGrid* field, int firstSeason, int secondSeason,
    float fraction, int x, int y, float fallback)
{
    if (!field || firstSeason < 0 || firstSeason > 3 || secondSeason < 0 || secondSeason > 3)
        return fallback;
    const float first = field->value(firstSeason * 3, x, y);
    const float second = field->value(secondSeason * 3, x, y);
    const auto valid = [](float value) { return std::isfinite(value) && value != -9999.9f; };
    const float alpha = std::clamp(fraction, 0.0f, 1.0f);
    if (alpha == 0.0f) return valid(first) ? first : fallback;
    if (alpha == 1.0f) return valid(second) ? second : fallback;
    if (!valid(first) || !valid(second)) return fallback;
    return first + alpha * (second - first);
}
}
