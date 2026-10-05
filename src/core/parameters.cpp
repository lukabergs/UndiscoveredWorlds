#include "core/parameters.hpp"

#include <charconv>
#include <limits>

namespace parameters
{
namespace
{
template <typename T>
bool parsenumber(std::string_view text, T& value)
{
    T parsed{};
    const char* last = text.data() + text.size();
    const std::from_chars_result result = std::from_chars(text.data(), last, parsed);
    if (result.ec != std::errc() || result.ptr != last)
        return false;
    value = parsed;
    return true;
}

template <typename T>
std::string numbertext(T value)
{
    char buffer[64];
    const std::to_chars_result result = std::to_chars(buffer, buffer + sizeof(buffer), value);
    return std::string(buffer, result.ptr);
}

struct Setter
{
    std::string_view text;

    bool operator()(bool* value) const
    {
        if (text == "true" || text == "on")
            *value = true;
        else if (text == "false" || text == "off")
            *value = false;
        else
            return false;
        return true;
    }

    bool operator()(short* value) const
    {
        int parsed = 0;
        if (!parsenumber(text, parsed) || parsed < std::numeric_limits<short>::min() || parsed > std::numeric_limits<short>::max())
            return false;
        *value = static_cast<short>(parsed);
        return true;
    }

    bool operator()(int* value) const { return parsenumber(text, *value); }
    bool operator()(float* value) const { return parsenumber(text, *value); }
    bool operator()(double* value) const { return parsenumber(text, *value); }
};

struct Printer
{
    std::string operator()(const bool* value) const { return *value ? "true" : "false"; }
    std::string operator()(const short* value) const { return numbertext(static_cast<int>(*value)); }
    std::string operator()(const int* value) const { return numbertext(*value); }
    std::string operator()(const float* value) const { return numbertext(*value); }
    std::string operator()(const double* value) const { return numbertext(*value); }
};
}

bool settext(const Entry& entry, std::string_view text, std::string* error)
{
    if (std::visit(Setter{ text }, entry.value))
        return true;
    if (error != nullptr)
        *error = "'" + std::string(text) + "' is not a value for " + entry.name;
    return false;
}

std::string text(const Entry& entry)
{
    return std::visit(Printer{}, entry.value);
}
}
