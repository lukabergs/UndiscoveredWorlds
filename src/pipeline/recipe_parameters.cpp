#include "pipeline/recipe_parameters.hpp"

#include <cctype>
#include <fstream>
#include <sstream>

#include "simulations/climate/climate_parameters.hpp"
#include "simulations/landscape/landscape_parameters.hpp"
#include "simulations/terrain/terrain_parameters.hpp"
#include "world/world_parameters.hpp"

namespace
{
const parameters::Entry* findentry(std::string_view name)
{
    for (const parameters::Entry& entry : recipeparameters())
    {
        if (name == entry.name)
            return &entry;
    }
    return nullptr;
}

bool setnamed(std::string_view name, std::string_view text, std::string* error)
{
    const parameters::Entry* entry = findentry(name);
    if (entry == nullptr)
    {
        if (error != nullptr)
            *error = "unknown parameter " + std::string(name);
        return false;
    }
    return parameters::settext(*entry, text, error);
}

// Reads nested JSON objects of numbers and booleans; a value's name is its keys joined by dots.
struct ValueReader
{
    const std::string& text;
    size_t position = 0;
    std::string* error;

    void skipspace()
    {
        while (position < text.size() && std::isspace(static_cast<unsigned char>(text[position])))
            position++;
    }

    bool fail(const std::string& message)
    {
        if (error != nullptr && error->empty())
            *error = message + " at character " + std::to_string(position);
        return false;
    }

    bool readobject(const std::string& prefix)
    {
        skipspace();
        if (position >= text.size() || text[position] != '{')
            return fail("expected an object");
        position++;
        skipspace();
        if (position < text.size() && text[position] == '}')
        {
            position++;
            return true;
        }
        while (true)
        {
            skipspace();
            if (position >= text.size() || text[position] != '"')
                return fail("expected a key");
            const size_t end = text.find('"', position + 1);
            if (end == std::string::npos)
                return fail("unterminated key");
            const std::string name = prefix + text.substr(position + 1, end - position - 1);
            position = end + 1;
            skipspace();
            if (position >= text.size() || text[position] != ':')
                return fail("expected ':'");
            position++;
            skipspace();
            if (position < text.size() && text[position] == '{')
            {
                if (!readobject(name + "."))
                    return false;
            }
            else
            {
                const size_t start = position;
                while (position < text.size() && text[position] != ',' && text[position] != '}' && !std::isspace(static_cast<unsigned char>(text[position])))
                    position++;
                if (!setnamed(name, std::string_view(text).substr(start, position - start), error))
                    return false;
            }
            skipspace();
            if (position < text.size() && text[position] == ',')
            {
                position++;
                continue;
            }
            if (position < text.size() && text[position] == '}')
            {
                position++;
                return true;
            }
            return fail("expected ',' or '}'");
        }
    }
};
}

const std::vector<parameters::Entry>& recipeparameters()
{
    static const std::vector<parameters::Entry> entries = []
    {
        std::vector<parameters::Entry> all;
        for (const std::vector<parameters::Entry>* module : { &tuning::world::parameterentries(), &tuning::terrain::parameterentries(), &tuning::climate::parameterentries(), &tuning::landscape::parameterentries() })
            all.insert(all.end(), module->begin(), module->end());
        return all;
    }();
    return entries;
}

bool setrecipeparameter(std::string_view assignment, std::string* error)
{
    const size_t equals = assignment.find('=');
    if (equals == std::string_view::npos || equals == 0)
    {
        if (error != nullptr)
            *error = "expected name=value, not '" + std::string(assignment) + "'";
        return false;
    }
    return setnamed(assignment.substr(0, equals), assignment.substr(equals + 1), error);
}

bool getrecipeparameter(std::string_view name, std::string& value)
{
    const parameters::Entry* entry = findentry(name);
    if (entry == nullptr)
        return false;
    value = parameters::text(*entry);
    return true;
}

bool checkrecipeparameter(std::string_view name, std::string_view value, std::string* error)
{
    const parameters::Entry* entry = findentry(name);
    if (entry == nullptr)
    {
        if (error != nullptr)
            *error = "unknown parameter " + std::string(name);
        return false;
    }
    const std::string previous = parameters::text(*entry);
    if (!parameters::settext(*entry, value, error))
    {
        if (error != nullptr && error->empty())
            *error = "invalid value " + std::string(value) + " for " + std::string(name);
        return false;
    }
    parameters::settext(*entry, previous);
    return true;
}

bool loadrecipeparameters(const std::filesystem::path& path, std::string* error)
{
    std::ifstream file(path, std::ios::binary);
    if (!file)
    {
        if (error != nullptr)
            *error = "cannot read " + path.string();
        return false;
    }
    std::ostringstream content;
    content << file.rdbuf();
    const std::string text = content.str();
    if (error != nullptr)
        error->clear();
    ValueReader reader{ text, 0, error };
    return reader.readobject("");
}

bool validaterecipeparameters(std::string* error)
{
    namespace atmosphere = tuning::climate::atmosphere;
    const char* problem = nullptr;
    if (atmosphere::katabaticShallowTenMetreDiagnostic && !atmosphere::separateKatabaticTransport)
        problem = "climate.atmosphere.katabatic_shallow_ten_metre_diagnostic requires separate_katabatic_transport";
    else if (atmosphere::katabaticShallowTenMetreDiagnostic &&
        !(atmosphere::katabaticShallowDepthMetres > 20.0f &&
            atmosphere::katabaticShallowDepthMetres <= atmosphere::surfaceBoundaryLayerMomentumDepthMetres))
        problem = "climate.atmosphere.katabatic_shallow_depth_metres must contain 10 m and fit within surface_boundary_layer_momentum_depth_metres";
    else if (!(tuning::world::ocean::landFraction >= 0.0f && tuning::world::ocean::landFraction < 1.0f))
        problem = "world.ocean.land_fraction must be at least 0 (0 keeps the tectonics' sea level) and below 1";
    if (problem != nullptr && error != nullptr)
        *error = problem;
    return problem == nullptr;
}

std::string recipeparametersjson()
{
    std::string json = "{";
    bool first = true;
    for (const parameters::Entry& entry : recipeparameters())
    {
        json += first ? "\n" : ",\n";
        json += "    \"" + std::string(entry.name) + "\": " + parameters::text(entry);
        first = false;
    }
    json += first ? "}" : "\n  }";
    return json;
}
