#include "core/world/world_fields.hpp"

#include <array>
#include <stdexcept>

namespace worldmodel
{
namespace
{
const std::array<FieldSpec, kFieldCount> kFields =
{ {
#define UW_WORLD_FIELD_SPEC(id, component, element, layout, units, description) \
    { FieldId::id, Component::component, ElementType::element, Layout::layout, #id, units, description },
    UW_WORLD_FIELD_LIST(UW_WORLD_FIELD_SPEC)
#undef UW_WORLD_FIELD_SPEC
} };

const std::array<TableSpec, kTableCount> kTables =
{ {
    { TableId::tectonic_boundary_segments, Component::geosphere, "tectonic_boundary_segments", "Plate boundary segments" },
    { TableId::tectonic_deforming_regions, Component::geosphere, "tectonic_deforming_regions", "Deforming regions" },
    { TableId::mineral_deposits, Component::geosphere, "mineral_deposits", "Mineral deposits formed through tectonic history" },
    { TableId::mineral_budget, Component::geosphere, "mineral_budget", "Gold budget of the mineral history" },
    { TableId::craters, Component::geosphere, "craters", "Crater centres in drawing order" },
    { TableId::settlements, Component::anthroposphere, "settlements", "Settlements" },
    { TableId::polities, Component::anthroposphere, "polities", "Polities" },
    { TableId::route_edges, Component::anthroposphere, "route_edges", "Route network edges" },
    { TableId::history_events, Component::anthroposphere, "history_events", "Historical events" },
    { TableId::weather_state, Component::atmosphere, "weather_state", "Serialized transient weather anomaly state" },
    { TableId::atmosphere_levels, Component::atmosphere, "levels", "Atmosphere column levels: name and nominal pressure" },
    { TableId::tectonic_snapshot, Component::geosphere, "tectonic_snapshot", "Record of the last tectonic update: time, scenario, engine diagnostics, statistics and plates" },
} };

constexpr std::array<const char*, 8> kComponents =
{ "planetbody", "geosphere", "hydrosphere", "cryosphere", "atmosphere", "biosphere", "anthroposphere", "presentation" };
}

const FieldSpec& fieldspec(FieldId id)
{
    const auto index = static_cast<std::size_t>(id);
    if (index >= kFields.size())
        throw std::out_of_range("Unknown world field");
    return kFields[index];
}

const TableSpec& tablespec(TableId id)
{
    const auto index = static_cast<std::size_t>(id);
    if (index >= kTables.size() || kTables[index].id != id)
        throw std::out_of_range("Unknown world table");
    return kTables[index];
}

const char* componentname(Component component)
{
    return kComponents.at(static_cast<std::size_t>(component));
}

const char* elementname(ElementType element)
{
    switch (element)
    {
    case ElementType::b8: return "|b1";
    case ElementType::u8: return "|u1";
    case ElementType::i16: return "<i2";
    case ElementType::u16: return "<u2";
    case ElementType::i32: return "<i4";
    case ElementType::f32: return "<f4";
    }
    return "";
}

std::size_t elementsize(ElementType element)
{
    switch (element)
    {
    case ElementType::b8:
    case ElementType::u8: return 1;
    case ElementType::i16:
    case ElementType::u16: return 2;
    case ElementType::i32:
    case ElementType::f32: return 4;
    }
    return 0;
}

const char* layoutname(Layout layout)
{
    switch (layout)
    {
    case Layout::scalar: return "scalar";
    case Layout::grid: return "grid";
    case Layout::seasonal: return "seasonal";
    case Layout::monthly: return "monthly";
    case Layout::zonalprofile: return "zonalprofile";
    case Layout::levelseasonal: return "levelseasonal";
    }
    return "";
}

std::string qualifiedname(FieldId id)
{
    const FieldSpec& spec = fieldspec(id);
    return std::string(componentname(spec.component)) + "." + spec.name;
}

std::string qualifiedname(TableId id)
{
    const TableSpec& spec = tablespec(id);
    return std::string(componentname(spec.component)) + "." + spec.name;
}

bool findfield(std::string_view qualified, FieldId& id)
{
    for (const FieldSpec& spec : kFields)
    {
        if (qualifiedname(spec.id) == qualified)
        {
            id = spec.id;
            return true;
        }
    }
    return false;
}

bool findtable(std::string_view qualified, TableId& id)
{
    for (const TableSpec& spec : kTables)
    {
        if (qualifiedname(spec.id) == qualified)
        {
            id = spec.id;
            return true;
        }
    }
    return false;
}

bool findcomponent(std::string_view name, Component& component)
{
    for (std::size_t index = 0; index < kComponents.size(); index++)
    {
        if (name == kComponents[index])
        {
            component = static_cast<Component>(index);
            return true;
        }
    }
    return false;
}
}
