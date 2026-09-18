#pragma once

#include <cstdint>
#include <string>

class planet;

struct SocialGenerationOptions
{
    bool enabled = false;

    enum class Mode
    {
        static_ex_nihilo,
        historical
    };

    Mode mode = Mode::static_ex_nihilo;
    bool usePrehistory = true;
    int historyYears = 1200;
};

enum class RouteMode : std::uint8_t
{
    land = 0,
    river = 1,
    sea = 2
};

struct Settlement
{
    int id = -1;
    int x = 0;
    int y = 0;
    std::string name;
    int urbanPopulation = 0;
    int ruralPopulation = 0;
    float carryingCapacity = 0.0f;
    float infrastructure = 0.0f;
    float marketStrength = 0.0f;
    float harbor = 0.0f;
    float riverAccess = 0.0f;
    int polityId = -1;
};

struct Polity
{
    int id = -1;
    std::string name;
    int capitalSettlementId = -1;
    float cohesion = 0.0f;
    int population = 0;
    float infrastructure = 0.0f;
    float wealth = 0.0f;
    float militaryPressure = 0.0f;
};

struct RouteEdge
{
    int fromSettlementId = -1;
    int toSettlementId = -1;
    RouteMode mode = RouteMode::land;
    float cost = 0.0f;
    float capacity = 0.0f;
    float traffic = 0.0f;
};

struct HistoryEvent
{
    int year = 0;
    std::string type;
    int primarySettlementId = -1;
    int primaryPolityId = -1;
    int secondarySettlementId = -1;
    int secondaryPolityId = -1;
    std::string summary;
};

void generatesocialworld(planet& world, const SocialGenerationOptions& options);
