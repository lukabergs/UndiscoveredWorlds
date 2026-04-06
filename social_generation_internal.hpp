#pragma once

#include <random>
#include <string>
#include <vector>

#include "social_generation.hpp"

class planet;

namespace socialgen
{
struct SocialProfile
{
    std::string profileKey = "human_default";
    float coastWeight = 1.0f;
    float harborWeight = 1.0f;
    float riverWeight = 1.0f;
    float confluenceWeight = 0.7f;
    float climateWeight = 1.0f;
    float agricultureWeight = 1.2f;
    float terrainPenalty = 1.0f;
    float routeWeight = 0.8f;
    float resourceWeight = 0.8f;
    int spacingRadius = 26;
    float urbanizationBias = 0.25f;
};

struct CommodityDefinition
{
    std::string key;
    std::string stage;
    float baseValue = 0.0f;
    float weight = 0.0f;
    std::string sourceType;
    std::string requiredKnowledge;
};

struct RecipeDefinition
{
    std::string outputKey;
    std::string inputKey;
    float inputAmount = 0.0f;
    std::string laborType;
    std::string knowledgeGate;
};

struct KnowledgeDefinition
{
    std::string key;
    float spreadFactor = 0.0f;
    float productivityBonus = 0.0f;
    float cohesionBonus = 0.0f;
    std::string industryUnlocks;
};

struct SocialTables
{
    std::vector<CommodityDefinition> commodities;
    std::vector<RecipeDefinition> recipes;
    std::vector<KnowledgeDefinition> knowledges;
};

bool loadsocialprofilecsv(const std::string& filepath, SocialProfile& profile);
void loadsocialtables(const std::string& rootpath, SocialTables& tables);
void computesitesandinfrastructure(planet& world, const SocialProfile& profile, std::mt19937_64& rng);
void buildroutes(planet& world, const SocialProfile& profile, std::mt19937_64& rng);
void buildpolities(planet& world, const SocialGenerationOptions& options, std::mt19937_64& rng);
void runtradeandknowledge(planet& world, const SocialTables& tables, std::mt19937_64& rng);
void runhistoricalmode(planet& world, const SocialGenerationOptions& options, std::mt19937_64& rng);
}
