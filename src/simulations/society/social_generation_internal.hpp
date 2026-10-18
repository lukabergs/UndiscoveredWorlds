#pragma once

#include <random>
#include <string>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "simulations/society/social_generation.hpp"

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

// Neighbour slots in the order the former (x, y) loops visited them (dx outer,
// dy inner, y growing south): ties in the queues and sorts keep that order.
inline constexpr int scanslots[8] = { 8, 7, 6, 1, 5, 2, 3, 4 };
inline bool orthogonalslot(int slot) { return (slot & 1) == 1; }

// The lat-lon export's columns W of the world's grid (columns on lat-lon, 4 N on the cube).
int exportcolumns(const mesh::CellGrid& grid);
// The distance between two cell centres in grid steps of the lat-lon grid
// (wrapped column offset and row offset, Euclidean).
float celldistance(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b);
// The same distance as the sum of the two offsets (the sea routes' landing test).
int cellmanhattan(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b);
// The squared Euclidean offset in grid steps (the capital spacing test).
int celldistancesquared(const mesh::CellGrid& grid, mesh::CellId a, mesh::CellId b);
// Places a settlement on a cell: its cell id, and x, y on the lat-lon export
// grid (the cell's column and row on lat-lon; on the cube the export cell of its centre).
void placesettlement(Settlement& settlement, const mesh::CellGrid& grid, mesh::CellId cell);

void computesitesandinfrastructure(planet& world, const SocialProfile& profile, std::mt19937_64& rng);
void buildroutes(planet& world, const SocialProfile& profile, std::mt19937_64& rng);
void buildpolities(planet& world, const SocialGenerationOptions& options, std::mt19937_64& rng);
void runtradeandknowledge(planet& world, const SocialTables& tables, std::mt19937_64& rng);
void runhistoricalmode(planet& world, const SocialGenerationOptions& options, std::mt19937_64& rng);
}
