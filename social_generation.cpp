#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_map>

#include "functions.hpp"
#include "planet.hpp"
#include "social_generation_internal.hpp"

using namespace std;

namespace
{
string trim(string value)
{
    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front())))
        value.erase(value.begin());

    while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back())))
        value.pop_back();

    return value;
}

vector<string> parsecsvline(const string& line)
{
    vector<string> fields;
    string current;
    bool inquotes = false;

    for (size_t i = 0; i < line.size(); i++)
    {
        const char ch = line[i];

        if (ch == '"')
        {
            if (inquotes && i + 1 < line.size() && line[i + 1] == '"')
            {
                current.push_back('"');
                i++;
            }
            else
            {
                inquotes = !inquotes;
            }
        }
        else if (ch == ',' && !inquotes)
        {
            fields.push_back(trim(current));
            current.clear();
        }
        else
        {
            current.push_back(ch);
        }
    }

    fields.push_back(trim(current));
    return fields;
}

float stofsafe(const string& value, float fallback)
{
    if (value.empty())
        return fallback;

    try
    {
        return stof(value);
    }
    catch (const exception&)
    {
        return fallback;
    }
}

int stoisafe(const string& value, int fallback)
{
    if (value.empty())
        return fallback;

    try
    {
        return stoi(value);
    }
    catch (const exception&)
    {
        return fallback;
    }
}

vector<vector<string>> loadcsvrows(const string& filepath)
{
    vector<vector<string>> rows;
    ifstream infile(filepath);

    if (!infile.is_open())
        return rows;

    string line;

    while (getline(infile, line))
    {
        string trimmed = trim(line);

        if (trimmed.empty() || trimmed[0] == '#')
            continue;

        rows.push_back(parsecsvline(trimmed));
    }

    return rows;
}

unordered_map<string, int> buildheaderindex(const vector<string>& headers)
{
    unordered_map<string, int> index;

    for (int i = 0; i < static_cast<int>(headers.size()); i++)
        index.emplace(headers[i], i);

    return index;
}

string getcsvvalue(const vector<string>& row, const unordered_map<string, int>& index, const string& key)
{
    const auto found = index.find(key);

    if (found == index.end())
        return "";

    const int column = found->second;

    if (column < 0 || column >= static_cast<int>(row.size()))
        return "";

    return row[column];
}
}

namespace socialgen
{
bool loadsocialprofilecsv(const string& filepath, SocialProfile& profile)
{
    const vector<vector<string>> rows = loadcsvrows(filepath);

    if (rows.size() < 2)
        return false;

    const unordered_map<string, int> headers = buildheaderindex(rows[0]);
    const vector<string>& row = rows[1];

    profile.profileKey = getcsvvalue(row, headers, "profile_key");
    profile.coastWeight = stofsafe(getcsvvalue(row, headers, "coast_weight"), profile.coastWeight);
    profile.harborWeight = stofsafe(getcsvvalue(row, headers, "harbor_weight"), profile.harborWeight);
    profile.riverWeight = stofsafe(getcsvvalue(row, headers, "river_weight"), profile.riverWeight);
    profile.confluenceWeight = stofsafe(getcsvvalue(row, headers, "confluence_weight"), profile.confluenceWeight);
    profile.climateWeight = stofsafe(getcsvvalue(row, headers, "climate_weight"), profile.climateWeight);
    profile.agricultureWeight = stofsafe(getcsvvalue(row, headers, "agriculture_weight"), profile.agricultureWeight);
    profile.terrainPenalty = stofsafe(getcsvvalue(row, headers, "terrain_penalty"), profile.terrainPenalty);
    profile.routeWeight = stofsafe(getcsvvalue(row, headers, "route_weight"), profile.routeWeight);
    profile.resourceWeight = stofsafe(getcsvvalue(row, headers, "resource_weight"), profile.resourceWeight);
    profile.spacingRadius = stoisafe(getcsvvalue(row, headers, "spacing_radius"), profile.spacingRadius);
    profile.urbanizationBias = stofsafe(getcsvvalue(row, headers, "urbanization_bias"), profile.urbanizationBias);

    profile.spacingRadius = std::max(6, profile.spacingRadius);
    profile.urbanizationBias = std::clamp(profile.urbanizationBias, 0.05f, 0.75f);
    return true;
}

void loadsocialtables(const string& rootpath, SocialTables& tables)
{
    tables = {};

    const string commoditiespath = (filesystem::path(rootpath) / "social_commodities.csv").string();
    const string recipespath = (filesystem::path(rootpath) / "social_recipes.csv").string();
    const string knowledgespath = (filesystem::path(rootpath) / "social_knowledges.csv").string();

    const auto commodityrows = loadcsvrows(commoditiespath);
    if (commodityrows.size() >= 2)
    {
        const unordered_map<string, int> headers = buildheaderindex(commodityrows[0]);

        for (size_t i = 1; i < commodityrows.size(); i++)
        {
            CommodityDefinition definition;
            definition.key = getcsvvalue(commodityrows[i], headers, "commodity_key");
            definition.stage = getcsvvalue(commodityrows[i], headers, "stage");
            definition.baseValue = stofsafe(getcsvvalue(commodityrows[i], headers, "base_value"), 1.0f);
            definition.weight = stofsafe(getcsvvalue(commodityrows[i], headers, "weight"), 1.0f);
            definition.sourceType = getcsvvalue(commodityrows[i], headers, "source_type");
            definition.requiredKnowledge = getcsvvalue(commodityrows[i], headers, "required_knowledge");

            if (!definition.key.empty())
                tables.commodities.push_back(std::move(definition));
        }
    }

    const auto reciperows = loadcsvrows(recipespath);
    if (reciperows.size() >= 2)
    {
        const unordered_map<string, int> headers = buildheaderindex(reciperows[0]);

        for (size_t i = 1; i < reciperows.size(); i++)
        {
            RecipeDefinition recipe;
            recipe.outputKey = getcsvvalue(reciperows[i], headers, "output_key");
            recipe.inputKey = getcsvvalue(reciperows[i], headers, "input_key");
            recipe.inputAmount = stofsafe(getcsvvalue(reciperows[i], headers, "input_amount"), 1.0f);
            recipe.laborType = getcsvvalue(reciperows[i], headers, "labor_type");
            recipe.knowledgeGate = getcsvvalue(reciperows[i], headers, "knowledge_gate");

            if (!recipe.outputKey.empty())
                tables.recipes.push_back(std::move(recipe));
        }
    }

    const auto knowledgerows = loadcsvrows(knowledgespath);
    if (knowledgerows.size() >= 2)
    {
        const unordered_map<string, int> headers = buildheaderindex(knowledgerows[0]);

        for (size_t i = 1; i < knowledgerows.size(); i++)
        {
            KnowledgeDefinition knowledge;
            knowledge.key = getcsvvalue(knowledgerows[i], headers, "knowledge_key");
            knowledge.spreadFactor = stofsafe(getcsvvalue(knowledgerows[i], headers, "spread_factor"), 0.0f);
            knowledge.productivityBonus = stofsafe(getcsvvalue(knowledgerows[i], headers, "productivity_bonus"), 0.0f);
            knowledge.cohesionBonus = stofsafe(getcsvvalue(knowledgerows[i], headers, "cohesion_bonus"), 0.0f);
            knowledge.industryUnlocks = getcsvvalue(knowledgerows[i], headers, "industry_unlocks");

            if (!knowledge.key.empty())
                tables.knowledges.push_back(std::move(knowledge));
        }
    }
}
}

void generatesocialworld(planet& world, const SocialGenerationOptions& options)
{
    if (!options.enabled)
        return;

    updatereport("Generating social world");

    socialgen::SocialProfile profile;
    socialgen::SocialTables tables;
    const string extrapath = (filesystem::current_path() / "extra").string();

    if (!socialgen::loadsocialprofilecsv((filesystem::path(extrapath) / "social_profiles.csv").string(), profile))
        updatereport("social_profiles.csv missing or invalid, using defaults");

    socialgen::loadsocialtables(extrapath, tables);

    std::mt19937_64 rng(deterministiccontextseed(world.seed(), 0x53534f43));

    updatereport("Social step: suitability and settlement siting");
    socialgen::computesitesandinfrastructure(world, profile, rng);

    updatereport("Social step: routes");
    socialgen::buildroutes(world, profile, rng);

    updatereport("Social step: polities");
    socialgen::buildpolities(world, options, rng);

    updatereport("Social step: trade and knowledge");
    socialgen::runtradeandknowledge(world, tables, rng);

    if (options.mode == SocialGenerationOptions::Mode::historical)
    {
        updatereport("Social step: historical simulation");
        socialgen::runhistoricalmode(world, options, rng);
    }
}
