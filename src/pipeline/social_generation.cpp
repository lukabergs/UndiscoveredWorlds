#include <filesystem>

#include "deterministic_random.hpp"
#include "generation_progress.hpp"
#include "planet.hpp"
#include "social_definitions.hpp"

using namespace std;

void generatesocialworld(planet& world, const SocialGenerationOptions& options)
{
    if (!options.enabled)
        return;

    updatereport("Generating social world");

    socialgen::SocialProfile profile;
    socialgen::SocialTables tables;
    const string extrapath = (filesystem::current_path() / "definitions" / "society").string();

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
