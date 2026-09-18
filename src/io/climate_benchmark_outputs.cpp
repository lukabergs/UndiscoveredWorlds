#include "climate_benchmark_outputs.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace
{
constexpr int seasoncount = 4;
constexpr std::array<const char*, seasoncount> seasonnames = {
    "jan", "apr", "jul", "oct"
};

bool isseasonal(climatebenchmarkmapkind kind)
{
    return kind != climatebenchmarkmapkind::koppen &&
        kind != climatebenchmarkmapkind::temperature &&
        kind != climatebenchmarkmapkind::precipitation &&
        kind != climatebenchmarkmapkind::precipitationtiff;
}

const char* seasonalmapname(climatebenchmarkmapkind kind)
{
    switch (kind)
    {
    case climatebenchmarkmapkind::surfacewindlic: return "s-wind-lic";
    case climatebenchmarkmapkind::surfacewindparticles: return "s-wind-part";
    case climatebenchmarkmapkind::upperwindlic: return "u-wind-lic";
    case climatebenchmarkmapkind::upperwindparticles: return "u-wind-part";
    case climatebenchmarkmapkind::surfacewindvectorerror: return "s-wind-err";
    case climatebenchmarkmapkind::upperwindvectorerror: return "u-wind-err";
    case climatebenchmarkmapkind::surfacewindspeed: return "s-wind-speed";
    case climatebenchmarkmapkind::upperwindspeed: return "u-wind-speed";
    case climatebenchmarkmapkind::surfacedivergence: return "s-div";
    case climatebenchmarkmapkind::moisturefluxconvergence: return "moist-conv";
    case climatebenchmarkmapkind::boundarymoistureflux: return "b-moist-flux";
    case climatebenchmarkmapkind::freemoistureflux: return "f-moist-flux";
    case climatebenchmarkmapkind::columnmoistureflux: return "moist-flux";
    case climatebenchmarkmapkind::verticalascent: return "ascent";
    case climatebenchmarkmapkind::surfacewindconsistency: return "s-wind-cons";
    case climatebenchmarkmapkind::upperwindconsistency: return "u-wind-cons";
    case climatebenchmarkmapkind::seasurfacetemperature: return "sst";
    case climatebenchmarkmapkind::oceancurrents: return "ocean-current";
    case climatebenchmarkmapkind::columnheating: return "heating";
    case climatebenchmarkmapkind::era5surfacewindlic: return "era5-s-wind-lic";
    case climatebenchmarkmapkind::era5surfacewindparticles: return "era5-s-wind-part";
    case climatebenchmarkmapkind::era5upperwindlic: return "era5-u-wind-lic";
    case climatebenchmarkmapkind::era5upperwindparticles: return "era5-u-wind-part";
    default: return "";
    }
}

bool parseseason(const std::string& value, int& season)
{
    for (int index = 0; index < seasoncount; index++)
    {
        if (value == seasonnames[index])
        {
            season = index;
            return true;
        }
    }

    return false;
}

bool parsekind(const std::string& value, climatebenchmarkmapkind& kind)
{
    constexpr std::array<climatebenchmarkmapkind, 23> kinds = {
        climatebenchmarkmapkind::surfacewindlic,
        climatebenchmarkmapkind::surfacewindparticles,
        climatebenchmarkmapkind::upperwindlic,
        climatebenchmarkmapkind::upperwindparticles,
        climatebenchmarkmapkind::surfacewindvectorerror,
        climatebenchmarkmapkind::upperwindvectorerror,
        climatebenchmarkmapkind::surfacewindspeed,
        climatebenchmarkmapkind::upperwindspeed,
        climatebenchmarkmapkind::surfacedivergence,
        climatebenchmarkmapkind::moisturefluxconvergence,
        climatebenchmarkmapkind::boundarymoistureflux,
        climatebenchmarkmapkind::freemoistureflux,
        climatebenchmarkmapkind::columnmoistureflux,
        climatebenchmarkmapkind::verticalascent,
        climatebenchmarkmapkind::surfacewindconsistency,
        climatebenchmarkmapkind::upperwindconsistency,
        climatebenchmarkmapkind::seasurfacetemperature,
        climatebenchmarkmapkind::oceancurrents,
        climatebenchmarkmapkind::columnheating,
        climatebenchmarkmapkind::era5surfacewindlic,
        climatebenchmarkmapkind::era5surfacewindparticles,
        climatebenchmarkmapkind::era5upperwindlic,
        climatebenchmarkmapkind::era5upperwindparticles
    };

    for (const climatebenchmarkmapkind candidate : kinds)
    {
        if (value == seasonalmapname(candidate))
        {
            kind = candidate;
            return true;
        }
    }

    return false;
}

std::string trim(std::string value)
{
    const auto notspace = [](unsigned char character) { return !std::isspace(character); };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), notspace));
    value.erase(std::find_if(value.rbegin(), value.rend(), notspace).base(), value.end());
    return value;
}
}

climatebenchmarkmapselection climatebenchmarkmapselection::defaultselection()
{
    climatebenchmarkmapselection selection;
    selection.add(climatebenchmarkmapkind::koppen);
    selection.add(climatebenchmarkmapkind::temperature);
    selection.add(climatebenchmarkmapkind::surfacewindlic, 0);
    selection.add(climatebenchmarkmapkind::surfacewindparticles, 0);
    selection.add(climatebenchmarkmapkind::surfacewindlic, 2);
    selection.add(climatebenchmarkmapkind::surfacewindparticles, 2);
    for (const int season : { 0, 2 })
    {
        selection.add(climatebenchmarkmapkind::surfacewindspeed, season);
        selection.add(climatebenchmarkmapkind::upperwindspeed, season);
        selection.add(climatebenchmarkmapkind::upperwindlic, season);
        selection.add(climatebenchmarkmapkind::upperwindparticles, season);
    }
    selection.add(climatebenchmarkmapkind::seasurfacetemperature, 0);
    selection.add(climatebenchmarkmapkind::seasurfacetemperature, 2);
    selection.add(climatebenchmarkmapkind::oceancurrents, 0);
    selection.add(climatebenchmarkmapkind::oceancurrents, 2);
    selection.add(climatebenchmarkmapkind::precipitation);
    selection.add(climatebenchmarkmapkind::precipitationtiff);
    return selection;
}

climatebenchmarkmapselection climatebenchmarkmapselection::allselection()
{
    climatebenchmarkmapselection selection;
    selection.add(climatebenchmarkmapkind::koppen);
    selection.add(climatebenchmarkmapkind::temperature);
    selection.add(climatebenchmarkmapkind::precipitation);
    selection.add(climatebenchmarkmapkind::precipitationtiff);

    for (int season = 0; season < seasoncount; season++)
    {
        for (int kind = static_cast<int>(climatebenchmarkmapkind::surfacewindlic);
             kind <= static_cast<int>(climatebenchmarkmapkind::era5upperwindparticles);
             kind++)
        {
            selection.add(static_cast<climatebenchmarkmapkind>(kind), season);
        }
    }

    return selection;
}

void climatebenchmarkmapselection::setwindmapcolumns(int columns)
{
    if (columns < 64 || columns > 2048 || columns % 2 != 0)
        throw std::invalid_argument("Wind map width must be even and between 64 and 2048");
    windmapcolumns_ = columns;
}

void climatebenchmarkmapselection::clear()
{
    requests_.clear();
}

void climatebenchmarkmapselection::add(climatebenchmarkmapkind kind, int season)
{
    if (!isseasonal(kind))
        season = -1;

    if (season < -1 || season >= seasoncount || includes(kind, season))
        return;

    requests_.push_back({ kind, season });
}

bool climatebenchmarkmapselection::includes(climatebenchmarkmapkind kind, int season) const
{
    return std::any_of(requests_.begin(), requests_.end(), [&](const auto& request)
    {
        return request.kind == kind && request.season == season;
    });
}

bool climatebenchmarkmapselection::requirescirculation() const
{
    return std::any_of(requests_.begin(), requests_.end(), [](const auto& request)
    {
        return isseasonal(request.kind) && !climatebenchmarkmapisreference(request.kind);
    });
}

bool climatebenchmarkmapselection::requiresreferencewinds() const
{
    return std::any_of(requests_.begin(), requests_.end(), [](const auto& request)
    {
        return climatebenchmarkmapisreference(request.kind) ||
            request.kind == climatebenchmarkmapkind::surfacewindvectorerror ||
            request.kind == climatebenchmarkmapkind::upperwindvectorerror;
    });
}

const std::vector<climatebenchmarkmaprequest>& climatebenchmarkmapselection::requests() const
{
    return requests_;
}

bool addclimatebenchmarkmapargument(
    climatebenchmarkmapselection& selection,
    const std::string& argument,
    std::string* failuremessage)
{
    std::stringstream stream(argument);
    std::string token;

    while (std::getline(stream, token, ','))
    {
        token = trim(token);
        if (token.empty())
            continue;

        if (token == "none")
        {
            selection.clear();
            continue;
        }
        if (token == "all")
        {
            const climatebenchmarkmapselection all = climatebenchmarkmapselection::allselection();
            for (const auto& request : all.requests())
                selection.add(request.kind, request.season);
            continue;
        }
        if (token == "koppen")
        {
            selection.add(climatebenchmarkmapkind::koppen);
            continue;
        }
        if (token == "temp")
        {
            selection.add(climatebenchmarkmapkind::temperature);
            continue;
        }
        if (token == "precip")
        {
            selection.add(climatebenchmarkmapkind::precipitation);
            continue;
        }
        if (token == "precip-tif")
        {
            selection.add(climatebenchmarkmapkind::precipitationtiff);
            continue;
        }

        const std::size_t separator = token.find('-');
        int season = -1;
        climatebenchmarkmapkind kind = climatebenchmarkmapkind::koppen;
        if (separator == std::string::npos ||
            !parseseason(token.substr(0, separator), season) ||
            !parsekind(token.substr(separator + 1), kind))
        {
            if (failuremessage != nullptr)
                *failuremessage = "Unknown benchmark map id: " + token;
            return false;
        }

        selection.add(kind, season);
    }

    return true;
}

const char* climatebenchmarkseasonname(int season)
{
    return season >= 0 && season < seasoncount ? seasonnames[season] : "";
}

std::string climatebenchmarkmapid(const climatebenchmarkmaprequest& request)
{
    switch (request.kind)
    {
    case climatebenchmarkmapkind::koppen: return "koppen";
    case climatebenchmarkmapkind::temperature: return "temp";
    case climatebenchmarkmapkind::precipitation: return "precip";
    case climatebenchmarkmapkind::precipitationtiff: return "precip-tif";
    default:
        return std::string(climatebenchmarkseasonname(request.season)) + "-" +
            seasonalmapname(request.kind);
    }
}

std::string climatebenchmarkmapfilename(
    const climatebenchmarkmaprequest& request,
    int runid)
{
    if (request.kind == climatebenchmarkmapkind::precipitationtiff)
        return std::to_string(runid) + "_precip.tif";

    std::string id = climatebenchmarkmapid(request);
    std::replace(id.begin(), id.end(), '-', '_');
    return std::to_string(runid) + "_" + id + ".png";
}

std::string climatebenchmarkmaprelativepath(const climatebenchmarkmaprequest& request, int runid)
{
    const std::string season = std::to_string(request.season + 1);
    std::string directory;
    switch (request.kind)
    {
    case climatebenchmarkmapkind::koppen: directory = "koppen"; break;
    case climatebenchmarkmapkind::temperature: directory = "air_temp/annual/s"; break;
    case climatebenchmarkmapkind::precipitation:
    case climatebenchmarkmapkind::precipitationtiff: directory = "rain/annual/s"; break;
    case climatebenchmarkmapkind::surfacewindlic: directory = "wind/lic/" + season + "/s"; break;
    case climatebenchmarkmapkind::upperwindlic: directory = "wind/lic/" + season + "/u"; break;
    case climatebenchmarkmapkind::surfacewindparticles: directory = "wind/particles/" + season + "/s"; break;
    case climatebenchmarkmapkind::upperwindparticles: directory = "wind/particles/" + season + "/u"; break;
    case climatebenchmarkmapkind::surfacewindvectorerror: directory = "wind/error/" + season + "/s"; break;
    case climatebenchmarkmapkind::upperwindvectorerror: directory = "wind/error/" + season + "/u"; break;
    case climatebenchmarkmapkind::surfacewindspeed: directory = "wind/speed/" + season + "/s"; break;
    case climatebenchmarkmapkind::upperwindspeed: directory = "wind/speed/" + season + "/u"; break;
    case climatebenchmarkmapkind::surfacedivergence: directory = "wind/divergence/" + season + "/s"; break;
    case climatebenchmarkmapkind::moisturefluxconvergence: directory = "moisture/convergence/" + season + "/column"; break;
    case climatebenchmarkmapkind::boundarymoistureflux: directory = "moisture/flux/" + season + "/b"; break;
    case climatebenchmarkmapkind::freemoistureflux: directory = "moisture/flux/" + season + "/f"; break;
    case climatebenchmarkmapkind::columnmoistureflux: directory = "moisture/flux/" + season + "/column"; break;
    case climatebenchmarkmapkind::verticalascent: directory = "wind/ascent/" + season + "/u"; break;
    case climatebenchmarkmapkind::surfacewindconsistency: directory = "wind/consistency/" + season + "/s"; break;
    case climatebenchmarkmapkind::upperwindconsistency: directory = "wind/consistency/" + season + "/u"; break;
    case climatebenchmarkmapkind::seasurfacetemperature: directory = "sea_temp/" + season + "/s"; break;
    case climatebenchmarkmapkind::oceancurrents: directory = "ocean/current/" + season + "/s"; break;
    case climatebenchmarkmapkind::columnheating: directory = "energy/heating/" + season + "/column"; break;
    default: return {}; // Reference previews have their own resolution-based cache.
    }
    return "climate/" + directory + "/" + std::to_string(runid) +
        (request.kind == climatebenchmarkmapkind::precipitationtiff ? ".tif" : ".png");
}

bool climatebenchmarkmapisreference(climatebenchmarkmapkind kind)
{
    return kind == climatebenchmarkmapkind::era5surfacewindlic ||
        kind == climatebenchmarkmapkind::era5surfacewindparticles ||
        kind == climatebenchmarkmapkind::era5upperwindlic ||
        kind == climatebenchmarkmapkind::era5upperwindparticles;
}

double climatebenchmarkweightedrelativeerror(
    const std::vector<long long>& simulationcounts,
    const std::vector<double>& referencecounts,
    int horizontalresolution,
    int referencehorizontalresolution)
{
    if (horizontalresolution <= 0 || referencehorizontalresolution <= 0 ||
        simulationcounts.size() != referencecounts.size())
    {
        return 0.0;
    }

    const double linearscale = static_cast<double>(referencehorizontalresolution) /
        static_cast<double>(horizontalresolution);
    const double areascale = linearscale * linearscale;
    double absoluteerror = 0.0;
    double simulationtotal = 0.0;

    for (std::size_t index = 0; index < simulationcounts.size(); index++)
    {
        const double scaledreference =
            referencecounts[index] / areascale;
        absoluteerror += std::abs(
            static_cast<double>(simulationcounts[index]) - scaledreference);
        simulationtotal += static_cast<double>(simulationcounts[index]);
    }

    return simulationtotal > 0.0 ? absoluteerror / simulationtotal : 0.0;
}
