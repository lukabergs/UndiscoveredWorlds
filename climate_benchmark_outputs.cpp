#include "climate_benchmark_outputs.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <sstream>

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
    case climatebenchmarkmapkind::surfacedivergence: return "s-div";
    case climatebenchmarkmapkind::moisturefluxconvergence: return "moist-conv";
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
    constexpr std::array<climatebenchmarkmapkind, 13> kinds = {
        climatebenchmarkmapkind::surfacewindlic,
        climatebenchmarkmapkind::surfacewindparticles,
        climatebenchmarkmapkind::upperwindlic,
        climatebenchmarkmapkind::upperwindparticles,
        climatebenchmarkmapkind::surfacewindvectorerror,
        climatebenchmarkmapkind::upperwindvectorerror,
        climatebenchmarkmapkind::surfacewindspeed,
        climatebenchmarkmapkind::surfacedivergence,
        climatebenchmarkmapkind::moisturefluxconvergence,
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
    selection.add(climatebenchmarkmapkind::surfacewindlic, 0);
    selection.add(climatebenchmarkmapkind::surfacewindparticles, 0);
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
