#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

#include "climate_atmosphere.hpp"
#include "climate_circulation_diagnostics.hpp"
#include "climate_energy.hpp"
#include "climate_grid.hpp"
#include "climate_hydrology.hpp"
#include "climate_ocean_dynamics.hpp"
#include "climate_physics.hpp"
#include "climate_weather.hpp"
#include "generation_tuning.hpp"
#include "planet.hpp"
#include "functions.hpp"

using namespace std;

namespace
{
using floatgrid = vector<vector<float>>;

template<typename T>
class rowmajorgrid
{
public:
    class columnview
    {
    public:
        columnview(T* values, int columns, int x)
            : values_(values), columns_(columns), x_(x) {}

        T& operator[](int y) const
        {
            return values_[static_cast<size_t>(y) * columns_ + x_];
        }

    private:
        T* values_;
        int columns_;
        int x_;
    };

    class constcolumnview
    {
    public:
        constcolumnview(const T* values, int columns, int x)
            : values_(values), columns_(columns), x_(x) {}

        const T& operator[](int y) const
        {
            return values_[static_cast<size_t>(y) * columns_ + x_];
        }

    private:
        const T* values_;
        int columns_;
        int x_;
    };

    rowmajorgrid() = default;
    rowmajorgrid(int columns, int rows, const T& value = T{})
        : columns_(columns), rows_(rows), values_(
            static_cast<size_t>(columns) * static_cast<size_t>(rows), value) {}

    columnview operator[](int x)
    {
        return columnview(values_.data(), columns_, x);
    }

    constcolumnview operator[](int x) const
    {
        return constcolumnview(values_.data(), columns_, x);
    }

    void fill(const T& value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }

private:
    int columns_ = 0;
    int rows_ = 0;
    vector<T> values_;
};

using hydrologyfloatgrid = rowmajorgrid<float>;
using hydrologybytegrid = rowmajorgrid<signed char>;

constexpr std::array<float, CLIMATESEASONCOUNT> seasonlatitudephase = { -1.0f, 0.0f, 1.0f, 0.0f };

struct circulationfloatcache
{
    planet* source = nullptr;
    int width = -1;
    int height = -1;
    std::array<floatgrid, CLIMATESEASONCOUNT> pressure;
    std::array<floatgrid, CLIMATESEASONCOUNT> upperheight;
    std::array<bool, CLIMATESEASONCOUNT> populated{};

    void reset(planet& world)
    {
        source = &world;
        width = world.width();
        height = world.height();
        populated.fill(false);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            pressure[season].clear();
            upperheight[season].clear();
        }
    }

    bool completeFor(const planet& world) const
    {
        return source == &world && width == world.width() && height == world.height() &&
            std::all_of(populated.begin(), populated.end(), [](bool value) { return value; });
    }

    void clear()
    {
        source = nullptr;
        width = -1;
        height = -1;
        populated.fill(false);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            pressure[season].clear();
            upperheight[season].clear();
        }
    }
};

circulationfloatcache circulationcache;

planet* processfieldworld = nullptr;
std::array<climatehydrology::SeasonalProcessFields, CLIMATESEASONCOUNT> processfields;
planet* oceanfieldworld = nullptr;
std::array<floatgrid, CLIMATESEASONCOUNT> oceansstfields;
std::array<floatgrid, CLIMATESEASONCOUNT> oceanskinfields, oceanicefields;
std::array<bool, CLIMATESEASONCOUNT> oceanaccepted{};
std::array<bool, CLIMATESEASONCOUNT> atmosphereaccepted{};
struct hydrologystoragecache
{
    planet* source = nullptr;
    int columns = 0, rows = 0;
    hydrologyfloatgrid boundary, free, soil, snow;
};
hydrologystoragecache climatestorage;

struct hydrologyforcingfloatcache
{
    planet* source = nullptr;
    int width = -1;
    int height = -1;
    std::array<floatgrid, CLIMATESEASONCOUNT> surfacewindu;
    std::array<floatgrid, CLIMATESEASONCOUNT> surfacewindv;
    std::array<floatgrid, CLIMATESEASONCOUNT> upperwindu;
    std::array<floatgrid, CLIMATESEASONCOUNT> upperwindv;
    std::array<bool, CLIMATESEASONCOUNT> populated{};

    void reset(planet& world)
    {
        source = &world;
        width = world.width();
        height = world.height();
        populated.fill(false);
    }

    bool completeFor(const planet& world) const
    {
        return source == &world && width == world.width() && height == world.height() &&
            std::all_of(populated.begin(), populated.end(), [](bool value) { return value; });
    }

    void clear()
    {
        source = nullptr;
        width = -1;
        height = -1;
        populated.fill(false);

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            surfacewindu[season].clear();
            surfacewindv[season].clear();
            upperwindu[season].clear();
            upperwindv[season].clear();
        }
    }
};

hydrologyforcingfloatcache hydrologyforcingcache;

class stationaryspectralprojector
{
public:
    stationaryspectralprojector(
        int longitudecellcount,
        int latitudecellcount,
        int maximumzonalwavenumber,
        int maximummeridionalwavenumber)
        : longitudecellcount_(longitudecellcount),
          latitudecellcount_(latitudecellcount),
          maximumzonalwavenumber_(std::clamp(
              maximumzonalwavenumber,
              0,
              std::max(0, longitudecellcount / 2 - 1))),
          maximummeridionalwavenumber_(std::clamp(
              maximummeridionalwavenumber,
              0,
              std::max(0, latitudecellcount - 1))),
          zonalcosine_(
              maximumzonalwavenumber_ + 1,
              vector<float>(longitudecellcount_, 0.0f)),
          zonalsine_(
              maximumzonalwavenumber_ + 1,
              vector<float>(longitudecellcount_, 0.0f)),
          meridionalcosine_(
              maximummeridionalwavenumber_ + 1,
              vector<float>(latitudecellcount_, 0.0f))
    {
        constexpr float pi = 3.14159265358979323846f;

        for (int wavenumber = 1; wavenumber <= maximumzonalwavenumber_; wavenumber++)
        {
            for (int x = 0; x < longitudecellcount_; x++)
            {
                const float phase = 2.0f * pi * static_cast<float>(wavenumber * x) /
                    static_cast<float>(longitudecellcount_);
                zonalcosine_[wavenumber][x] = std::cos(phase);
                zonalsine_[wavenumber][x] = std::sin(phase);
            }
        }

        for (int wavenumber = 1; wavenumber <= maximummeridionalwavenumber_; wavenumber++)
        {
            for (int y = 0; y < latitudecellcount_; y++)
            {
                const float phase = pi * static_cast<float>(wavenumber) *
                    (static_cast<float>(y) + 0.5f) /
                    static_cast<float>(latitudecellcount_);
                meridionalcosine_[wavenumber][y] = std::cos(phase);
            }
        }
    }

    void project(floatgrid& field) const
    {
        if (longitudecellcount_ <= 0 || latitudecellcount_ <= 0 ||
            field.size() != static_cast<size_t>(longitudecellcount_) ||
            field.front().size() != static_cast<size_t>(latitudecellcount_))
        {
            return;
        }

        const float inverselongitudecellcount =
            1.0f / static_cast<float>(longitudecellcount_);
        vector<float> zonalmean(latitudecellcount_, 0.0f);
        vector<vector<float>> cosinecoefficient(
            maximumzonalwavenumber_ + 1,
            vector<float>(latitudecellcount_, 0.0f));
        vector<vector<float>> sinecoefficient(
            maximumzonalwavenumber_ + 1,
            vector<float>(latitudecellcount_, 0.0f));

        parallelforrows(0, latitudecellcount_ - 1, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                double mean = 0.0;

                for (int x = 0; x < longitudecellcount_; x++)
                    mean += field[x][y];

                zonalmean[y] = static_cast<float>(mean * inverselongitudecellcount);

                for (int wavenumber = 1; wavenumber <= maximumzonalwavenumber_; wavenumber++)
                {
                    double cosinesum = 0.0;
                    double sinesum = 0.0;

                    for (int x = 0; x < longitudecellcount_; x++)
                    {
                        cosinesum += static_cast<double>(field[x][y]) *
                            zonalcosine_[wavenumber][x];
                        sinesum += static_cast<double>(field[x][y]) *
                            zonalsine_[wavenumber][x];
                    }

                    cosinecoefficient[wavenumber][y] = static_cast<float>(
                        2.0 * cosinesum * inverselongitudecellcount);
                    sinecoefficient[wavenumber][y] = static_cast<float>(
                        2.0 * sinesum * inverselongitudecellcount);
                }
            }
        });

        projectmeridionally(zonalmean);

        for (int wavenumber = 1; wavenumber <= maximumzonalwavenumber_; wavenumber++)
        {
            projectmeridionally(cosinecoefficient[wavenumber]);
            projectmeridionally(sinecoefficient[wavenumber]);
        }

        parallelforrows(0, latitudecellcount_ - 1, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x < longitudecellcount_; x++)
                {
                    float projected = zonalmean[y];

                    for (int wavenumber = 1;
                        wavenumber <= maximumzonalwavenumber_;
                        wavenumber++)
                    {
                        projected += cosinecoefficient[wavenumber][y] *
                                zonalcosine_[wavenumber][x] +
                            sinecoefficient[wavenumber][y] *
                                zonalsine_[wavenumber][x];
                    }

                    field[x][y] = projected;
                }
            }
        });
    }

private:
    void projectmeridionally(vector<float>& values) const
    {
        if (values.size() != static_cast<size_t>(latitudecellcount_))
            return;

        const double inversecellcount = 1.0 / static_cast<double>(latitudecellcount_);
        vector<float> coefficient(maximummeridionalwavenumber_ + 1, 0.0f);
        double mean = 0.0;

        for (float value : values)
            mean += value;

        coefficient[0] = static_cast<float>(mean * inversecellcount);

        for (int wavenumber = 1;
            wavenumber <= maximummeridionalwavenumber_;
            wavenumber++)
        {
            double sum = 0.0;

            for (int y = 0; y < latitudecellcount_; y++)
                sum += static_cast<double>(values[y]) * meridionalcosine_[wavenumber][y];

            coefficient[wavenumber] = static_cast<float>(2.0 * sum * inversecellcount);
        }

        for (int y = 0; y < latitudecellcount_; y++)
        {
            float projected = coefficient[0];

            for (int wavenumber = 1;
                wavenumber <= maximummeridionalwavenumber_;
                wavenumber++)
            {
                projected += coefficient[wavenumber] * meridionalcosine_[wavenumber][y];
            }

            values[y] = projected;
        }
    }

    int longitudecellcount_ = 0;
    int latitudecellcount_ = 0;
    int maximumzonalwavenumber_ = 0;
    int maximummeridionalwavenumber_ = 0;
    vector<vector<float>> zonalcosine_;
    vector<vector<float>> zonalsine_;
    vector<vector<float>> meridionalcosine_;
};

int wrapx(int x, int width)
{
    if (x < 0 || x > width)
        return wrap(x, width);

    return x;
}

float latitudeforrow(int y, int height)
{
    constexpr double radiansToDegrees = 57.2957795130823208768;
    return static_cast<float>(climategrid::latitudeCentreRadians(
        y, height + 1, climategrid::LatitudeLayout::poleInclusive) *
        radiansToDegrees);
}

double rowareaweight(int y, int height)
{
    if (height <= 0)
        return 0.0;

    constexpr double pi = 3.14159265358979323846;
    const double latitude = static_cast<double>(latitudeforrow(y, height)) * pi / 180.0;
    const double halfstep = 0.5 * pi / static_cast<double>(height);
    const double north = std::min(0.5 * pi, latitude + halfstep);
    const double south = std::max(-0.5 * pi, latitude - halfstep);

    return (std::sin(north) - std::sin(south)) / (2.0 * std::sin(halfstep));
}

void smoothconvergencefootprint(
    const hydrologyfloatgrid& source,
    hydrologyfloatgrid& destination,
    hydrologyfloatgrid& scratch,
    const vector<float>& rowareas,
    int width,
    int height,
    float mixingfraction,
    int passes)
{
    const float mixing = std::clamp(mixingfraction, 0.0f, 1.0f);

    for (int pass = 0; pass < std::max(1, passes); pass++)
    {
        const hydrologyfloatgrid& current = pass == 0 ? source : destination;

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    scratch[x][y] =
                        (1.0f - mixing) * current[x][y] +
                        0.5f * mixing *
                            (current[wrapx(x - 1, width)][y] +
                                current[wrapx(x + 1, width)][y]);
                }
            }
        });

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                const float centrearea = rowareas[y];

                for (int x = 0; x <= width; x++)
                {
                    const float centre = scratch[x][y];
                    float smoothed = centre;

                    if (y > 0)
                    {
                        smoothed += mixing * rowareas[y - 1] /
                            std::max(1.0e-6f, centrearea + rowareas[y - 1]) *
                            (scratch[x][y - 1] - centre);
                    }

                    if (y < height)
                    {
                        smoothed += mixing * rowareas[y + 1] /
                            std::max(1.0e-6f, centrearea + rowareas[y + 1]) *
                            (scratch[x][y + 1] - centre);
                    }

                    destination[x][y] = smoothed;
                }
            }
        });
    }
}

void exchangeeddymoisture(
    float& first,
    float& second,
    double firstarea,
    double secondarea,
    float exchangefraction);

void damphighwavenumbers(
    hydrologyfloatgrid& field,
    hydrologyfloatgrid& laplacian,
    hydrologyfloatgrid& next,
    const vector<float>& rowareas,
    const vector<float>& meridionalfactors,
    int width,
    int height,
    float strength,
    int passes,
    bool nonnegative)
{
    (void)laplacian;
    (void)nonnegative;
    const float damping = std::clamp(strength, 0.0f, 0.20f);

    for (int pass = 0; pass < std::max(0, passes); pass++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float zonalmean = 0.5f * (
                        field[wrapx(x - 1, width)][y] +
                        field[wrapx(x + 1, width)][y]);
                    next[x][y] = field[x][y] +
                        damping * (zonalmean - field[x][y]);
                }
            }
        });

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    field[x][y] = next[x][y];
            }
        });

        parallelforrows(0, width, [&](int startcolumn, int endcolumn)
        {
            for (int x = startcolumn; x <= endcolumn; x++)
            {
                for (int parity = 0; parity < 2; parity++)
                {
                    for (int y = parity; y < height; y += 2)
                    {
                        const float edgefactor = std::min(
                            meridionalfactors[y],
                            meridionalfactors[y + 1]);
                        exchangeeddymoisture(
                            field[x][y],
                            field[x][y + 1],
                            rowareas[y],
                            rowareas[y + 1],
                            0.5f * damping * edgefactor);
                    }
                }
            }
        });
    }
}

float samplehydrologygrid(
    const hydrologyfloatgrid& field,
    int sourceWidth,
    int sourceHeight,
    int targetX,
    int targetY,
    int targetWidth,
    int targetHeight)
{
    const float sourceX = (static_cast<float>(targetX) + 0.5f) *
        static_cast<float>(sourceWidth + 1) /
        static_cast<float>(targetWidth + 1) - 0.5f;
    const int baseX = static_cast<int>(std::floor(sourceX));
    const float fractionX = sourceX - std::floor(sourceX);
    const int x0 = wrapx(baseX, sourceWidth);
    const int x1 = wrapx(baseX + 1, sourceWidth);
    const float sourceY = targetHeight > 0
        ? static_cast<float>(targetY * sourceHeight) / static_cast<float>(targetHeight)
        : 0.0f;
    const int y0 = std::clamp(static_cast<int>(std::floor(sourceY)), 0, sourceHeight);
    const int y1 = std::min(sourceHeight, y0 + 1);
    const float fractionY = sourceY - std::floor(sourceY);
    const float north = field[x0][y0] * (1.0f - fractionX) + field[x1][y0] * fractionX;
    const float south = field[x0][y1] * (1.0f - fractionX) + field[x1][y1] * fractionX;
    return north * (1.0f - fractionY) + south * fractionY;
}

float eddyexchangefraction(
    float diffusivitym2s,
    double edgelengthmetres,
    double celldistancemetres,
    double firstcellareametres2,
    double secondcellareametres2,
    float timestepseconds)
{
    if (diffusivitym2s <= 0.0f || edgelengthmetres <= 0.0 || celldistancemetres <= 0.0 ||
        firstcellareametres2 <= 0.0 || secondcellareametres2 <= 0.0 || timestepseconds <= 0.0f)
    {
        return 0.0f;
    }

    const double coupling = static_cast<double>(diffusivitym2s) *
        edgelengthmetres / celldistancemetres;
    const double decayrate = coupling *
        (1.0 / firstcellareametres2 + 1.0 / secondcellareametres2);
    return static_cast<float>(std::clamp(
        1.0 - std::exp(-decayrate * static_cast<double>(timestepseconds)),
        0.0,
        1.0));
}

void exchangeeddymoisture(
    float& first,
    float& second,
    double firstarea,
    double secondarea,
    float exchangefraction)
{
    if (exchangefraction <= 0.0f || firstarea <= 0.0 || secondarea <= 0.0)
        return;

    const double equilibrium =
        (static_cast<double>(first) * firstarea + static_cast<double>(second) * secondarea) /
        (firstarea + secondarea);
    first = static_cast<float>(first + (equilibrium - first) * exchangefraction);
    second = static_cast<float>(second + (equilibrium - second) * exchangefraction);
}

void createeddyexchangefractions(
    int width,
    int height,
    const vector<float>& latitudes,
    const vector<float>& areaweights,
    const vector<float>& meridionaltransportfactors,
    const hydrologyfloatgrid& surfacewindu,
    const hydrologyfloatgrid& surfacewindv,
    const hydrologyfloatgrid& upperwindu,
    const hydrologyfloatgrid& upperwindv,
    float timestepseconds,
    hydrologyfloatgrid& zonalfractions,
    hydrologyfloatgrid& meridionalfractions)
{
    constexpr double pi = 3.14159265358979323846;
    const double radius = tuning::climate::moistureadvection::referencePlanetRadiusMetres;
    const double equatorialcellwidth = 2.0 * pi * radius / static_cast<double>(width + 1);
    const double meridionalcellheight = pi * radius / static_cast<double>(height + 1);
    const double equatorialcellarea = equatorialcellwidth * meridionalcellheight;
    hydrologyfloatgrid diffusivity(width + 1, height + 1, 0.0f);
    vector<double> rowareas(height + 1, 0.0);
    vector<double> zonalcelldistances(height + 1, 0.0);

    for (int y = 0; y <= height; y++)
    {
        constexpr double degreestoradians = pi / 180.0;
        const float latitude = latitudes[y];
        rowareas[y] = equatorialcellarea * areaweights[y];
        zonalcelldistances[y] = equatorialcellwidth *
            std::max(0.05, std::abs(std::cos(static_cast<double>(latitude) * degreestoradians)));

        for (int x = 0; x <= width; x++)
        {
            diffusivity[x][y] = std::min(
                tuning::climate::moistureadvection::maximumTransientEddyDiffusivityM2S,
                tuning::climate::moistureadvection::backgroundMoistureDiffusivityM2S +
                    climatephysics::transientEddyDiffusivityM2S(
                        surfacewindu[x][y],
                        surfacewindv[x][y],
                        upperwindu[x][y],
                        upperwindv[x][y],
                        latitude,
                        tuning::climate::moistureadvection::transientEddyMixingLengthMetres,
                        tuning::climate::moistureadvection::maximumTransientEddyDiffusivityM2S,
                        tuning::climate::moistureadvection::transientEddyMinimumLatitudeDegrees,
                        tuning::climate::moistureadvection::transientEddyFullStrengthLatitudeDegrees));
        }
    }

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                const int east = x < width ? x + 1 : 0;
                const float edgediffusivity = 0.5f * (diffusivity[x][y] + diffusivity[east][y]);
                zonalfractions[x][y] = eddyexchangefraction(
                    edgediffusivity,
                    meridionalcellheight,
                    zonalcelldistances[y],
                    rowareas[y],
                    rowareas[y],
                    timestepseconds);
            }
        }
    });

    parallelforrows(0, width, [&](int startcolumn, int endcolumn)
    {
        for (int x = startcolumn; x <= endcolumn; x++)
        {
            for (int y = 0; y < height; y++)
            {
                const float edgediffusivity = 0.5f * (diffusivity[x][y] + diffusivity[x][y + 1]);
                const double edgelength = equatorialcellwidth *
                    0.5 * (areaweights[y] + areaweights[y + 1]);
                meridionalfractions[x][y] = eddyexchangefraction(
                    edgediffusivity,
                    edgelength,
                    meridionalcellheight,
                    rowareas[y],
                    rowareas[y + 1],
                    timestepseconds) *
                    std::min(
                        meridionaltransportfactors[y],
                        meridionaltransportfactors[y + 1]);
            }
        }
    });
}

void applyeddyexchange(
    int width,
    int height,
    const vector<float>& areaweights,
    hydrologyfloatgrid& moisture,
    const hydrologyfloatgrid& zonalfractions,
    const hydrologyfloatgrid& meridionalfractions)
{
    vector<double> rowareas(height + 1, 0.0);

    for (int y = 0; y <= height; y++)
        rowareas[y] = areaweights[y];

    auto applyzonal = [&](bool halfstep)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int parity = 0; parity < 2; parity++)
                {
                    for (int x = parity; x < width; x += 2)
                    {
                        float fraction = zonalfractions[x][y];

                        if (halfstep)
                            fraction = 1.0f - std::sqrt(std::max(0.0f, 1.0f - fraction));

                        exchangeeddymoisture(
                            moisture[x][y],
                            moisture[x + 1][y],
                            rowareas[y],
                            rowareas[y],
                            fraction);
                    }
                }

                float wrapfraction = zonalfractions[width][y];

                if (halfstep)
                    wrapfraction = 1.0f - std::sqrt(std::max(0.0f, 1.0f - wrapfraction));

                exchangeeddymoisture(
                    moisture[width][y],
                    moisture[0][y],
                    rowareas[y],
                    rowareas[y],
                    wrapfraction);
            }
        });
    };

    applyzonal(true);

    parallelforrows(0, width, [&](int startcolumn, int endcolumn)
    {
        for (int parity = 0; parity < 2; parity++)
        {
            for (int y = 1 + parity; y < height - 1; y += 2)
            {
                for (int x = startcolumn; x <= endcolumn; x++)
                {
                    exchangeeddymoisture(
                        moisture[x][y],
                        moisture[x][y + 1],
                        rowareas[y],
                        rowareas[y + 1],
                        meridionalfractions[x][y]);
                }
            }
        }
    });

    applyzonal(true);
}

float coastalweight(int distance, int maxdistance)
{
    if (distance <= 0 || distance > maxdistance)
        return 0.0f;

    return 1.0f - (static_cast<float>(distance - 1) / static_cast<float>(maxdistance));
}

bool landindir(planet& world, int x, int y, int dx, int dy, int maxdistance, int& nearestdistance)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    nearestdistance = maxdistance + 1;

    for (int distance = 1; distance <= maxdistance; distance++)
    {
        const int xx = wrapx(x + dx * distance, width);
        const int yy = y + dy * distance;

        if (yy < 0 || yy > height)
            break;

        if (world.nom(xx, yy) > sealevel)
        {
            nearestdistance = distance;
            return true;
        }
    }

    return false;
}

float sampleoceanfield(const floatgrid& field, planet& world, int x, int y)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    if (y >= 0 && y <= height && world.nom(x, y) <= sealevel)
        return field[x][y];

    for (int radius = 1; radius <= 2; radius++)
    {
        for (int dy = -radius; dy <= radius; dy++)
        {
            const int yy = y + dy;

            if (yy < 0 || yy > height)
                continue;

            for (int dx = -radius; dx <= radius; dx++)
            {
                const int xx = wrapx(x + dx, width);

                if (world.nom(xx, yy) <= sealevel)
                    return field[xx][yy];
            }
        }
    }

    return field[wrapx(x, width)][std::clamp(y, 0, height)];
}

void smoothallfield(planet& world, floatgrid& field, int iterations);

float samplewrappedfield(const floatgrid& field, planet& world, float x, float y)
{
    const int width = world.width();
    const int height = world.height();
    const float span = static_cast<float>(width + 1);

    while (x < 0.0f)
        x += span;

    while (x >= span)
        x -= span;

    y = std::clamp(y, 0.0f, static_cast<float>(height));

    const int x0 = static_cast<int>(std::floor(x));
    const int y0 = static_cast<int>(std::floor(y));
    const int x1 = wrapx(x0 + 1, width);
    const int y1 = std::min(y0 + 1, height);
    const float fracx = x - static_cast<float>(x0);
    const float fracy = y - static_cast<float>(y0);
    const float v00 = field[x0][y0];
    const float v10 = field[x1][y0];
    const float v01 = field[x0][y1];
    const float v11 = field[x1][y1];

    return
        v00 * (1.0f - fracx) * (1.0f - fracy) +
        v10 * fracx * (1.0f - fracy) +
        v01 * (1.0f - fracx) * fracy +
        v11 * fracx * fracy;
}

void smoothseasonalfield(planet& world, floatgrid& field, int iterations)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();

    floatgrid scratch = field;

    for (int iteration = 0; iteration < iterations; iteration++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        scratch[x][y] = 0.0f;
                        continue;
                    }

                    float total = 0.0f;
                    float weighttotal = 0.0f;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        const int yy = y + dy;

                        if (yy < 0 || yy > height)
                            continue;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int xx = wrapx(x + dx, width);

                            if (world.nom(xx, yy) > sealevel)
                                continue;

                            const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                            total += field[xx][yy] * weight;
                            weighttotal += weight;
                        }
                    }

                    scratch[x][y] = (weighttotal > 0.0f) ? total / weighttotal : field[x][y];
                }
            }
        });

        field.swap(scratch);
    }
}

void applytopographicwindeffects(planet& world, const floatgrid& macroterrain, floatgrid& windu, floatgrid& windv)
{
    const int width = world.width();
    const int height = world.height();
    const float maxvectorwind = tuning::climate::atmosphere::maxVectorWind;

    auto computegradient = [&](int x, int y)
    {
        const int xwest = wrapx(x - 1, width);
        const int xeast = wrapx(x + 1, width);
        const int ynorth = (y > 0) ? y - 1 : y;
        const int ysouth = (y < height) ? y + 1 : y;

        const float gradx = (macroterrain[xeast][y] - macroterrain[xwest][y]) / 2.0f;
        const float grady = (macroterrain[x][ysouth] - macroterrain[x][ynorth]) / 2.0f;

        return std::pair<float, float>(gradx, grady);
    };

    for (int iteration = 0; iteration < tuning::climate::atmosphere::topographyIterations; iteration++)
    {
        floatgrid nextu = windu;
        floatgrid nextv = windv;

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float u = windu[x][y];
                    const float v = windv[x][y];
                    const float speed = std::sqrt(u * u + v * v);

                    if (speed < tuning::climate::atmosphere::topographyMinimumWindSpeed)
                        continue;

                    const auto [gradx, grady] = computegradient(x, y);
                    const float gradmag = std::sqrt(gradx * gradx + grady * grady);
                    const float terrainhere = macroterrain[x][y];

                    if (gradmag < tuning::climate::atmosphere::topographyMinimumRelief && terrainhere < tuning::climate::atmosphere::topographyMinimumRelief)
                        continue;

                    const float nx = gradx / std::max(gradmag, 0.0001f);
                    const float ny = grady / std::max(gradmag, 0.0001f);
                    const float tx = -ny;
                    const float ty = nx;
                    const float dirx = u / speed;
                    const float diry = v / speed;

                    float lookaheadrise = 0.0f;

                    for (int step = 1; step <= tuning::climate::atmosphere::topographyLookaheadDistance; step++)
                    {
                        const float samplex = static_cast<float>(x) + dirx * static_cast<float>(step);
                        const float sampley = static_cast<float>(y) + diry * static_cast<float>(step);
                        const float rise = samplewrappedfield(macroterrain, world, samplex, sampley) - terrainhere;

                        if (rise > lookaheadrise)
                            lookaheadrise = rise;
                    }

                    const float positivecontour = samplewrappedfield(macroterrain, world,
                        static_cast<float>(x) + tx * tuning::climate::atmosphere::topographySideSampleDistance,
                        static_cast<float>(y) + ty * tuning::climate::atmosphere::topographySideSampleDistance);
                    const float negativecontour = samplewrappedfield(macroterrain, world,
                        static_cast<float>(x) - tx * tuning::climate::atmosphere::topographySideSampleDistance,
                        static_cast<float>(y) - ty * tuning::climate::atmosphere::topographySideSampleDistance);

                    const float steeringdirection = (positivecontour <= negativecontour) ? 1.0f : -1.0f;
                    const float gradientbarrier = std::clamp(gradmag / tuning::climate::atmosphere::topographyGradientScale, 0.0f, 1.0f);
                    const float lookaheadbarrier = std::clamp(lookaheadrise / tuning::climate::atmosphere::topographyLookaheadRiseScale, 0.0f, 1.0f);
                    const float barrier = std::max(gradientbarrier, lookaheadbarrier);

                    if (barrier <= 0.0f)
                        continue;

                    float crossridge = u * nx + v * ny;
                    float alongridge = u * tx + v * ty;

                    if (crossridge > 0.0f)
                    {
                        const float blockedfraction = barrier * (1.0f - tuning::climate::atmosphere::blockedComponentFactor);
                        const float blockedcross = crossridge * blockedfraction;

                        crossridge = crossridge - blockedcross;
                        alongridge = alongridge + blockedcross * tuning::climate::atmosphere::topographyDeflectionFactor * steeringdirection;
                        crossridge = crossridge * (1.0f - barrier * tuning::climate::atmosphere::topographyChannelFactor);
                    }
                    else if (crossridge < 0.0f)
                    {
                        crossridge = crossridge * (1.0f + barrier * tuning::climate::atmosphere::topographyDownslopeAcceleration);
                    }

                    float newu = tx * alongridge + nx * crossridge;
                    float newv = ty * alongridge + ny * crossridge;
                    const float newspeed = std::sqrt(newu * newu + newv * newv);

                    if (newspeed > 0.0f)
                    {
                        const float roughnessdrag = 1.0f - barrier * tuning::climate::atmosphere::topographySpeedReduction;
                        const float dragfactor = std::max(0.0f, roughnessdrag);
                        const float cappedspeed = std::min(maxvectorwind, newspeed * dragfactor);
                        const float speedscale = cappedspeed / newspeed;

                        newu = newu * speedscale;
                        newv = newv * speedscale;
                    }

                    nextu[x][y] = std::clamp(newu, -maxvectorwind, maxvectorwind);
                    nextv[x][y] = std::clamp(newv, -maxvectorwind, maxvectorwind);
                }
            }
        });

        windu.swap(nextu);
        windv.swap(nextv);
        smoothallfield(world, windu, 1);
        smoothallfield(world, windv, 1);
    }
}

void storeterrainverticalmotion(planet& world, int season, const floatgrid& macroterrain, const floatgrid& windu, const floatgrid& windv)
{
    const int width = world.width();
    const int height = world.height();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const int ynorth = (y > 0) ? y - 1 : y;
            const int ysouth = (y < height) ? y + 1 : y;

            for (int x = 0; x <= width; x++)
            {
                const float u = windu[x][y];
                const float v = windv[x][y];
                const float speed = std::sqrt(u * u + v * v);

                if (speed < tuning::climate::atmosphere::topographyMinimumWindSpeed)
                {
                    world.setseasonaluplift(season, x, y, 0);
                    world.setseasonalsubsidence(season, x, y, 0);
                    continue;
                }

                const int xwest = wrapx(x - 1, width);
                const int xeast = wrapx(x + 1, width);
                const float latitude = latitudeforrow(y, height);
                const auto spacing = climateatmosphere::cellSpacingMetres(
                    latitude,
                    width + 1,
                    height + 1,
                    tuning::climate::atmosphere::referencePlanetRadiusMetres);
                const float elevationchangeeast =
                    (macroterrain[xeast][y] - macroterrain[xwest][y]) / 2.0f;
                const float elevationchangesouth =
                    (macroterrain[x][ysouth] - macroterrain[x][ynorth]) / 2.0f;
                const float gradmag = std::sqrt(
                    elevationchangeeast * elevationchangeeast +
                    elevationchangesouth * elevationchangesouth);
                const float terrainhere = macroterrain[x][y];

                if (gradmag < tuning::climate::atmosphere::topographyMinimumRelief && terrainhere < tuning::climate::atmosphere::topographyMinimumRelief)
                {
                    world.setseasonaluplift(season, x, y, 0);
                    world.setseasonalsubsidence(season, x, y, 0);
                    continue;
                }

                const float slopeeast = elevationchangeeast / spacing.zonalMetres;
                const float slopesouth = elevationchangesouth / spacing.meridionalMetres;
                const float verticalvelocity = u * slopeeast + v * slopesouth;
                const float directioneast = u / speed;
                const float directionsouth = v / speed;
                float windwardrise = 0.0f;
                float leewarddrop = 0.0f;

                for (int step = 1;
                    step <= tuning::climate::atmosphere::topographyLookaheadDistance;
                    step++)
                {
                    const float samplex = static_cast<float>(x) +
                        directioneast * static_cast<float>(step);
                    const float sampley = static_cast<float>(y) +
                        directionsouth * static_cast<float>(step);
                    const float elevationchange =
                        samplewrappedfield(macroterrain, world, samplex, sampley) - terrainhere;
                    windwardrise = std::max(windwardrise, elevationchange);
                    leewarddrop = std::max(leewarddrop, -elevationchange);
                }

                const float inversecellcrossingtime = std::sqrt(
                    (u / spacing.zonalMetres) * (u / spacing.zonalMetres) +
                    (v / spacing.meridionalMetres) * (v / spacing.meridionalMetres));
                const float cellcrossingtime = inversecellcrossingtime > 0.0f ?
                    1.0f / inversecellcrossingtime : 0.0f;
                float parceldisplacement = verticalvelocity * cellcrossingtime;

                if (parceldisplacement > 0.0f)
                    parceldisplacement = std::max(parceldisplacement, windwardrise);
                else if (parceldisplacement < 0.0f)
                    parceldisplacement = std::min(parceldisplacement, -leewarddrop);

                parceldisplacement = std::clamp(
                    parceldisplacement,
                    -tuning::climate::atmosphere::maximumOrographicParcelDisplacementMetres,
                    tuning::climate::atmosphere::maximumOrographicParcelDisplacementMetres);
                const float uplift = std::max(0.0f, parceldisplacement);
                const float subsidence = std::max(0.0f, -parceldisplacement);

                world.setseasonaluplift(season, x, y, static_cast<int>(std::round(uplift * tuning::climate::atmosphere::topographyVerticalMotionStorageScale)));
                world.setseasonalsubsidence(season, x, y, static_cast<int>(std::round(subsidence * tuning::climate::atmosphere::topographyVerticalMotionStorageScale)));
            }
        }
    });
}
}

void createoceancurrentmap(planet& world)
{
    const int outputcolumns = world.width() + 1;
    const int outputrows = world.height() + 1;
    const int sealevel = world.sealevel();
    const auto dimensions = climatehydrology::climateGridDimensions(
        outputcolumns,
        outputrows,
        tuning::climate::oceancurrents::internalHorizontalCells);
    const int columns = dimensions.columns;
    const int rows = dimensions.rows;
    const std::size_t outputcellcount = static_cast<std::size_t>(outputcolumns) * outputrows;
    const std::size_t cellcount = static_cast<std::size_t>(columns) * rows;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        climateocean::OceanConfig config;
        std::vector<float> outputland(outputcellcount, 0.0f);
        std::vector<float> outputdepth(outputcellcount, 0.0f);
        std::vector<float> outputeastwind(outputcellcount, 0.0f);
        std::vector<float> outputsouthwind(outputcellcount, 0.0f);
        std::vector<float> outputtemperature(outputcellcount, 0.0f);
        std::vector<float> outputsst(outputcellcount, 0.0f);
        std::vector<float> outputice(outputcellcount, 0.0f), outputfluxreference(outputcellcount, 0.0f);
        for (int y = 0; y < outputrows; y++)
        {
            for (int x = 0; x < outputcolumns; x++)
            {
                const std::size_t cell = static_cast<std::size_t>(y) * outputcolumns + x;
                const bool land = world.nom(x, y) > sealevel;
                outputland[cell] = land ? 1.0f : 0.0f;
                outputdepth[cell] = land
                    ? 0.0f
                    : static_cast<float>(std::max(50, sealevel - world.nom(x, y)));
                const bool cached = hydrologyforcingcache.source == &world && hydrologyforcingcache.populated[season];
                outputeastwind[cell] = cached ? hydrologyforcingcache.surfacewindu[season][x][y] : static_cast<float>(world.seasonaluwind(season, x, y));
                outputsouthwind[cell] = cached ? hydrologyforcingcache.surfacewindv[season][x][y] : static_cast<float>(world.seasonalvwind(season, x, y));
                outputtemperature[cell] = static_cast<float>(world.seasonaltemp(season, x, y));
                // The base climate supplies ice-skin temperature, not supercooled
                // liquid. Prescribe a seasonal 1 m ice reservoir where frozen.
                const float freezing = config.freezingTemperatureC;
                outputsst[cell] = land ? outputtemperature[cell] : std::max(freezing, outputtemperature[cell]);
                outputice[cell] = !land && outputtemperature[cell] < freezing ? 1.0f : 0.0f;
                outputfluxreference[cell] = oceanfieldworld == &world && !oceanskinfields[season].empty()
                    ? oceanskinfields[season][x][y] : outputtemperature[cell];
            }
        }

        const auto remapinput = [&](const std::vector<float>& field)
        {
            return climategrid::remapField(
                outputcolumns,
                outputrows,
                climategrid::LatitudeLayout::poleInclusive,
                field,
                columns,
                rows,
                climategrid::LatitudeLayout::cellCentred);
        };
        climateocean::OceanForcing forcing;
        const std::vector<float> landfraction = remapinput(outputland);
        forcing.landMask.resize(cellcount, 0);
        for (std::size_t cell = 0; cell < cellcount; cell++)
            forcing.landMask[cell] = landfraction[cell] >= 0.5f ? 1 : 0;
        forcing.bathymetryMetres = remapinput(outputdepth);
        forcing.eastWindMps = remapinput(outputeastwind);
        forcing.southWindMps = remapinput(outputsouthwind);
        forcing.atmosphericTemperatureC = remapinput(outputtemperature);
        forcing.initialSstC = remapinput(outputsst);
        for (std::size_t cell = 0; cell < cellcount; ++cell)
            if (!forcing.landMask[cell]) forcing.initialSstC[cell] = std::max(config.freezingTemperatureC, forcing.initialSstC[cell]);
        forcing.initialIceThicknessMetres = remapinput(outputice);
        forcing.surfaceHeatFluxReferenceTemperatureC = remapinput(outputfluxreference);
        if (processfieldworld == &world && processfields[season].durationSeconds > 0.0)
        {
            const auto& fields = processfields[season];
            forcing.surfaceHeatFluxWm2 = climategrid::remapField(fields.columns, fields.rows,
                climategrid::LatitudeLayout::cellCentred, fields.surfaceNetHeatingWm2,
                columns, rows, climategrid::LatitudeLayout::cellCentred);
        }

        config.rotationDirection = world.rotation() ? 1.0f : -1.0f;
        config.rotationRatePerSecond = tuning::climate::atmosphere::rotationRatePerSecond;
        config.planetRadiusMetres = tuning::climate::atmosphere::referencePlanetRadiusMetres;
        config.airDensityKgM3 = tuning::climate::atmosphere::surfaceAirDensityKgM3;
        config.dragCoefficient = tuning::climate::atmosphere::oceanMomentumDragCoefficient;
        config.waterDensityKgM3 = tuning::climate::oceancurrents::waterDensityKgM3;
        config.barotropicDragPerSecond =
            tuning::climate::oceancurrents::barotropicDragPerSecond;
        config.mixedLayerDepthMetres =
            tuning::climate::oceancurrents::mixedLayerDepthMetres;
        config.heatDiffusivityM2S = tuning::climate::oceancurrents::heatDiffusivityM2S;
        config.surfaceHeatExchangeWm2K =
            tuning::climate::oceancurrents::surfaceHeatExchangeWm2K;
        config.streamfunctionIterations =
            tuning::climate::oceancurrents::streamfunctionIterations;
        config.heatStepsPerIteration =
            tuning::climate::oceancurrents::heatStepsPerCouplingIteration;
        config.couplingIterations = tuning::climate::oceancurrents::couplingIterations;
        config.underRelaxation =
            tuning::climate::oceancurrents::couplingUnderRelaxation;
        config.convergenceTolerance = tuning::climate::oceancurrents::couplingTolerance;
        config.sstWindFeedbackMpsPerK =
            tuning::climate::oceancurrents::sstWindFeedbackMpsPerK;
        config.maximumCurrentMps = tuning::climate::oceancurrents::maximumCurrentMps;
        // Atmosphere feedback is handled by the shared outer climate solve,
        // never by the ocean's legacy SST-to-pressure surrogate.
        config.oneWay = true;
        climateocean::OceanState ocean = climateocean::solveWindDrivenOcean(
            columns, rows, forcing, config);
        oceanaccepted[season] = climateocean::usableOceanState(ocean, cellcount);
        if (!oceanaccepted[season])
        {
            ocean.converged = false;
            std::cout << "Ocean fallback: season=" << season << " iterations=" << ocean.couplingIterations
                << " residual=" << ocean.relativeResidual << " basin_residual=" << ocean.streamfunctionRelativeResidual << '\n';
            ocean.sstC = forcing.initialSstC;
            ocean.surfaceSkinTemperatureC = forcing.atmosphericTemperatureC;
            ocean.iceThicknessMetres = forcing.initialIceThicknessMetres;
            ocean.eastCurrentMps.assign(cellcount, 0.0f);
            ocean.southCurrentMps.assign(cellcount, 0.0f);
            ocean.ekmanUpwellingMps.assign(cellcount, 0.0f);
            ocean.coupledEastWindMps = forcing.eastWindMps;
            ocean.coupledSouthWindMps = forcing.southWindMps;
            ocean.coupledPressureAnomalyHpa.assign(cellcount, 0.0f);
        }

        climatevalidation::captureoceanfields(world, season, columns, rows, ocean);
        const auto remapoutput = [&](const std::vector<float>& field)
        {
            return climategrid::remapField(
                columns,
                rows,
                climategrid::LatitudeLayout::cellCentred,
                field,
                outputcolumns,
                outputrows,
                climategrid::LatitudeLayout::poleInclusive);
        };
        const std::vector<float> eastcurrent = remapoutput(ocean.eastCurrentMps);
        const std::vector<float> southcurrent = remapoutput(ocean.southCurrentMps);
        auto sst = remapoutput(ocean.sstC);
        for (std::size_t cell = 0; cell < outputcellcount; ++cell)
            if (outputland[cell] == 0.0f) sst[cell] = std::max(config.freezingTemperatureC, sst[cell]);
        const auto skin = remapoutput(ocean.surfaceSkinTemperatureC);
        auto icecover = ocean.iceThicknessMetres;
        for (auto& h : icecover) h = std::clamp(h / 0.10f, 0.0f, 1.0f);
        const auto ice = remapoutput(icecover);
        oceanfieldworld = &world;
        oceansstfields[season].assign(outputcolumns, vector<float>(outputrows));
        oceanskinfields[season].assign(outputcolumns, vector<float>(outputrows));
        oceanicefields[season].assign(outputcolumns, vector<float>(outputrows));
        for (int y = 0; y < outputrows; ++y)
            for (int x = 0; x < outputcolumns; ++x)
            {
                const auto cell = y * outputcolumns + x;
                oceansstfields[season][x][y] = sst[cell];
                oceanskinfields[season][x][y] = skin[cell];
                oceanicefields[season][x][y] = ice[cell];
            }
        const std::vector<float> coupledeastwind = remapoutput(ocean.coupledEastWindMps);
        const std::vector<float> coupledsouthwind = remapoutput(ocean.coupledSouthWindMps);
        const auto pressurefeedback = remapoutput(ocean.coupledPressureAnomalyHpa);
        for (int y = 0; y < outputrows; y++)
        {
            for (int x = 0; x < outputcolumns; x++)
            {
                const std::size_t cell = static_cast<std::size_t>(y) * outputcolumns + x;
                if (world.nom(x, y) > sealevel)
                {
                    world.setseasonalcurrentu(season, x, y, 0);
                    world.setseasonalcurrentv(season, x, y, 0);
                    world.setseasonalsst(season, x, y, 0);
                    continue;
                }
                world.setseasonalcurrentu(season, x, y, static_cast<int>(std::round(
                    eastcurrent[cell] * tuning::climate::oceancurrents::
                        currentStorageCentimetresPerSecond)));
                world.setseasonalcurrentv(season, x, y, static_cast<int>(std::round(
                    southcurrent[cell] * tuning::climate::oceancurrents::
                        currentStorageCentimetresPerSecond)));
                world.setseasonalsst(season, x, y, static_cast<int>(std::round(sst[cell])));
                world.setseasonaluwind(season, x, y, static_cast<int>(std::round(
                    coupledeastwind[cell])));
                world.setseasonalvwind(season, x, y, static_cast<int>(std::round(
                    coupledsouthwind[cell])));
                world.setseasonalpressure(season, x, y, world.seasonalpressure(season, x, y) +
                    static_cast<int>(std::round(pressurefeedback[cell])));
                if (hydrologyforcingcache.source == &world &&
                    hydrologyforcingcache.populated[season])
                {
                    hydrologyforcingcache.surfacewindu[season][x][y] = coupledeastwind[cell];
                    hydrologyforcingcache.surfacewindv[season][x][y] = coupledsouthwind[cell];
                }
            }
        }
        if (isworldgendebugrunactive() && hydrologyforcingcache.source == &world && hydrologyforcingcache.populated[season])
            climatevalidation::capturecirculationwindfields(world, season,
                hydrologyforcingcache.surfacewindu[season], hydrologyforcingcache.surfacewindv[season],
                hydrologyforcingcache.upperwindu[season], hydrologyforcingcache.upperwindv[season]);
        std::cout
            << "Coupled ocean season " << season
            << ": iterations=" << ocean.couplingIterations
            << " residual=" << ocean.relativeResidual
            << " streamfunction_residual=" << ocean.streamfunctionRelativeResidual
            << " converged=" << ocean.converged
            << " max_divergence=" << ocean.maximumTransportDivergenceMps
            << " heat_residual_j=" << ocean.heatBudgetResidualJ
            << " relative_heat_residual=" << ocean.relativeHeatBudgetResidual << '\n';
    }
    if (!tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
        for (int y = 0; y < outputrows; ++y)
            for (int x = 0; x < outputcolumns; ++x)
            {
                int frozenSeasons = 0;
                for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
                    frozenSeasons += oceanicefields[season][x][y] >= 0.5f;
                world.setseaice(x, y, world.sea(x, y)
                    ? (frozenSeasons == CLIMATESEASONCOUNT ? 2 : (frozenSeasons > 0 ? 1 : 0)) : 0);
            }
}

[[maybe_unused]] void createprescribedoceancurrentmap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const int coastalsearchdistance = tuning::climate::oceancurrents::coastalSearchDistance;
    const float equatorialband = tuning::climate::oceancurrents::equatorialBand;
    const float midlatitudeband = tuning::climate::oceancurrents::midLatitudeBand;
    const float polarband = tuning::climate::oceancurrents::polarBand;
    const float countercurrentband = tuning::climate::oceancurrents::counterCurrentBand;
    const float retainedbasestrength = tuning::climate::oceancurrents::retainedBaseStrength;
    const float smoothingblend = tuning::climate::oceancurrents::smoothingBlend;
    const float blockedcomponentfactor = tuning::climate::oceancurrents::blockedComponentFactor;
    const float maxcurrentspeed = tuning::climate::oceancurrents::equatorialSpeed;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const float latitudeshift = seasonlatitudephase[season] * world.tilt() * tuning::climate::oceancurrents::seasonalShiftFactor;

        floatgrid baseu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid basev(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid currentu(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid currentv(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                        continue;

                    const float latitude = latitudeforrow(y, height);
                    const float effectivelatitude = latitude - latitudeshift;
                    const float absolutelatitude = std::fabs(effectivelatitude);

                    float u = 0.0f;
                    float v = 0.0f;

                    if (absolutelatitude <= countercurrentband)
                        u = tuning::climate::oceancurrents::counterCurrentSpeed;
                    else if (absolutelatitude < equatorialband)
                        u = -tuning::climate::oceancurrents::equatorialSpeed;
                    else if (absolutelatitude < midlatitudeband)
                        u = tuning::climate::oceancurrents::midLatitudeSpeed;
                    else
                        u = -tuning::climate::oceancurrents::polarSpeed;

                    const float polewardsign = (latitude >= 0.0f) ? -1.0f : 1.0f;

                    int westdistance = 0;
                    int eastdistance = 0;

                    const bool westland = landindir(world, x, y, -1, 0, coastalsearchdistance, westdistance);
                    const bool eastland = landindir(world, x, y, 1, 0, coastalsearchdistance, eastdistance);

                    if (absolutelatitude >= countercurrentband && absolutelatitude < midlatitudeband)
                    {
                        if (westland)
                            v += polewardsign * tuning::climate::oceancurrents::westernBoundarySpeed * coastalweight(westdistance, coastalsearchdistance);

                        if (eastland)
                            v -= polewardsign * tuning::climate::oceancurrents::easternBoundarySpeed * coastalweight(eastdistance, coastalsearchdistance);
                    }
                    else if (absolutelatitude >= midlatitudeband && absolutelatitude < polarband)
                    {
                        if (westland)
                            v -= polewardsign * tuning::climate::oceancurrents::subpolarBoundarySpeed * coastalweight(westdistance, coastalsearchdistance);

                        if (eastland)
                            v += polewardsign * tuning::climate::oceancurrents::subpolarBoundarySpeed * coastalweight(eastdistance, coastalsearchdistance);
                    }

                    baseu[x][y] = u;
                    basev[x][y] = v;
                    currentu[x][y] = u;
                    currentv[x][y] = v;
                }
            }
        });

        floatgrid nextu = currentu;
        floatgrid nextv = currentv;

        for (int iteration = 0; iteration < tuning::climate::oceancurrents::smoothingIterations; iteration++)
        {
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        if (world.nom(x, y) > sealevel)
                        {
                            nextu[x][y] = 0.0f;
                            nextv[x][y] = 0.0f;
                            continue;
                        }

                        float avgu = 0.0f;
                        float avgv = 0.0f;
                        float weighttotal = 0.0f;

                        for (int dy = -1; dy <= 1; dy++)
                        {
                            const int yy = y + dy;

                            if (yy < 0 || yy > height)
                                continue;

                            for (int dx = -1; dx <= 1; dx++)
                            {
                                const int xx = wrapx(x + dx, width);

                                if (world.nom(xx, yy) > sealevel)
                                    continue;

                                const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                                avgu += currentu[xx][yy] * weight;
                                avgv += currentv[xx][yy] * weight;
                                weighttotal += weight;
                            }
                        }

                        if (weighttotal > 0.0f)
                        {
                            avgu = avgu / weighttotal;
                            avgv = avgv / weighttotal;
                        }

                        float blendedu = currentu[x][y] * (1.0f - smoothingblend) + avgu * smoothingblend;
                        float blendedv = currentv[x][y] * (1.0f - smoothingblend) + avgv * smoothingblend;

                        blendedu = blendedu * retainedbasestrength + baseu[x][y] * (1.0f - retainedbasestrength);
                        blendedv = blendedv * retainedbasestrength + basev[x][y] * (1.0f - retainedbasestrength);

                        const bool eastblocked = world.nom(wrapx(x + 1, width), y) > sealevel;
                        const bool westblocked = world.nom(wrapx(x - 1, width), y) > sealevel;
                        const bool northblocked = (y == 0) || world.nom(x, y - 1) > sealevel;
                        const bool southblocked = (y == height) || world.nom(x, y + 1) > sealevel;

                        if (eastblocked && blendedu > 0.0f)
                            blendedu = blendedu * blockedcomponentfactor;

                        if (westblocked && blendedu < 0.0f)
                            blendedu = blendedu * blockedcomponentfactor;

                        if (northblocked && blendedv < 0.0f)
                            blendedv = blendedv * blockedcomponentfactor;

                        if (southblocked && blendedv > 0.0f)
                            blendedv = blendedv * blockedcomponentfactor;

                        nextu[x][y] = std::clamp(blendedu, -maxcurrentspeed, maxcurrentspeed);
                        nextv[x][y] = std::clamp(blendedv, -maxcurrentspeed, maxcurrentspeed);
                    }
                }
            });

            currentu.swap(nextu);
            currentv.swap(nextv);
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        world.setseasonalcurrentu(season, x, y, 0);
                        world.setseasonalcurrentv(season, x, y, 0);
                        continue;
                    }

                    world.setseasonalcurrentu(season, x, y, static_cast<int>(std::round(currentu[x][y])));
                    world.setseasonalcurrentv(season, x, y, static_cast<int>(std::round(currentv[x][y])));
                }
            }
        });
    }
}

void createsurfacetemperaturemap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        world.setseasonalsst(season, x, y, 0);
                        world.setseasonalevaporation(season, x, y, 0);
                        world.setseasonalmoisture(season, x, y, 0);
                        continue;
                    }
                    const float sst = tuning::climate::oceancurrents::oneWayDiagnosticsOnly
                        ? static_cast<float>(world.seasonaltemp(season, x, y))
                        : oceanskinfields[season][x][y];
                    const float eastwind = static_cast<float>(world.seasonaluwind(
                        season, x, y));
                    const float southwind = static_cast<float>(world.seasonalvwind(
                        season, x, y));
                    const float evaporation = climatephysics::bulkAerodynamicEvaporationMm(
                        sst,
                        0.0f,
                        std::max(1.0f, std::hypot(eastwind, southwind)),
                        0.75f,
                        tuning::climate::circulation::secondsPerDay,
                        tuning::climate::atmosphere::oceanMomentumDragCoefficient);
                    world.setseasonalevaporation(
                        season, x, y, static_cast<int>(std::round(evaporation)));
                    world.setseasonalmoisture(
                        season, x, y, static_cast<int>(std::round(evaporation)));
                }
            }
        });
    }
}

[[maybe_unused]] void createprescribedsurfacetemperaturemap(planet& world)
{
    const int width = world.width();
    const int height = world.height();
    const int sealevel = world.sealevel();
    const int coastalsearchdistance = tuning::climate::oceancurrents::coastalSearchDistance;
    const float maxcurrentspeed = tuning::climate::oceancurrents::equatorialSpeed;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        floatgrid basetemperatures(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) <= sealevel)
                        basetemperatures[x][y] = static_cast<float>(world.seasonaltemp(season, x, y));
                }
            }
        });

        smoothseasonalfield(world, basetemperatures, tuning::climate::sst::smoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    if (world.nom(x, y) > sealevel)
                    {
                        world.setseasonalsst(season, x, y, 0);
                        world.setseasonalevaporation(season, x, y, 0);
                        world.setseasonalmoisture(season, x, y, 0);
                        continue;
                    }

                    const float currentu = static_cast<float>(world.seasonalcurrentu(season, x, y));
                    const float currentv = static_cast<float>(world.seasonalcurrentv(season, x, y));
                    const float magnitude = std::sqrt(currentu * currentu + currentv * currentv);

                    const int sourcex = wrapx(x - static_cast<int>(std::round((currentu / maxcurrentspeed) * tuning::climate::sst::advectionSampleDistance)), width);
                    const int sourcey = std::clamp(y - static_cast<int>(std::round((currentv / maxcurrentspeed) * tuning::climate::sst::advectionSampleDistance)), 0, height);

                    float sst = basetemperatures[x][y];
                    const float sourcetemperature = sampleoceanfield(basetemperatures, world, sourcex, sourcey);
                    sst = sst + (sourcetemperature - sst) * tuning::climate::sst::advectionBlend;

                    const float latitude = latitudeforrow(y, height);
                    const float polewardflow = ((latitude >= 0.0f) ? -currentv : currentv) / maxcurrentspeed;
                    const float equatorwardflow = ((latitude >= 0.0f) ? currentv : -currentv) / maxcurrentspeed;

                    int westdistance = 0;
                    int eastdistance = 0;

                    const bool westland = landindir(world, x, y, -1, 0, coastalsearchdistance, westdistance);
                    const bool eastland = landindir(world, x, y, 1, 0, coastalsearchdistance, eastdistance);

                    if (westland)
                        sst += std::max(0.0f, polewardflow) * tuning::climate::sst::westernBoundaryWarming * coastalweight(westdistance, coastalsearchdistance);

                    if (eastland)
                        sst -= std::max(0.0f, equatorwardflow) * tuning::climate::sst::easternBoundaryCooling * coastalweight(eastdistance, coastalsearchdistance);

                    sst = std::clamp(sst, tuning::climate::sst::minimumSst, tuning::climate::sst::maximumSst);

                    const float evaporation = std::max(0.0f, (sst + 10.0f) * tuning::climate::sst::evaporationScale + magnitude * tuning::climate::sst::evaporationCurrentBoost);

                    world.setseasonalsst(season, x, y, static_cast<int>(std::round(sst)));
                    world.setseasonalevaporation(season, x, y, static_cast<int>(std::round(evaporation)));
                    world.setseasonalmoisture(season, x, y, static_cast<int>(std::round(evaporation)));
                }
            }
        });
    }
}

namespace
{
float pressuresurfacetemperature(planet& world, int season, int x, int y)
{
    if (tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
        return static_cast<float>(world.seasonaltemp(season, x, y));
    if (world.sea(x, y) == 1 && oceanfieldworld == &world &&
        oceansstfields[season].size() == static_cast<std::size_t>(world.width() + 1))
        return oceanskinfields[season][x][y];
    if (world.sea(x, y) == 1 && world.seasonalsst(season, x, y) != 0)
        return static_cast<float>(world.seasonalsst(season, x, y));

    return static_cast<float>(world.seasonaltemp(season, x, y));
}

void smoothallfield(planet& world, floatgrid& field, int iterations)
{
    const int width = world.width();
    const int height = world.height();
    floatgrid scratch = field;

    for (int iteration = 0; iteration < iterations; iteration++)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    float total = 0.0f;
                    float weighttotal = 0.0f;

                    for (int dy = -1; dy <= 1; dy++)
                    {
                        const int yy = y + dy;

                        if (yy < 0 || yy > height)
                            continue;

                        for (int dx = -1; dx <= 1; dx++)
                        {
                            const int xx = wrapx(x + dx, width);
                            const float weight = (dx == 0 && dy == 0) ? 2.0f : 1.0f;
                            total += field[xx][yy] * weight;
                            weighttotal += weight;
                        }
                    }

                    scratch[x][y] = (weighttotal > 0.0f) ? total / weighttotal : field[x][y];
                }
            }
        });

        field.swap(scratch);
    }
}

floatgrid buildcontinentalityfield(planet& world, int smoothingiterations, float exponent)
{
    const int width = world.width();
    const int height = world.height();
    floatgrid continentality(width + 1, vector<float>(height + 1, 0.0f));

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
                continentality[x][y] = world.sea(x, y) == 1 ? 0.0f : 1.0f;
        }
    });

    smoothallfield(world, continentality, smoothingiterations);

    if (std::fabs(exponent - 1.0f) > 0.001f)
    {
        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    continentality[x][y] = std::pow(std::clamp(continentality[x][y], 0.0f, 1.0f), exponent);
            }
        });
    }

    return continentality;
}

climateatmosphere::CirculationPrecisionStageDiagnostics measurecirculationprecision(
    planet& world,
    const floatgrid& pressure,
    const floatgrid& upperheight,
    const floatgrid& continentality,
    const floatgrid& macroterrain)
{
    const int width = world.width();
    const int height = world.height();
    double weighttotal = 0.0;
    double floatpressuregradientsquared = 0.0;
    double roundedpressuregradientsquared = 0.0;
    double pressuregradientdifferencesquared = 0.0;
    double floatwindsquared = 0.0;
    double roundedwindsquared = 0.0;
    double winddifferencesquared = 0.0;
    double floatheightgradientsquared = 0.0;
    double roundedheightgradientsquared = 0.0;
    double heightgradientdifferencesquared = 0.0;
    double floatupperwindsquared = 0.0;
    double roundedupperwindsquared = 0.0;
    double upperwinddifferencesquared = 0.0;

    for (int y = 0; y <= height; y++)
    {
        const int ynorth = y > 0 ? y - 1 : y;
        const int ysouth = y < height ? y + 1 : y;
        const float latitude = latitudeforrow(y, height);
        const double weight = rowareaweight(y, height);
        const auto spacing = climateatmosphere::cellSpacingMetres(
            latitude,
            width + 1,
            height + 1,
            tuning::climate::atmosphere::referencePlanetRadiusMetres);

        for (int x = 0; x <= width; x++)
        {
            const int xwest = wrapx(x - 1, width);
            const int xeast = wrapx(x + 1, width);
            const float floatpressureeast =
                (pressure[xeast][y] - pressure[xwest][y]) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.zonalMetres);
            const float floatpressurenorth =
                (pressure[x][ynorth] - pressure[x][ysouth]) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.meridionalMetres);
            const float roundedpressureeast =
                (std::round(pressure[xeast][y]) - std::round(pressure[xwest][y])) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.zonalMetres);
            const float roundedpressurenorth =
                (std::round(pressure[x][ynorth]) - std::round(pressure[x][ysouth])) *
                tuning::climate::atmosphere::pressurePascalsPerHectopascal /
                (2.0f * spacing.meridionalMetres);
            const float continental = std::clamp(continentality[x][y], 0.0f, 1.0f);
            const float relief = std::clamp(macroterrain[x][y] / 4500.0f, 0.0f, 1.0f);
            const float basedragcoefficient =
                tuning::climate::atmosphere::oceanMomentumDragCoefficient +
                continental * (tuning::climate::atmosphere::landMomentumDragCoefficient -
                    tuning::climate::atmosphere::oceanMomentumDragCoefficient);
            const float dragcoefficient = basedragcoefficient + relief *
                (tuning::climate::atmosphere::highReliefMomentumDragCoefficient -
                    basedragcoefficient);
            const auto floatwind = climateatmosphere::steadyQuadraticDragCoriolisWind(
                -floatpressureeast / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                -floatpressurenorth / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                latitude,
                dragcoefficient,
                tuning::climate::atmosphere::surfaceBoundaryLayerMomentumDepthMetres,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const auto roundedwind = climateatmosphere::steadyQuadraticDragCoriolisWind(
                -roundedpressureeast / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                -roundedpressurenorth / tuning::climate::atmosphere::surfaceAirDensityKgM3,
                latitude,
                dragcoefficient,
                tuning::climate::atmosphere::surfaceBoundaryLayerMomentumDepthMetres,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const float floatheighteast =
                (upperheight[xeast][y] - upperheight[xwest][y]) /
                (2.0f * spacing.zonalMetres);
            const float floatheightnorth =
                (upperheight[x][ynorth] - upperheight[x][ysouth]) /
                (2.0f * spacing.meridionalMetres);
            const float roundedheighteast =
                (std::round(upperheight[xeast][y]) - std::round(upperheight[xwest][y])) /
                (2.0f * spacing.zonalMetres);
            const float roundedheightnorth =
                (std::round(upperheight[x][ynorth]) - std::round(upperheight[x][ysouth])) /
                (2.0f * spacing.meridionalMetres);
            const auto floatupperwind = climateatmosphere::steadyRayleighCoriolisWind(
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * floatheighteast,
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * floatheightnorth,
                latitude,
                tuning::climate::circulation::upperLayerDragTimeSeconds,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const auto roundedupperwind = climateatmosphere::steadyRayleighCoriolisWind(
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * roundedheighteast,
                -tuning::climate::atmosphere::gravityMetresPerSecondSquared * roundedheightnorth,
                latitude,
                tuning::climate::circulation::upperLayerDragTimeSeconds,
                tuning::climate::atmosphere::rotationRatePerSecond,
                world.rotation() ? 1.0f : -1.0f);
            const double pressureeastdifference = roundedpressureeast - floatpressureeast;
            const double pressurenorthdifference = roundedpressurenorth - floatpressurenorth;
            const double windeastdifference = roundedwind.eastMetresPerSecond - floatwind.eastMetresPerSecond;
            const double windsouthdifference = roundedwind.southMetresPerSecond - floatwind.southMetresPerSecond;
            const double heighteastdifference = roundedheighteast - floatheighteast;
            const double heightnorthdifference = roundedheightnorth - floatheightnorth;
            const double upperwindeastdifference = roundedupperwind.eastMetresPerSecond - floatupperwind.eastMetresPerSecond;
            const double upperwindsouthdifference = roundedupperwind.southMetresPerSecond - floatupperwind.southMetresPerSecond;

            weighttotal += weight;
            floatpressuregradientsquared += weight *
                (floatpressureeast * floatpressureeast + floatpressurenorth * floatpressurenorth);
            roundedpressuregradientsquared += weight *
                (roundedpressureeast * roundedpressureeast + roundedpressurenorth * roundedpressurenorth);
            pressuregradientdifferencesquared += weight *
                (pressureeastdifference * pressureeastdifference + pressurenorthdifference * pressurenorthdifference);
            floatwindsquared += weight *
                (floatwind.eastMetresPerSecond * floatwind.eastMetresPerSecond +
                    floatwind.southMetresPerSecond * floatwind.southMetresPerSecond);
            roundedwindsquared += weight *
                (roundedwind.eastMetresPerSecond * roundedwind.eastMetresPerSecond +
                    roundedwind.southMetresPerSecond * roundedwind.southMetresPerSecond);
            winddifferencesquared += weight *
                (windeastdifference * windeastdifference + windsouthdifference * windsouthdifference);
            floatheightgradientsquared += weight *
                (floatheighteast * floatheighteast + floatheightnorth * floatheightnorth);
            roundedheightgradientsquared += weight *
                (roundedheighteast * roundedheighteast + roundedheightnorth * roundedheightnorth);
            heightgradientdifferencesquared += weight *
                (heighteastdifference * heighteastdifference + heightnorthdifference * heightnorthdifference);
            floatupperwindsquared += weight *
                (floatupperwind.eastMetresPerSecond * floatupperwind.eastMetresPerSecond +
                    floatupperwind.southMetresPerSecond * floatupperwind.southMetresPerSecond);
            roundedupperwindsquared += weight *
                (roundedupperwind.eastMetresPerSecond * roundedupperwind.eastMetresPerSecond +
                    roundedupperwind.southMetresPerSecond * roundedupperwind.southMetresPerSecond);
            upperwinddifferencesquared += weight *
                (upperwindeastdifference * upperwindeastdifference +
                    upperwindsouthdifference * upperwindsouthdifference);
        }
    }

    const double inverseweight = weighttotal > 0.0 ? 1.0 / weighttotal : 0.0;
    climateatmosphere::CirculationPrecisionStageDiagnostics diagnostics;
    diagnostics.areaWeightedFloatPressureGradientRmsPaPerMetre =
        std::sqrt(floatpressuregradientsquared * inverseweight);
    diagnostics.areaWeightedRoundedPressureGradientRmsPaPerMetre =
        std::sqrt(roundedpressuregradientsquared * inverseweight);
    diagnostics.areaWeightedPressureGradientDifferenceRmsPaPerMetre =
        std::sqrt(pressuregradientdifferencesquared * inverseweight);
    diagnostics.areaWeightedFloatSurfaceWindRmsMetresPerSecond =
        std::sqrt(floatwindsquared * inverseweight);
    diagnostics.areaWeightedRoundedSurfaceWindRmsMetresPerSecond =
        std::sqrt(roundedwindsquared * inverseweight);
    diagnostics.areaWeightedSurfaceWindDifferenceRmsMetresPerSecond =
        std::sqrt(winddifferencesquared * inverseweight);
    diagnostics.areaWeightedFloatUpperHeightGradientRms =
        std::sqrt(floatheightgradientsquared * inverseweight);
    diagnostics.areaWeightedRoundedUpperHeightGradientRms =
        std::sqrt(roundedheightgradientsquared * inverseweight);
    diagnostics.areaWeightedUpperHeightGradientDifferenceRms =
        std::sqrt(heightgradientdifferencesquared * inverseweight);
    diagnostics.areaWeightedFloatUpperWindRmsMetresPerSecond =
        std::sqrt(floatupperwindsquared * inverseweight);
    diagnostics.areaWeightedRoundedUpperWindRmsMetresPerSecond =
        std::sqrt(roundedupperwindsquared * inverseweight);
    diagnostics.areaWeightedUpperWindDifferenceRmsMetresPerSecond =
        std::sqrt(upperwinddifferencesquared * inverseweight);
    return diagnostics;
}
}

void createpressuremap(planet& world)
{
    climatestorage = {};
    processfieldworld = nullptr;
    processfields = {};
    oceanfieldworld = nullptr;
    oceansstfields = {};
    oceanskinfields = {};
    oceanicefields = {};
    oceanaccepted.fill(false);
    atmosphereaccepted.fill(false);
    climatephysics::clearClimateCouplingDiagnostics();
    const int width = world.width();
    const int height = world.height();
    circulationcache.reset(world);
    const float thicknessresponse = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
        tuning::climate::circulation::surfaceReferencePressurePa,
        tuning::climate::circulation::upperReferencePressurePa);

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        floatgrid surface(width + 1, vector<float>(height + 1, 0.0f));

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                    surface[x][y] = pressuresurfacetemperature(world, season, x, y);
            }
        });

        double temperaturesum = 0.0;
        double weighttotal = 0.0;

        for (int y = 0; y <= height; y++)
        {
            constexpr double pi = 3.14159265358979323846;
            const double weight = std::max(
                0.0,
                std::cos(static_cast<double>(latitudeforrow(y, height)) * pi / 180.0));

            for (int x = 0; x <= width; x++)
            {
                temperaturesum += weight * static_cast<double>(surface[x][y]);
                weighttotal += weight;
            }
        }

        const float globalmean = weighttotal > 0.0 ?
            static_cast<float>(temperaturesum / weighttotal) : 0.0f;
        floatgrid stationarytemperature = surface;
        smoothallfield(
            world,
            stationarytemperature,
            tuning::climate::circulation::thermalHeightSmoothingIterations);
        floatgrid pressure(width + 1, vector<float>(height + 1, 0.0f));
        floatgrid upperheight(width + 1, vector<float>(height + 1, 0.0f));
        vector<float> zonaltemperature(height + 1, 0.0f);
        vector<float> stationaryzonaltemperature(height + 1, 0.0f);

        for (int y = 0; y <= height; y++)
        {
            double rowsum = 0.0;
            double stationaryrowsum = 0.0;

            for (int x = 0; x <= width; x++)
            {
                rowsum += static_cast<double>(surface[x][y]);
                stationaryrowsum += static_cast<double>(stationarytemperature[x][y]);
            }

            zonaltemperature[y] = static_cast<float>(rowsum / static_cast<double>(width + 1));
            stationaryzonaltemperature[y] = static_cast<float>(
                stationaryrowsum / static_cast<double>(width + 1));
        }

        for (int iteration = 0;
             iteration < tuning::climate::circulation::zonalTemperatureSmoothingIterations;
             iteration++)
        {
            vector<float> smoothed = zonaltemperature;

            for (int y = 0; y <= height; y++)
            {
                float total = 2.0f * zonaltemperature[y];
                float weights = 2.0f;

                if (y > 0)
                {
                    total += zonaltemperature[y - 1];
                    weights += 1.0f;
                }

                if (y < height)
                {
                    total += zonaltemperature[y + 1];
                    weights += 1.0f;
                }

                smoothed[y] = total / weights;
            }

            zonaltemperature.swap(smoothed);
        }

        const int thermalrow = static_cast<int>(std::distance(
            zonaltemperature.begin(),
            std::max_element(zonaltemperature.begin(), zonaltemperature.end())));
        const int polarcaprows = std::max(1, height / 12);
        double polartemperaturetotal = 0.0;
        int polartemperaturecount = 0;

        for (int y = 0; y <= polarcaprows; y++)
        {
            polartemperaturetotal += zonaltemperature[y];
            polartemperaturetotal += zonaltemperature[height - y];
            polartemperaturecount += 2;
        }

        const float polartemperature = polartemperaturecount > 0 ?
            static_cast<float>(polartemperaturetotal / polartemperaturecount) : globalmean;
        const float temperaturecontrast = std::max(
            1.0f,
            zonaltemperature[thermalrow] - polartemperature);
        const float gravity = tuning::climate::atmosphere::gravityMetresPerSecondSquared *
            std::max(0.05f, world.gravity());
        const float hadleyhalfwidth = std::clamp(
            climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
                temperaturecontrast,
                tuning::climate::circulation::troposphereHeightMetres,
                tuning::climate::circulation::referenceTemperatureK,
                gravity,
                tuning::climate::atmosphere::rotationRatePerSecond,
                tuning::climate::moistureadvection::referencePlanetRadiusMetres),
            tuning::climate::circulation::minimumHadleyHalfWidthDegrees,
            tuning::climate::circulation::maximumHadleyHalfWidthDegrees);
        const float maximumthermalequatorshift = hadleyhalfwidth *
            tuning::climate::circulation::maximumThermalEquatorShiftHadleyFraction;
        const float thermalequator = std::clamp(
            latitudeforrow(thermalrow, height),
            -maximumthermalequatorshift,
            maximumthermalequatorshift);
        const float pressureamplitude = std::clamp(
            tuning::climate::circulation::surfaceReferencePressurePa / 100.0f *
                temperaturecontrast / tuning::climate::circulation::referenceTemperatureK *
                tuning::climate::circulation::overturningMassRedistributionEfficiency,
            tuning::climate::circulation::minimumOverturningPressureAmplitudeHpa,
            tuning::climate::circulation::maximumOverturningPressureAmplitudeHpa);
        vector<float> zonalpressure(height + 1, 0.0f);
        double pressuretotal = 0.0;
        double pressureweighttotal = 0.0;

        for (int y = 0; y <= height; y++)
        {
            constexpr double pi = 3.14159265358979323846;
            const float latitude = latitudeforrow(y, height);
            const double areaweight = std::max(
                0.0,
                std::cos(static_cast<double>(latitude) * pi / 180.0));
            zonalpressure[y] = climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(
                latitude,
                thermalequator,
                hadleyhalfwidth,
                pressureamplitude);
            pressuretotal += areaweight * zonalpressure[y];
            pressureweighttotal += areaweight;
        }

        const float meanpressure = pressureweighttotal > 0.0 ?
            static_cast<float>(pressuretotal / pressureweighttotal) : 0.0f;

        for (float& rowpressure : zonalpressure)
            rowpressure -= meanpressure;

        vector<float> localstationarytemperature(
            static_cast<size_t>(width + 1) * static_cast<size_t>(height + 1),
            0.0f);
        for (int y = 0; y <= height; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                localstationarytemperature[
                    static_cast<size_t>(y) * (width + 1) + x] =
                    stationarytemperature[x][y] - stationaryzonaltemperature[y];
            }
        }
        cout
            << "Overturning pressure season " << season
            << ": thermal_equator=" << thermalequator
            << " C_contrast=" << temperaturecontrast
            << " hadley_half_width=" << hadleyhalfwidth
            << " pressure_amplitude_hpa=" << pressureamplitude << '\n';

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    pressure[x][y] = zonalpressure[y];
                    upperheight[x][y] = 0.0f;
                }
            }
        });

        smoothallfield(world, upperheight, tuning::climate::circulation::thermalHeightSmoothingIterations);

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    world.setseasonalpressure(season, x, y, static_cast<int>(std::round(pressure[x][y])));
                    world.setseasonalupperheight(season, x, y, static_cast<int>(std::round(upperheight[x][y])));
                }
            }
        });

        circulationcache.pressure[season] = std::move(pressure);
        circulationcache.upperheight[season] = std::move(upperheight);
        circulationcache.populated[season] = true;
    }
}

void updatehorsebeltsfrompressure(planet& world)
{
    const int width = world.width();
    const int height = world.height();

    auto rowfromlatitude = [height](float latitude)
    {
        return std::clamp(static_cast<int>(std::round((90.0f - latitude) * static_cast<float>(height) / 180.0f)), 0, height);
    };

    const int defaultnorthouter = rowfromlatitude(32.0f);
    const int defaultnorthinner = rowfromlatitude(26.0f);
    const int defaultsouthinner = rowfromlatitude(-26.0f);
    const int defaultsouthouter = rowfromlatitude(-32.0f);
    const int minband = std::max(2, height / 60);

    vector<int> northouter(width + 1, defaultnorthouter);
    vector<int> northinner(width + 1, defaultnorthinner);
    vector<int> southinner(width + 1, defaultsouthinner);
    vector<int> southouter(width + 1, defaultsouthouter);

    for (int x = 0; x <= width; x++)
    {
        int northpeakpoleward = height;
        int northpeakequatorward = 0;
        int southpeakequatorward = height;
        int southpeakpoleward = 0;
        bool foundnorth = false;
        bool foundsouth = false;

        for (int season = 0; season < CLIMATESEASONCOUNT; season++)
        {
            float bestnorth = -1000000.0f;
            float bestsouth = -1000000.0f;
            int bestnorthy = defaultnorthinner;
            int bestsouthy = defaultsouthinner;

            for (int y = 0; y <= height; y++)
            {
                const float latitude = latitudeforrow(y, height);
                float pressure = 0.0f;

                for (int dx = -1; dx <= 1; dx++)
                {
                    const int xx = wrapx(x + dx, width);
                    pressure += static_cast<float>(world.seasonalpressure(season, xx, y));
                }

                pressure = pressure / 3.0f;

                if (latitude >= 15.0f && latitude <= 45.0f)
                {
                    if (pressure > bestnorth)
                    {
                        bestnorth = pressure;
                        bestnorthy = y;
                    }
                }

                if (latitude <= -15.0f && latitude >= -45.0f)
                {
                    if (pressure > bestsouth)
                    {
                        bestsouth = pressure;
                        bestsouthy = y;
                    }
                }
            }

            northpeakpoleward = std::min(northpeakpoleward, bestnorthy);
            northpeakequatorward = std::max(northpeakequatorward, bestnorthy);
            southpeakequatorward = std::min(southpeakequatorward, bestsouthy);
            southpeakpoleward = std::max(southpeakpoleward, bestsouthy);
            foundnorth = true;
            foundsouth = true;
        }

        if (foundnorth)
        {
            if (northpeakequatorward - northpeakpoleward < minband)
            {
                const int centre = (northpeakequatorward + northpeakpoleward) / 2;
                northpeakpoleward = std::max(rowfromlatitude(45.0f), centre - minband / 2);
                northpeakequatorward = std::min(rowfromlatitude(15.0f), northpeakpoleward + minband);
            }

            northouter[x] = northpeakpoleward;
            northinner[x] = northpeakequatorward;
        }

        if (foundsouth)
        {
            if (southpeakpoleward - southpeakequatorward < minband)
            {
                const int centre = (southpeakpoleward + southpeakequatorward) / 2;
                southpeakequatorward = std::max(rowfromlatitude(-15.0f), centre - minband / 2);
                southpeakpoleward = std::min(rowfromlatitude(-45.0f), southpeakequatorward + minband);
            }

            southinner[x] = southpeakequatorward;
            southouter[x] = southpeakpoleward;
        }
    }

    auto smoothband = [width](vector<int>& band)
    {
        vector<int> smoothed = band;

        for (int x = 0; x <= width; x++)
        {
            int total = 0;
            int count = 0;

            for (int dx = -2; dx <= 2; dx++)
            {
                total += band[wrapx(x + dx, width)];
                count++;
            }

            smoothed[x] = total / count;
        }

        band.swap(smoothed);
    };

    smoothband(northouter);
    smoothband(northinner);
    smoothband(southinner);
    smoothband(southouter);

    for (int x = 0; x <= width; x++)
    {
        if (northinner[x] <= northouter[x])
            northinner[x] = std::min(height, northouter[x] + minband);

        if (southouter[x] <= southinner[x])
            southouter[x] = std::min(height, southinner[x] + minband);

        world.sethorse(x, 1, northouter[x]);
        world.sethorse(x, 2, northinner[x]);
        world.sethorse(x, 3, southinner[x]);
        world.sethorse(x, 4, southouter[x]);
    }
}

void convergeclimatecoupling(planet& world, vector<vector<int>>& fractal)
{
    using namespace tuning::climate::circulation;
    if (processfieldworld != &world) return;
    const auto snapshot = [&]()
    {
        std::array<vector<double>, 4> values;
        for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
        {
            const auto& fields = processfields[season];
            for (int y = 0; y < fields.rows; ++y)
                for (int x = 0; x < fields.columns; ++x)
                {
                    const int cell = y * fields.columns + x;
                    const int xx = std::min(world.width(), static_cast<int>((x + 0.5) * (world.width() + 1) / fields.columns));
                    const int yy = std::clamp(static_cast<int>(std::round((y + 0.5) * world.height() / fields.rows)), 0, world.height());
                    const double weight = std::sqrt(std::max(0.0, std::cos(3.141592653589793 * ((y + 0.5) / fields.rows - 0.5))));
                    values[0].push_back(weight * hydrologyforcingcache.surfacewindu[season][xx][yy]);
                    values[0].push_back(weight * hydrologyforcingcache.surfacewindv[season][xx][yy]);
                    values[0].push_back(weight * hydrologyforcingcache.upperwindu[season][xx][yy]);
                    values[0].push_back(weight * hydrologyforcingcache.upperwindv[season][xx][yy]);
                    if (world.sea(xx, yy)) values[1].push_back(weight * pressuresurfacetemperature(world, season, xx, yy));
                    for (int layer = 0; layer < 2; ++layer)
                        values[2].push_back(weight * (fields.radiativeHeatingWm2[layer][cell] + fields.latentHeatingWm2[layer][cell] +
                            (layer == 0 ? fields.sensibleHeatingWm2[cell] : 0.0f)));
                    values[2].push_back(weight * fields.surfaceNetHeatingWm2[cell]);
                    values[3].push_back(weight * world.seasonalrainfloat(season, xx, yy));
                }
        }
        return values;
    };
    const auto relative = [](const vector<double>& old, const vector<double>& now)
    {
        if (old.size() != now.size()) return std::numeric_limits<double>::infinity();
        double change = 0.0, scale = 0.0;
        for (std::size_t i = 0; i < old.size(); ++i)
        {
            if (!std::isfinite(old[i]) || !std::isfinite(now[i])) return std::numeric_limits<double>::infinity();
            change += std::pow(now[i] - old[i], 2);
            scale += std::max(old[i] * old[i], now[i] * now[i]);
        }
        return std::sqrt(change / std::max(static_cast<double>(std::max<std::size_t>(1, old.size())), scale));
    };
    auto old = snapshot();
    auto previousforcing = processfields;
    vector<double> previousinput, previousresidual;
    double relaxation = couplingHeatingRelaxation;
    for (int iteration = 1; iteration <= maximumCouplingIterations; ++iteration)
    {
        if (iteration > 1)
        {
            auto residual = snapshot()[2];
            for (std::size_t i = 0; i < residual.size(); ++i) residual[i] -= previousinput[i];
            relaxation = climatephysics::adaptiveCouplingRelaxation(relaxation, previousresidual, residual,
                couplingMinimumHeatingRelaxation, couplingHeatingRelaxation);
            previousresidual = std::move(residual);
            for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
            {
                auto& fields = processfields[season];
                const auto blend = [&](vector<float>& now, const vector<float>& previous) {
                    for (std::size_t i = 0; i < now.size(); ++i)
                        now[i] = previous[i] + static_cast<float>(relaxation) * (now[i] - previous[i]); };
                for (int layer = 0; layer < 2; ++layer)
                {
                    blend(fields.radiativeHeatingWm2[layer], previousforcing[season].radiativeHeatingWm2[layer]);
                    blend(fields.latentHeatingWm2[layer], previousforcing[season].latentHeatingWm2[layer]);
                }
                blend(fields.sensibleHeatingWm2, previousforcing[season].sensibleHeatingWm2);
                blend(fields.surfaceNetHeatingWm2, previousforcing[season].surfaceNetHeatingWm2);
            }
        }
        previousforcing = processfields;
        // Test the generated heating against the forcing actually supplied to
        // this iteration, not against the previous unrelaxed diagnostic.
        old[2] = snapshot()[2];
        previousinput = old[2];
        createvectorwindmap(world);
        updatehorsebeltsfrompressure(world);
        createoceancurrentmap(world);
        createsurfacetemperaturemap(world);
        refreshadvectedrainfall(world, fractal);
        const auto now = snapshot();
        climatephysics::ClimateCouplingDiagnostics d;
        d.iteration = iteration;
        d.heatingRelaxation = iteration == 1 ? 1.0 : relaxation;
        d.windChange = relative(old[0], now[0]); d.sstChange = relative(old[1], now[1]);
        d.heatingChange = relative(old[2], now[2]); d.rainfallChange = relative(old[3], now[3]);
        std::array<double, 3> heatingChange{}, heatingScale{};
        std::size_t peakHeatingIndex = 0;
        double peakHeatingChange = 0.0;
        for (std::size_t i = 0; i < old[2].size(); ++i)
        {
            const double change = now[2][i] - old[2][i];
            heatingChange[i % 3] += change * change;
            heatingScale[i % 3] += std::max(old[2][i] * old[2][i], now[2][i] * now[2][i]);
            if (std::abs(change) > std::abs(peakHeatingChange))
            {
                peakHeatingChange = change;
                peakHeatingIndex = i;
            }
        }
        const int forcingColumns = std::max(1, processfields[0].columns);
        const int forcingRows = std::max(1, processfields[0].rows);
        const std::size_t peakCell = peakHeatingIndex / 3;
        d.innerSolvesAccepted = std::all_of(oceanaccepted.begin(), oceanaccepted.end(), [](bool accepted) { return accepted; }) &&
            std::all_of(atmosphereaccepted.begin(), atmosphereaccepted.end(), [](bool accepted) { return accepted; }) &&
            climatephysics::lastHydrologySpinupDiagnostics().converged;
        d.converged = climatephysics::climateCouplingConverged(d, minimumCouplingIterations, couplingRelativeTolerance);
        climatephysics::appendClimateCouplingDiagnostics(d);
        std::cout << "Climate coupling iteration=" << iteration << " wind_change=" << d.windChange << " sst_change=" << d.sstChange
            << " heating_change=" << d.heatingChange << " rainfall_change=" << d.rainfallChange
            << " heating_relaxation=" << d.heatingRelaxation
            << " heating_lower=" << std::sqrt(heatingChange[0] / std::max(1.0, heatingScale[0]))
            << " heating_upper=" << std::sqrt(heatingChange[1] / std::max(1.0, heatingScale[1]))
            << " heating_surface=" << std::sqrt(heatingChange[2] / std::max(1.0, heatingScale[2]))
            << " heating_peak_component=" << peakHeatingIndex % 3
            << " heating_peak_season=" << peakCell / (forcingColumns * forcingRows)
            << " heating_peak_lat=" << 90.0 - 180.0 * (peakCell / forcingColumns % forcingRows + 0.5) / forcingRows
            << " heating_peak_lon=" << -180.0 + 360.0 * (peakCell % forcingColumns + 0.5) / forcingColumns
            << " heating_peak_area_scaled_wm2=" << peakHeatingChange
            << " inner_accepted=" << d.innerSolvesAccepted << " converged=" << d.converged << '\n';
        if (d.converged) break;
        old = now;
    }
}

void createvectorwindmap(planet& world)
{
    using namespace tuning::climate;
    const int width = world.width(), height = world.height();
    const int columns = climatehydrology::climateGridDimensions(width + 1, height + 1,
        circulation::stationaryWaveLongitudeCells).columns;
    const int rows = columns / 2;
    if (columns < 6) return;
    const auto grid = climategrid::makeSphericalGrid(columns, rows, atmosphere::referencePlanetRadiusMetres);
    const auto input = climategrid::makeConservativeRemap(width + 1, height + 1,
        climategrid::LatitudeLayout::poleInclusive, columns, rows, climategrid::LatitudeLayout::cellCentred);
    const auto sample = [&](auto accessor)
    {
        vector<float> result(columns * rows, 0.0f);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                double sum = 0.0;
                for (const auto& wy : input.latitude[y])
                    for (const auto& wx : input.longitude[x])
                        sum += wy.fraction * wx.fraction * accessor(wx.source, wy.source);
                result[y * columns + x] = static_cast<float>(sum);
            }
        return result;
    };
    const auto output = [&](const vector<float>& field)
    {
        const auto flat = climategrid::remapField(columns, rows, climategrid::LatitudeLayout::cellCentred,
            field, width + 1, height + 1, climategrid::LatitudeLayout::poleInclusive);
        floatgrid result(width + 1, vector<float>(height + 1));
        for (int y = 0; y <= height; ++y)
            for (int x = 0; x <= width; ++x) result[x][y] = flat[y * (width + 1) + x];
        return result;
    };
    const auto land = sample([&](int x, int y) { return world.nom(x, y) > world.sealevel() ? 1.0f : 0.0f; });
    const auto terrain = sample([&](int x, int y) { return std::max(0, world.nom(x, y) - world.sealevel()); });
    const auto terrainoutput = output(terrain);
    const auto landoutput = output(land);
    hydrologyforcingcache.reset(world);
    for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
    {
        const auto temperature = sample([&](int x, int y) { return pressuresurfacetemperature(world, season, x, y); });
        const auto rain = sample([&](int x, int y) { return std::max(0.0f, world.seasonalrainfloat(season, x, y)); });
        vector<float> rowtemperature(rows, 0.0f), zonalpressure(rows, 0.0f);
        double globaltemperature = 0.0, totalarea = 0.0;
        for (int y = 0; y < rows; ++y)
        {
            for (int x = 0; x < columns; ++x) rowtemperature[y] += temperature[y * columns + x] / columns;
            globaltemperature += grid.cellAreasSquareMetres[y] * rowtemperature[y];
            totalarea += grid.cellAreasSquareMetres[y];
        }
        globaltemperature /= totalarea;
        const int warmrow = static_cast<int>(std::max_element(rowtemperature.begin(), rowtemperature.end()) - rowtemperature.begin());
        const float contrast = std::max(1.0f, rowtemperature[warmrow] - 0.5f * (rowtemperature.front() + rowtemperature.back()));
        climateatmosphere::ModeSeparatedCirculationConfig config;
        config.surfaceBoundaryLayerDepthMetres = atmosphere::surfaceBoundaryLayerMomentumDepthMetres;
        config.gravityMetresPerSecondSquared *= std::max(0.05f, world.gravity());
        config.rotationDirection = world.rotation() ? 1.0f : -1.0f;
        config.enabled = {circulation::enableZonalMode, circulation::enableStationaryMode,
            circulation::enableSurfaceMode, circulation::enableUpperMode};
        config.surfaceToUpperCoupling = circulation::surfaceToUpperModeCoupling;
        config.upperDragTimeSeconds = circulation::upperLayerDragTimeSeconds;
        config.maximumIterations = circulation::stationaryWaveMaximumIterations;
        config.relativeTolerance = circulation::stationaryWaveRelativeTolerance;
        const float hadleywidth = std::clamp(climateatmosphere::heldHouHadleyEdgeLatitudeDegrees(
            contrast, circulation::troposphereHeightMetres, static_cast<float>(globaltemperature + 273.15),
            config.gravityMetresPerSecondSquared, config.rotationRatePerSecond, config.planetRadiusMetres),
            circulation::minimumHadleyHalfWidthDegrees, circulation::maximumHadleyHalfWidthDegrees);
        const float thermalshift = hadleywidth * circulation::maximumThermalEquatorShiftHadleyFraction;
        const float thermalequator = std::clamp(static_cast<float>(grid.latitudeCentresRadians[warmrow] * 180.0 / 3.141592653589793),
            -thermalshift, thermalshift);
        const float amplitude = std::clamp(circulation::surfaceReferencePressurePa / 100.0f *
            contrast / circulation::referenceTemperatureK * circulation::overturningMassRedistributionEfficiency,
            circulation::minimumOverturningPressureAmplitudeHpa, circulation::maximumOverturningPressureAmplitudeHpa);
        config.zonalUpperHeightMetres.resize(rows);
        const float hydrostatic = climateatmosphere::hypsometricHeightResponseMetresPerKelvin(
            circulation::surfaceReferencePressurePa, circulation::upperReferencePressurePa, 287.05f,
            config.gravityMetresPerSecondSquared);
        double pressuremean = 0.0;
        for (int y = 0; y < rows; ++y)
        {
            zonalpressure[y] = climateatmosphere::axisymmetricOverturningPressureAnomalyHpa(
                static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793),
                thermalequator, hadleywidth, amplitude);
            pressuremean += grid.cellAreasSquareMetres[y] * zonalpressure[y] / totalarea;
            config.zonalUpperHeightMetres[y] = hydrostatic * (rowtemperature[y] - static_cast<float>(globaltemperature));
        }
        for (float& value : zonalpressure) value -= static_cast<float>(pressuremean);

        // Stable dry/moist lapse-rate blend is diagnosed from simulated rain.
        double meanrain = 0.0;
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
                meanrain += rain[y * columns + x] * grid.cellAreasSquareMetres[y] / (totalarea * columns);
        const float humidity = std::clamp(static_cast<float>(meanrain) / lapse::rainfallHumidityScale, 0.0f, 1.0f);
        const float lapserate = (lapse::dryLapseRate + humidity * (lapse::moistLapseRate - lapse::dryLapseRate)) / 1000.0f;
        const float stability = std::max(0.001f, climateatmosphere::diagnoseBruntVaisalaFrequency(
            static_cast<float>(globaltemperature + 273.15), lapserate, config.gravityMetresPerSecondSquared));
        const auto surfaceparameters = climateatmosphere::diagnoseStationaryParameters(stability,
            circulation::surfaceModeDepthMetres, config.gravityMetresPerSecondSquared, config.planetRadiusMetres,
            config.rotationRatePerSecond, columns, rows, circulation::stationaryWaveNondimensionalDamping,
            circulation::diagnosedForcingScaleMetres);
        const auto upperparameters = climateatmosphere::diagnoseStationaryParameters(stability,
            circulation::upperModeDepthMetres, config.gravityMetresPerSecondSquared, config.planetRadiusMetres,
            config.rotationRatePerSecond, columns, rows, circulation::stationaryWaveNondimensionalDamping,
            circulation::diagnosedForcingScaleMetres);
        const float pressurepermetre = config.airDensityKgM3 * config.gravityMetresPerSecondSquared / 100.0f;
        config.surfaceEquivalentPressureDepthHpa = pressurepermetre * surfaceparameters.equivalentDepthMetres;
        config.upperEquivalentPressureDepthHpa = pressurepermetre * upperparameters.equivalentDepthMetres;
        config.surfaceDampingTimeSeconds = surfaceparameters.dampingTimeSeconds;
        config.upperDampingTimeSeconds = upperparameters.dampingTimeSeconds;
        config.maximumZonalWavenumber = surfaceparameters.maximumZonalWavenumber;
        config.maximumMeridionalWavenumber = surfaceparameters.maximumMeridionalWavenumber;
        config.upperMaximumZonalWavenumber = upperparameters.maximumZonalWavenumber;
        config.upperMaximumMeridionalWavenumber = upperparameters.maximumMeridionalWavenumber;
        config.surfaceDragTimesSeconds.resize(columns * rows);
        config.surfaceDragCoefficients.resize(columns * rows);
        for (int cell = 0; cell < columns * rows; ++cell)
        {
            const float basedrag = atmosphere::oceanMomentumDragCoefficient + land[cell] *
                (atmosphere::landMomentumDragCoefficient - atmosphere::oceanMomentumDragCoefficient);
            const float relief = std::clamp(terrain[cell] / 4500.0f, 0.0f, 1.0f);
            const float drag = basedrag + relief * (atmosphere::highReliefMomentumDragCoefficient - basedrag);
            config.surfaceDragCoefficients[cell] = drag;
            config.surfaceDragTimesSeconds[cell] = atmosphere::surfaceBoundaryLayerMomentumDepthMetres /
                (drag * circulation::stationaryWaveLinearizationWindMps);
        }
        vector<float> zero(columns * rows, 0.0f);
        auto backgroundconfig = config;
        backgroundconfig.enabled = {true, false, true, true};
        const auto zonal = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonalpressure, zero, zero, backgroundconfig);
        auto surfaceforcing = climateatmosphere::mechanicalTopographicPressureForcingHpa(columns, rows,
            terrain, zonal.surfaceEastWindMps, zonal.surfaceSouthWindMps, 1.0f,
            circulation::mechanicalTopographicTerrainScaleMetres, circulation::mechanicalTopographicPressureAmplitudeHpa,
            circulation::mechanicalTopographicMinimumWindMps, circulation::mechanicalTopographicFullStrengthWindMps,
            circulation::mechanicalTopographicMinimumLatitudeDegrees, circulation::mechanicalTopographicFullStrengthLatitudeDegrees);
        // Hydrostatic boundary-layer pressure response to the simulated surface
        // temperature, with anomalies tapering linearly to zero aloft. Net
        // column heating alone omits this temperature-gradient-driven flow
        // (Lindzen & Nigam, 1987). Keep the existing diabatic mode as well.
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const int cell = y * columns + x;
                surfaceforcing[cell] += climateatmosphere::thermalModePressureAnomalyHpa(
                    0.5f * (temperature[cell] - rowtemperature[y]), config.airDensityKgM3,
                    circulation::surfaceReferencePressurePa, circulation::surfaceThermalTopPressurePa);
            }
        config.upperOrographicHeightMetres = climateatmosphere::upperOrographicHeightForcing(columns, rows,
            terrain, zonal.upperEastWindMps, stability, hydrostatic * static_cast<float>(globaltemperature + 273.15),
            config.upperDampingTimeSeconds, config);
        constexpr std::array<float, CLIMATESEASONCOUNT> days = {15.0f, 105.0f, 196.0f, 288.0f};
        const float declination = climateenergy::solarDeclinationRadians(days[season], world.tilt());
        const float distance = climateenergy::orbitalDistanceFactor(days[season], world.eccentricity(), world.perihelion());
        vector<float> lowerheating(columns * rows), upperheating(columns * rows);
        for (int y = 0; y < rows; ++y)
            for (int x = 0; x < columns; ++x)
            {
                const int cell = y * columns + x;
                const float albedo = static_cast<float>(temperature[cell] < -2.0f ? energybalance::snowAlbedo :
                    energybalance::oceanAlbedo + land[cell] * (energybalance::landAlbedo - energybalance::oceanAlbedo));
                climateatmosphere::ColumnHeatingInput initial;
                initial.incomingSolarWm2 = climateenergy::dailyMeanInsolationWm2(
                    static_cast<float>(grid.latitudeCentresRadians[y] * 180.0 / 3.141592653589793), declination, distance);
                initial.surfaceAlbedo = albedo;
                initial.surfaceTemperatureK = temperature[cell] + 273.15;
                initial.airTemperatureK = {initial.surfaceTemperatureK,
                    initial.surfaceTemperatureK - moistureadvection::freeTroposphereEnvironmentalLapseC};
                initial.longwaveOpticalDepth = {circulation::lowerLongwaveOpticalDepth, circulation::upperLongwaveOpticalDepth};
                initial.shortwaveOpticalDepth = {circulation::lowerShortwaveOpticalDepth, circulation::upperShortwaveOpticalDepth};
                const auto heat = climateatmosphere::diagnoseColumnHeating(initial);
                lowerheating[cell] = static_cast<float>(heat.totalWm2[0]);
                upperheating[cell] = static_cast<float>(heat.totalWm2[1]);
            }
        if (processfieldworld == &world && processfields[season].durationSeconds > 0.0)
        {
            const auto& fields = processfields[season];
            for (int layer = 0; layer < 2; ++layer)
            {
                vector<float> net(fields.columns * fields.rows);
                for (std::size_t cell = 0; cell < net.size(); ++cell)
                    net[cell] = fields.radiativeHeatingWm2[layer][cell] +
                        circulation::stationaryLatentHeatProjectionFraction * fields.latentHeatingWm2[layer][cell] +
                        (layer == 0 ? fields.sensibleHeatingWm2[cell] : 0.0f);
                (layer == 0 ? lowerheating : upperheating) = climategrid::remapField(fields.columns, fields.rows,
                    climategrid::LatitudeLayout::cellCentred, net, columns, rows, climategrid::LatitudeLayout::cellCentred);
            }
        }
        // Q / (cp * dp/g) is a temperature tendency. Damping time and the
        // hypsometric response convert it to each mode's equilibrium forcing.
        const float columnmass = circulation::surfaceReferencePressurePa / config.gravityMetresPerSecondSquared;
        config.surfaceHeatingPressureResponseHpaPerWm2 = -pressurepermetre * hydrostatic * config.surfaceDampingTimeSeconds /
            (1004.0f * columnmass * moistureadvection::boundaryLayerCapacityFraction);
        config.upperHeatingHeightResponseMetresPerWm2 = hydrostatic * config.upperDampingTimeSeconds /
            (1004.0f * columnmass * moistureadvection::freeTroposphereCapacityFraction);
        config.upperStationaryHeatingWm2 = upperheating;
        auto result = climateatmosphere::solveModeSeparatedCirculation(columns, rows, zonalpressure,
            lowerheating, surfaceforcing, config);
        atmosphereaccepted[season] = result.surfaceStationarySolver.converged && result.upperStationarySolver.converged;
        if (circulation::enablePrognosticClimateAtmosphere && config.enabled.stationary &&
            result.surfaceStationarySolver.converged && result.upperStationarySolver.converged)
        {
            climateweather::ShallowWaterConfig dynamics;
            dynamics.layerCount = 2;
            dynamics.rotationDirection = config.rotationDirection;
            dynamics.gravityMetresPerSecondSquared = config.gravityMetresPerSecondSquared;
            dynamics.lowerMeanDepthMetres = surfaceparameters.equivalentDepthMetres;
            dynamics.upperMeanDepthMetres = upperparameters.equivalentDepthMetres;
            dynamics.lowerDragTimeSeconds = config.surfaceDragTimeSeconds;
            dynamics.upperDragTimeSeconds = config.upperDragTimeSeconds;
            dynamics.heightRelaxationTimeSeconds = config.surfaceDampingTimeSeconds;
            auto state = climateweather::makeState(columns, rows, 2, static_cast<std::uint64_t>(world.seed()) + season + 1);
            climateweather::ShallowWaterForcing forcing;
            forcing.equilibriumHeightMetres = {zero, zero};
            forcing.backgroundEastWindMps = {zonal.surfaceEastWindMps, zonal.upperEastWindMps};
            forcing.backgroundSouthWindMps = {zonal.surfaceSouthWindMps, zonal.upperSouthWindMps};
            for (int cell = 0; cell < columns * rows; ++cell)
            {
                forcing.equilibriumHeightMetres[0][cell] = result.surfaceStationarySolver.equilibriumPressureAnomalyHpa[cell] / pressurepermetre;
                forcing.equilibriumHeightMetres[1][cell] = result.upperStationarySolver.equilibriumPressureAnomalyHpa[cell] / pressurepermetre;
            }
            bool stable = true;
            for (int step = 0; step < circulation::prognosticClimateSpinupSteps && stable; ++step)
            {
                const auto diagnostics = climateweather::advance(state, dynamics, forcing, circulation::prognosticClimateStepSeconds);
                stable = diagnostics.finite && diagnostics.bounded;
            }
            const auto samples = stable ? climateweather::generateWeatherSequence(state, dynamics, forcing,
                circulation::prognosticClimateSamples, circulation::prognosticClimateStepSeconds) : vector<climateweather::ShallowWaterState>{};
            if (samples.size() == static_cast<std::size_t>(circulation::prognosticClimateSamples))
            {
                const auto surface = climateweather::calculateStatistics(samples, 0);
                const auto upper = climateweather::calculateStatistics(samples, 1);
                for (int y = 0; y < rows; ++y)
                {
                    double surfacemean = 0.0, uppermean = 0.0;
                    for (int x = 0; x < columns; ++x)
                    {
                        surfacemean += surface.meanHeightAnomalyMetres[y * columns + x] / columns;
                        uppermean += upper.meanHeightAnomalyMetres[y * columns + x] / columns;
                    }
                    for (int x = 0; x < columns; ++x)
                    {
                        const int cell = y * columns + x;
                        if (config.enabled.surface)
                            result.surfacePressureAnomalyHpa[cell] += circulation::prognosticClimateBlend *
                                (pressurepermetre * (surface.meanHeightAnomalyMetres[cell] - static_cast<float>(surfacemean)) -
                                    result.surfaceStationarySolver.pressureAnomalyHpa[cell]);
                        if (config.enabled.upper)
                            result.upperHeightAnomalyMetres[cell] += circulation::prognosticClimateBlend *
                                (upper.meanHeightAnomalyMetres[cell] - static_cast<float>(uppermean) -
                                    result.upperStationarySolver.pressureAnomalyHpa[cell] / pressurepermetre);
                    }
                }
                climateatmosphere::diagnoseModeWinds(columns, rows, config, result);
            }
            else
            {
                const auto amplitude = [](const vector<float>& field)
                {
                    float maximum = 0.0f;
                    for (float value : field) maximum = std::max(maximum, std::abs(value));
                    return maximum;
                };
                std::cout << "Prognostic atmosphere fallback: season=" << season
                    << " samples=" << samples.size() << " stable_spinup=" << stable
                    << " surface_forcing_m=" << amplitude(forcing.equilibriumHeightMetres[0])
                    << " upper_forcing_m=" << amplitude(forcing.equilibriumHeightMetres[1])
                    << " background_upper_mps=" << amplitude(zonal.upperEastWindMps) << '\n';
            }
        }
        auto pressure = output(result.surfacePressureAnomalyHpa);
        auto upperheight = output(result.upperHeightAnomalyMetres);
        auto windu = output(result.surfaceEastWindMps), windv = output(result.surfaceSouthWindMps);
        auto upperu = output(result.upperEastWindMps), upperv = output(result.upperSouthWindMps);
        auto vertical = output(result.ascentHpaPerDay);
        // Output remapping is the last dynamics operation. No raster-scale
        // pressure relaxation may erase the accepted zonal or stationary mode.
        storeterrainverticalmotion(world, season, terrainoutput, windu, windv);
        hydrologyforcingcache.surfacewindu[season] = windu;
        hydrologyforcingcache.surfacewindv[season] = windv;
        hydrologyforcingcache.upperwindu[season] = upperu;
        hydrologyforcingcache.upperwindv[season] = upperv;
        hydrologyforcingcache.populated[season] = true;
        climateatmosphere::CirculationPrecisionDiagnostics precision;
        precision.base = measurecirculationprecision(world, pressure, upperheight, landoutput, terrainoutput);
        precision.final = precision.base;
        climateatmosphere::setLastCirculationPrecisionDiagnostics(season, precision);
        if (isworldgendebugrunactive())
            climatevalidation::capturecirculationwindfields(world, season, windu, windv, upperu, upperv);
        for (int y = 0; y <= height; ++y)
            for (int x = 0; x <= width; ++x)
            {
                const auto windstorage = [](float value) { return static_cast<int>(std::round(std::clamp(value,
                    -atmosphere::maxVectorWind, atmosphere::maxVectorWind))); };
                world.setseasonalpressure(season, x, y, static_cast<int>(std::round(pressure[x][y])));
                world.setseasonalupperheight(season, x, y, static_cast<int>(std::round(upperheight[x][y])));
                world.setseasonaluwind(season, x, y, windstorage(windu[x][y]));
                world.setseasonalvwind(season, x, y, windstorage(windv[x][y]));
                world.setseasonalupperuwind(season, x, y, windstorage(upperu[x][y]));
                world.setseasonaluppervwind(season, x, y, windstorage(upperv[x][y]));
                world.setseasonalverticalvelocity(season, x, y, static_cast<int>(std::round(std::clamp(vertical[x][y],
                    -circulation::maximumVerticalVelocity, circulation::maximumVerticalVelocity) * circulation::verticalVelocityStorageScale)));
            }
        std::cout << "Mode circulation season " << season << ": grid=" << columns << 'x' << rows
            << " surface_residual=" << result.surfaceStationarySolver.relativeResidual
            << " upper_residual=" << result.upperStationarySolver.relativeResidual
            << " surface_converged=" << result.surfaceStationarySolver.converged
            << " upper_converged=" << result.upperStationarySolver.converged
            << " stagnated=" << (result.surfaceStationarySolver.stagnated || result.upperStationarySolver.stagnated)
            << " N=" << stability << " depth_m=" << surfaceparameters.equivalentDepthMetres
            << " hadley_deg=" << hadleywidth << " thermal_equator_deg=" << thermalequator
            << " bandwidth=" << config.maximumZonalWavenumber << '/' << config.maximumMeridionalWavenumber
            << " kinetic_jm2=" << result.areaWeightedKineticEnergyJm2
            << " drag_wm2=" << result.areaWeightedDragDissipationWm2
            << " mass_kgm2=" << result.areaWeightedMassAnomalyKgM2
            << " momentum_exchange_residual=" << result.maximumMomentumExchangeResidual
            << " row_transfer_hpa=" << result.maximumStationaryRowMeanHpa << '\n';
    }

    circulationcache.clear();

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                float averageu = 0.0f;
                float averageupperu = 0.0f;

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    averageu += static_cast<float>(world.seasonaluwind(season, x, y));
                    averageupperu += static_cast<float>(world.seasonalupperuwind(season, x, y));
                }

                averageu = averageu / static_cast<float>(CLIMATESEASONCOUNT);
                averageupperu = averageupperu / static_cast<float>(CLIMATESEASONCOUNT);

                int scalarwind = static_cast<int>(std::round(averageu / tuning::climate::atmosphere::scalarWindDivisor));
                scalarwind = std::clamp(scalarwind, -10, 10);

                if (std::fabs(averageu) < tuning::climate::atmosphere::minimumScalarZonalWind || scalarwind == 0)
                {
                    const float transportu = averageu + (averageupperu - averageu) *
                        tuning::climate::moistureadvection::upperWindTransportFraction;

                    if (transportu > 0.25f)
                        world.setwind(x, y, 101);
                    else if (transportu < -0.25f)
                        world.setwind(x, y, 99);
                    else
                        world.setwind(x, y, 0);
                }
                else
                {
                    if (scalarwind > 0)
                        scalarwind = std::max(1, scalarwind);
                    else
                        scalarwind = std::min(-1, scalarwind);

                    world.setwind(x, y, scalarwind);
                }
            }
        }
    });

    for (int y = 0; y <= height; y++)
        world.setwind(0, y, world.wind(width, y));
}

void createadvectedrainfall(planet& world, vector<vector<int>>& inland, vector<vector<int>>& fractal)
{
    (void)inland;
    (void)fractal;
    const int outputwidth = world.width();
    const int outputheight = world.height();
    const climatehydrology::ClimateGridDimensions climatedimensions =
        climatehydrology::climateGridDimensions(
            outputwidth + 1,
            outputheight + 1,
            tuning::climate::moistureadvection::internalClimateHorizontalCells);
    const int width = climatedimensions.columns - 1;
    const int height = climatedimensions.rows - 1;
    const int sealevel = world.sealevel();
    const float gravitymultiplier = std::max(0.05f, world.gravity());
    const bool hasfloatingwinds = hydrologyforcingcache.completeFor(world);
    const int baselineiterations = tuning::climate::moistureadvection::iterations;
    const int referenceiterations = tuning::climateresolution::scaleDistance(
        baselineiterations, width, height);
    const float condensationconversiontimeseconds =
        tuning::climate::moistureadvection::condensationConversionTimeDays *
        tuning::climate::circulation::secondsPerDay;
    vector<float> rowareaweights(height + 1, 0.0f);
    vector<float> climatelatitudes(height + 1, 0.0f);
    vector<float> polartransportfactors(height + 1, 1.0f);
    vector<float> polarconvectionfactors(height + 1, 1.0f);

    for (int y = 0; y <= height; y++)
    {
        climatelatitudes[y] = climatehydrology::climateCellLatitudeDegrees(
            y,
            height + 1);
        rowareaweights[y] = climatehydrology::climateCellAreaWeight(
            y,
            height + 1);
        polartransportfactors[y] = climatehydrology::polarTaperFactor(
            climatelatitudes[y],
            tuning::climate::moistureadvection::polarMeridionalTransportTaperStartDegrees,
            tuning::climate::moistureadvection::polarMeridionalTransportTaperEndDegrees);
        polarconvectionfactors[y] = climatehydrology::polarTaperFactor(
            climatelatitudes[y],
            tuning::climate::moistureadvection::polarConvectionTaperStartDegrees,
            tuning::climate::moistureadvection::polarConvectionTaperEndDegrees);
    }

    const auto inputremap = climategrid::makeConservativeRemap(
        outputwidth + 1, outputheight + 1, climategrid::LatitudeLayout::poleInclusive,
        width + 1, height + 1, climategrid::LatitudeLayout::cellCentred);
    auto aggregateworldfield = [&](int x, int y, auto accessor)
    {
        double total = 0.0;
        for (const auto& yw : inputremap.latitude[y])
        {
            for (const auto& xw : inputremap.longitude[x])
                total += xw.fraction * yw.fraction * accessor(xw.source, yw.source);
        }
        return static_cast<float>(total);
    };
    auto aggregateworldmaximum = [&](int x, int y, auto accessor)
    {
        float maximum = 0.0f;
        for (const auto& yw : inputremap.latitude[y])
        {
            for (const auto& xw : inputremap.longitude[x])
                maximum = std::max(maximum, static_cast<float>(accessor(xw.source, yw.source)));
        }

        return maximum;
    };

    hydrologyfloatgrid landfraction(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid meanelevation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid peakelevation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid icefraction(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid coastnormalu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid coastnormalv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid coastalstrength(width + 1, height + 1, 0.0f);
    hydrologybytegrid seacell(width + 1, height + 1, 0);

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            for (int x = 0; x <= width; x++)
            {
                landfraction[x][y] = aggregateworldfield(
                    x,
                    y,
                    [&](int xx, int yy) { return world.sea(xx, yy) == 0 ? 1.0f : 0.0f; });
                meanelevation[x][y] = aggregateworldfield(
                    x,
                    y,
                    [&](int xx, int yy)
                    {
                        return static_cast<float>(std::max(0, world.map(xx, yy) - sealevel));
                    });
                peakelevation[x][y] = aggregateworldmaximum(
                    x,
                    y,
                    [&](int xx, int yy)
                    {
                        return static_cast<float>(std::max(0, world.map(xx, yy) - sealevel));
                    });
                icefraction[x][y] = aggregateworldfield(
                    x,
                    y,
                    [&](int xx, int yy) { return world.seaice(xx, yy) != 0 ? 1.0f : 0.0f; });
                seacell[x][y] = landfraction[x][y] < 0.5f ? 1 : 0;
            }
        }
    });

    parallelforrows(0, height, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const int north = std::max(0, y - 1);
            const int south = std::min(height, y + 1);

            for (int x = 0; x <= width; x++)
            {
                const float gradienteast = 0.5f *
                    (landfraction[wrapx(x + 1, width)][y] -
                        landfraction[wrapx(x - 1, width)][y]);
                const float gradientsouth = 0.5f *
                    (landfraction[x][south] - landfraction[x][north]);
                const float magnitude = std::sqrt(
                    gradienteast * gradienteast + gradientsouth * gradientsouth);

                if (magnitude > 1.0e-5f)
                {
                    coastnormalu[x][y] = gradienteast / magnitude;
                    coastnormalv[x][y] = gradientsouth / magnitude;
                }

                coastalstrength[x][y] = std::clamp(
                    magnitude * 2.0f + 4.0f * landfraction[x][y] *
                        (1.0f - landfraction[x][y]),
                    0.0f,
                    1.0f);
            }
        }
    });

    std::cout
        << "Climate hydrology grid output=" << outputwidth + 1 << 'x' << outputheight + 1
        << " internal=" << width + 1 << 'x' << height + 1
        << " weather_phases="
        << tuning::climate::moistureadvection::weatherPhaseCount
        << '\n';

    hydrologyfloatgrid soilmoisture(
        width + 1,
        height + 1,
        tuning::climate::moistureadvection::landSoilMoistureCapacity *
            tuning::climate::moistureadvection::initialSoilMoistureFraction);
    hydrologyfloatgrid snowwater(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid boundarymoisture(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid freemoisture(width + 1, height + 1, 0.0f);
    // Continue spin-up across coupled iterations, but never across worlds.
    // Reinitializing these stores discarded the very equilibrium being solved.
    if (climatestorage.source == &world && climatestorage.columns == width + 1 && climatestorage.rows == height + 1)
    {
        boundarymoisture = climatestorage.boundary;
        freemoisture = climatestorage.free;
        soilmoisture = climatestorage.soil;
        snowwater = climatestorage.snow;
    }
    hydrologyfloatgrid rawannualrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid rawmaximummonthlyrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid quarterrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid quarterconvergence(width + 1, height + 1, 0.0f);
    hydrologybytegrid rawmaximumrainmonth(
        width + 1,
        height + 1,
        static_cast<signed char>(-1));
    std::array<hydrologyfloatgrid, CLIMATESEASONCOUNT> storedseasonalrain;
    std::array<hydrologyfloatgrid, CLIMATESEASONCOUNT> storedseasonalmoisture;
    std::array<hydrologyfloatgrid, CLIMATESEASONCOUNT> storedseasonalconvergence;

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        storedseasonalrain[season] = hydrologyfloatgrid(width + 1, height + 1, 0.0f);
        storedseasonalmoisture[season] = hydrologyfloatgrid(width + 1, height + 1, 0.0f);
        storedseasonalconvergence[season] = hydrologyfloatgrid(width + 1, height + 1, 0.0f);
    }
    std::array<climatephysics::WaterBudget, CLIMATESEASONCOUNT> finalareaweightedbudgets{};
    std::array<climatephysics::WaterBudget, CLIMATESEASONCOUNT> finalbudgets{};
    std::array<climatephysics::CondensationActivityDiagnostics, CLIMATESEASONCOUNT>
        finalcondensationactivity{};
    std::array<climatephysics::PrecipitationProcessDiagnostics, CLIMATESEASONCOUNT>
        finalprocessdiagnostics{};

    struct monthlysolverresult
    {
        climatephysics::WaterBudget budget;
        climatephysics::WaterBudget areaweightedbudget;
        climatephysics::CondensationActivityDiagnostics condensation;
        climatephysics::PrecipitationProcessDiagnostics processes;
        std::vector<std::uint8_t> weatherstate;
    };
    std::vector<std::uint8_t> finalweatherstate;
    climateweather::ShallowWaterState continuedweather;
    std::array<vector<climateweather::ShallowWaterState>, CLIMATESEASONCOUNT> seasonalweather;
    std::array<vector<double>, CLIMATESEASONCOUNT> seasonalweights;
    std::array<climatehydrology::SeasonalProcessFields, CLIMATESEASONCOUNT> seasonalprocess;
    hydrologyfloatgrid sensibleheat(width + 1, height + 1, 0.0f);

    hydrologyfloatgrid totalrain(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid totalseaevaporation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid totallandevaporation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid runoff(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid uplift(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid descent(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid dynamicvertical(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid saturationcapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid boundarycapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacepressure(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid temperature(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacewindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacewindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid upperwindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid upperwindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfacewindspeed(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasewindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasewindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseupperwindu(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseupperwindv(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid previousphasewindu, previousphasewindv, previousphaseupperwindu, previousphaseupperwindv;
    hydrologyfloatgrid phasetemperature(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasesurfacetemperature(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasedynamicvertical(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phasesaturationcapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid phaseboundarycapacity(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid persistentheating(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid cloudmemory(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid surfaceevaporation(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid nextboundarymoisture(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid nextfreemoisture(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid boundaryfluxtendency(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid freefluxtendency(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid convectiveconvergence(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid convergencescratch(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid highpasslaplacian(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid highpassscratch(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid totalconvergence(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid zonaleddyexchange(width + 1, height + 1, 0.0f);
    hydrologyfloatgrid meridionaleddyexchange(width + 1, height + 1, 0.0f);

    auto runsolver = [&](int month)
    {
        const climatehydrology::CalendarMonth calendar =
            climatehydrology::calendarMonth(month);
        const int quarter = month / 3;
        auto& exchanges = seasonalprocess[quarter];
        int startday = 0;
        for (int m = 0; m < month; ++m) startday += climatehydrology::calendarMonth(m).days;
        const int iterations = std::max(
            static_cast<int>(std::ceil(calendar.days / tuning::climate::moistureadvection::maximumClimateStepDays)),
            static_cast<int>(std::round(
                static_cast<float>(referenceiterations * calendar.days) /
                static_cast<float>(baselineiterations))));
        const float timestepseconds =
            static_cast<float>(calendar.days) *
            tuning::climate::moistureadvection::advectionTimeStepSeconds /
            static_cast<float>(iterations);
        totalrain.fill(0.0f);
        totalseaevaporation.fill(0.0f);
        totallandevaporation.fill(0.0f);
        runoff.fill(0.0f);
        totalconvergence.fill(0.0f);
        climatephysics::WaterBudget budget;
        climatephysics::WaterBudget areaweightedbudget;
        vector<climatephysics::CondensationActivityDiagnostics> rowcondensation(height + 1);
        vector<climatephysics::PrecipitationProcessDiagnostics> rowprocesses(height + 1);

        if (month % 3 == 0)
        {
            quarterrain.fill(0.0f);
            quarterconvergence.fill(0.0f);
        }

        for (int y = 0; y <= height; y++)
        {
            const double areaweight = rowareaweights[y];

            for (int x = 0; x <= width; x++)
            {
                const float atmosphericwater =
                    boundarymoisture[x][y] + freemoisture[x][y];
                budget.initialAtmosphericStorage += atmosphericwater;
                areaweightedbudget.initialAtmosphericStorage +=
                    areaweight * atmosphericwater;

                if (seacell[x][y] == 0)
                {
                    budget.initialSoilStorage += soilmoisture[x][y];
                    areaweightedbudget.initialSoilStorage += areaweight * soilmoisture[x][y];
                    budget.initialSnowStorage += snowwater[x][y];
                    areaweightedbudget.initialSnowStorage += areaweight * snowwater[x][y];
                }
            }
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    auto aggregatefield = [&](int season, auto accessor)
                    {
                        return aggregateworldfield(
                            x,
                            y,
                            [&](int xx, int yy) { return accessor(season, xx, yy); });
                    };
                    auto interpolatefield = [&](auto accessor)
                    {
                        return climatehydrology::interpolateSeasonal(
                            aggregatefield(calendar.firstSeason, accessor),
                            aggregatefield(calendar.secondSeason, accessor),
                            calendar.interpolation);
                    };
                    auto interpolatecachedwind = [&](
                        const std::array<floatgrid, CLIMATESEASONCOUNT>& field)
                    {
                        return climatehydrology::interpolateSeasonal(
                            aggregateworldfield(
                                x,
                                y,
                                [&](int xx, int yy)
                                {
                                    return field[calendar.firstSeason][xx][yy];
                                }),
                            aggregateworldfield(
                                x,
                                y,
                                [&](int xx, int yy)
                                {
                                    return field[calendar.secondSeason][xx][yy];
                                }),
                            calendar.interpolation);
                    };
                    const float u = hasfloatingwinds
                        ? interpolatecachedwind(hydrologyforcingcache.surfacewindu)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonaluwind(season, xx, yy);
                            });
                    const float v = hasfloatingwinds
                        ? interpolatecachedwind(hydrologyforcingcache.surfacewindv)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonalvwind(season, xx, yy);
                            });
                    const float upperu = hasfloatingwinds
                        ? interpolatecachedwind(hydrologyforcingcache.upperwindu)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonalupperuwind(season, xx, yy);
                            });
                    const float upperv = hasfloatingwinds
                        ? interpolatecachedwind(hydrologyforcingcache.upperwindv)
                        : interpolatefield(
                            [&](int season, int xx, int yy)
                            {
                                return world.seasonaluppervwind(season, xx, yy);
                            });
                    const float upperwinddeparturemultiplier = std::clamp(
                        (tuning::climate::moistureadvection::
                            polarUpperWindBlendTaperEndDegrees -
                            std::fabs(climatelatitudes[y])) /
                            std::max(
                                0.1f,
                                tuning::climate::moistureadvection::
                                    polarUpperWindBlendTaperEndDegrees -
                                tuning::climate::moistureadvection::
                                    polarUpperWindBlendTaperStartDegrees),
                        0.0f,
                        1.0f) *
                        tuning::climate::moistureadvection::
                            freeTroposphereUpperWindDepartureFraction;
                    const bool sea = seacell[x][y] != 0;
                    if (sea && oceanfieldworld == &world && !tuning::climate::oceancurrents::oneWayDiagnosticsOnly)
                        icefraction[x][y] = interpolatefield([&](int season, int xx, int yy)
                            { return oceanicefields[season][xx][yy]; });
                    const float localtemperature = interpolatefield(
                        [&](int season, int xx, int yy)
                        {
                            return pressuresurfacetemperature(world, season, xx, yy);
                        });
                    const float elevation = meanelevation[x][y];
                    const float pressureanomaly = interpolatefield(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalpressure(season, xx, yy);
                        });

                    temperature[x][y] = localtemperature;
                    surfacepressure[x][y] = std::max(
                        100.0f,
                        climatephysics::surfacePressureHpa(
                            elevation,
                            gravitymultiplier) + pressureanomaly);
                    const float boundarytransportwindfraction = 1.0f -
                        landfraction[x][y] *
                            (1.0f - tuning::climate::moistureadvection::
                                landBoundaryLayerTransportWindFraction);
                    surfacewindu[x][y] = u * boundarytransportwindfraction;
                    surfacewindv[x][y] = v * boundarytransportwindfraction;
                    upperwindu[x][y] = u +
                        (upperu - u) * upperwinddeparturemultiplier;
                    upperwindv[x][y] = v +
                        (upperv - v) * upperwinddeparturemultiplier;
                    saturationcapacity[x][y] =
                        climatephysics::saturationColumnWaterAtPressure(
                            localtemperature,
                            surfacepressure[x][y],
                            gravitymultiplier);
                    boundarycapacity[x][y] = saturationcapacity[x][y] *
                        tuning::climate::moistureadvection::boundaryLayerCapacityFraction;
                    surfacewindspeed[x][y] = std::max(
                        tuning::climate::moistureadvection::minimumSurfaceWind,
                        std::sqrt(u * u + v * v));
                    auto interpolateorographicextreme = [&](auto accessor)
                    {
                        auto retainextreme = [&](int season)
                        {
                            const float mean = aggregatefield(season, accessor);
                            const float maximum = aggregateworldmaximum(
                                x,
                                y,
                                [&](int xx, int yy) { return accessor(season, xx, yy); });
                            return mean + (maximum - mean) *
                                tuning::climate::moistureadvection::
                                    subgridOrographicExtremeRetention;
                        };
                        return climatehydrology::interpolateSeasonal(
                            retainextreme(calendar.firstSeason),
                            retainextreme(calendar.secondSeason),
                            calendar.interpolation);
                    };
                    uplift[x][y] = interpolateorographicextreme(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonaluplift(season, xx, yy);
                        }) / tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
                    descent[x][y] = interpolateorographicextreme(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalsubsidence(season, xx, yy);
                        }) / tuning::climate::atmosphere::topographyVerticalMotionStorageScale;
                    dynamicvertical[x][y] = interpolatefield(
                        [&](int season, int xx, int yy)
                        {
                            return world.seasonalverticalvelocity(season, xx, yy);
                        }) /
                        tuning::climate::circulation::verticalVelocityStorageScale;
                }
            }
        });

        createeddyexchangefractions(
            width,
            height,
            climatelatitudes,
            rowareaweights,
            polartransportfactors,
            surfacewindu,
            surfacewindv,
            upperwindu,
            upperwindv,
            timestepseconds * static_cast<float>(
                tuning::climate::moistureadvection::freeTroposphereTransportCadence),
            zonaleddyexchange,
            meridionaleddyexchange);

        auto advectlayer = [&] (
            const hydrologyfloatgrid& source,
            const hydrologyfloatgrid& windu,
            const hydrologyfloatgrid& windv,
            hydrologyfloatgrid& destination,
            float transporttimestepseconds)
        {
            const int columns = width + 1;
            const int rows = height + 1;
            const std::size_t cellcount = static_cast<std::size_t>(columns) * rows;
            std::vector<float> sourcefield(cellcount, 0.0f);
            std::vector<float> zonalwind(cellcount, 0.0f);
            std::vector<float> meridionalwind(cellcount, 0.0f);
            const bool boundary = &source == &boundarymoisture;
            const auto& startu = boundary ? previousphasewindu : previousphaseupperwindu;
            const auto& startv = boundary ? previousphasewindv : previousphaseupperwindv;
            climatehydrology::MpdataOptions options;
            options.endZonalWindMps.resize(cellcount);
            options.endMeridionalWindMps.resize(cellcount);
            for (int y = 0; y < rows; y++)
            {
                for (int x = 0; x < columns; x++)
                {
                    const std::size_t cell = static_cast<std::size_t>(y) * columns + x;
                    sourcefield[cell] = source[x][y];
                    zonalwind[cell] = startu[x][y];
                    meridionalwind[cell] = startv[x][y];
                    options.endZonalWindMps[cell] = windu[x][y];
                    options.endMeridionalWindMps[cell] = windv[x][y];
                }
            }
            options.maximumCourantPerSubstep =
                tuning::climate::moistureadvection::transportMaximumCourant;
            options.correctivePasses =
                tuning::climate::moistureadvection::transportCorrectivePasses;
            options.monotone = true;
            std::vector<float> transported;
            const auto flux = climatehydrology::advectSphericalTracerMpdata(
                columns,
                rows,
                sourcefield,
                zonalwind,
                meridionalwind,
                transporttimestepseconds,
                tuning::climate::moistureadvection::referencePlanetRadiusMetres,
                options,
                transported);
            const int layer = &source == &boundarymoisture ? 0 : 1;
            for (std::size_t cell = 0; cell < cellcount; ++cell)
            {
                exchanges.eastIntegratedFlux[layer][cell] += flux.eastIntegratedFlux[cell];
                exchanges.southIntegratedFlux[layer][cell] += flux.southIntegratedFlux[cell];
            }
            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                        destination[x][y] = transported[
                            static_cast<std::size_t>(y) * columns + x];
                }
            });
        };

        const int weathercolumns = std::min(
            width + 1,
            tuning::climate::moistureadvection::weatherAnomalyHorizontalCells);
        const int weatherrows = std::max(2, weathercolumns / 2);
        climateweather::ShallowWaterConfig weatherconfig;
        weatherconfig.layerCount = 1;
        weatherconfig.rotationDirection = world.rotation() ? 1.0f : -1.0f;
        weatherconfig.gravityMetresPerSecondSquared *= std::max(0.05f, world.gravity());
        weatherconfig.lowerMeanDepthMetres =
            tuning::climate::moistureadvection::weatherAnomalyEquivalentDepthMetres;
        weatherconfig.lowerDragTimeSeconds =
            tuning::climate::moistureadvection::weatherAnomalyDragTimeDays *
            tuning::climate::circulation::secondsPerDay;
        weatherconfig.heightRelaxationTimeSeconds =
            tuning::climate::moistureadvection::weatherAnomalyRelaxationTimeDays *
            tuning::climate::circulation::secondsPerDay;
        weatherconfig.stochasticHeightForcingMetresPerSecond =
            tuning::climate::moistureadvection::weatherAnomalyStochasticHeightMps;
        auto weatherstate = climateweather::makeState(
            weathercolumns,
            weatherrows,
            1,
            static_cast<std::uint64_t>(world.seed()) ^
                (static_cast<std::uint64_t>(month + 1) * 0x9e3779b97f4a7c15ULL));
        if (continuedweather.columns == weathercolumns && continuedweather.rows == weatherrows)
            weatherstate = continuedweather;
        climateweather::ShallowWaterForcing weatherforcing;
        const auto weatherbackground = [&](const hydrologyfloatgrid& field)
        {
            vector<float> flat((width + 1) * (height + 1));
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x) flat[y * (width + 1) + x] = field[x][y];
            return climategrid::remapField(width + 1, height + 1, climategrid::LatitudeLayout::cellCentred,
                flat, weathercolumns, weatherrows, climategrid::LatitudeLayout::cellCentred);
        };
        weatherforcing.backgroundEastWindMps = {weatherbackground(surfacewindu)};
        weatherforcing.backgroundSouthWindMps = {weatherbackground(surfacewindv)};
        bool evolvingweather =
            tuning::climate::moistureadvection::enableEvolvingWeatherAnomalies;
        const auto sampleweather = [&](const std::vector<float>& field, int x, int y)
        {
            const float sourcex = (static_cast<float>(x) + 0.5f) * weathercolumns /
                static_cast<float>(width + 1) - 0.5f;
            const float sourcey = std::clamp((static_cast<float>(y) + 0.5f) * weatherrows /
                static_cast<float>(height + 1) - 0.5f, 0.0f, static_cast<float>(weatherrows - 1));
            const int basex = static_cast<int>(std::floor(sourcex));
            const int x0 = ((basex % weathercolumns) + weathercolumns) % weathercolumns;
            const int x1 = (x0 + 1) % weathercolumns;
            const int y0 = std::clamp(
                static_cast<int>(std::floor(sourcey)), 0, weatherrows - 1);
            const int y1 = std::min(weatherrows - 1, y0 + 1);
            const float fractionx = sourcex - std::floor(sourcex);
            const float fractiony = std::clamp(
                sourcey - std::floor(sourcey), 0.0f, 1.0f);
            const auto at = [&](int xx, int yy)
            {
                return field[static_cast<std::size_t>(yy) * weathercolumns + xx];
            };
            const float north = at(x0, y0) * (1.0f - fractionx) +
                at(x1, y0) * fractionx;
            const float south = at(x0, y1) * (1.0f - fractionx) +
                at(x1, y1) * fractionx;
            return north * (1.0f - fractiony) + south * fractiony;
        };

        for (int iteration = 0; iteration < iterations; iteration++)
        {
            const auto startweather = weatherstate;
            exchanges.durationSeconds += timestepseconds;
            const float dayofyear = startday + (iteration + 0.5f) * calendar.days / iterations;
            const float solardeclination = climateenergy::solarDeclinationRadians(dayofyear, world.tilt());
            const float solardistance = climateenergy::orbitalDistanceFactor(dayofyear, world.eccentricity(), world.perihelion());
            if (evolvingweather)
            {
                const auto weatherdiagnostics = climateweather::advance(
                    weatherstate,
                    weatherconfig,
                    weatherforcing,
                    timestepseconds);
                evolvingweather = weatherdiagnostics.finite && weatherdiagnostics.bounded;
            }
            const auto buildphase = [&](const climateweather::ShallowWaterState& phaseweather, int phaseindex)
            {
            const climatehydrology::WeatherPhase weatherphase =
                climatehydrology::deterministicWeatherPhase(
                    phaseindex,
                    tuning::climate::moistureadvection::weatherPhaseCount,
                    tuning::climate::moistureadvection::enablePrescribedSynopticPerturbations
                        ? tuning::climate::moistureadvection::synopticWindRotationDegrees : 0.0f,
                    tuning::climate::moistureadvection::daytimeLandTemperatureAnomalyC,
                    tuning::climate::moistureadvection::nighttimeLandTemperatureAnomalyC,
                    tuning::climate::moistureadvection::daytimeSeaTemperatureAnomalyC,
                    tuning::climate::moistureadvection::nighttimeSeaTemperatureAnomalyC);
            const float rotationcos = std::cos(weatherphase.windRotationRadians);
            const float rotationsin = std::sin(weatherphase.windRotationRadians);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    constexpr float pi = 3.14159265358979323846f;
                    const float latitude = climatelatitudes[y] * pi / 180.0f;
                    const float latitudetaper = std::max(0.0f, std::cos(latitude));

                    for (int x = 0; x <= width; x++)
                    {
                        const float longitude = 2.0f * pi * static_cast<float>(x) /
                            static_cast<float>(width + 1);
                        const float synopticu = evolvingweather
                            ? std::clamp(
                                sampleweather(phaseweather.layers[0].eastWindMps, x, y),
                                -tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps,
                                tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps)
                            : (tuning::climate::moistureadvection::enablePrescribedSynopticPerturbations ? 1.0f : 0.0f) * tuning::climate::moistureadvection::synopticWindPerturbationMps *
                                std::sin(
                                    2.0f * longitude + weatherphase.synopticPhaseRadians) *
                                latitudetaper;
                        const float synopticv = evolvingweather
                            ? std::clamp(
                                sampleweather(phaseweather.layers[0].southWindMps, x, y),
                                -tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps,
                                tuning::climate::moistureadvection::
                                    maximumWeatherAnomalyWindMps)
                            : (tuning::climate::moistureadvection::enablePrescribedSynopticPerturbations ? 1.0f : 0.0f) * tuning::climate::moistureadvection::synopticWindPerturbationMps *
                                std::cos(
                                    3.0f * longitude - weatherphase.synopticPhaseRadians) *
                                latitudetaper;
                        const float coastalmagnitude = weatherphase.coastalDirection >= 0.0f
                            ? tuning::climate::moistureadvection::daytimeOnshoreCoastalWindMps
                            : tuning::climate::moistureadvection::nighttimeOffshoreCoastalWindMps;
                        const float coastalwind = (tuning::climate::moistureadvection::enablePrescribedCoastalDiurnalCycle ? 1.0f : 0.0f) * weatherphase.coastalDirection *
                            coastalmagnitude * coastalstrength[x][y];
                        const float rotatedu =
                            surfacewindu[x][y] * rotationcos - surfacewindv[x][y] * rotationsin;
                        const float rotatedv =
                            surfacewindu[x][y] * rotationsin + surfacewindv[x][y] * rotationcos;
                        phasewindu[x][y] = rotatedu + synopticu +
                            coastalwind * coastnormalu[x][y];
                        phasewindv[x][y] = (
                            rotatedv + synopticv + coastalwind * coastnormalv[x][y]) *
                            polartransportfactors[y];
                        phaseupperwindu[x][y] =
                            upperwindu[x][y] * rotationcos - upperwindv[x][y] * rotationsin +
                            0.5f * synopticu;
                        phaseupperwindv[x][y] = (
                            upperwindu[x][y] * rotationsin + upperwindv[x][y] * rotationcos +
                            0.5f * synopticv) * polartransportfactors[y];
                        const float diurnalanomaly = tuning::climate::moistureadvection::enablePrescribedCoastalDiurnalCycle
                            ? weatherphase.seaTemperatureAnomalyC * (1.0f - landfraction[x][y]) +
                                weatherphase.landTemperatureAnomalyC * landfraction[x][y] : 0.0f;
                        // Atmospheric latent/cloud adjustments do not directly
                        // change the prescribed land surface or ocean ice skin.
                        phasesurfacetemperature[x][y] = temperature[x][y] + diurnalanomaly;
                        phasetemperature[x][y] = phasesurfacetemperature[x][y] +
                            persistentheating[x][y] -
                            cloudmemory[x][y] *
                                tuning::climate::moistureadvection::maximumCloudRadiativeCoolingC;
                        phasedynamicvertical[x][y] =
                            dynamicvertical[x][y] * polartransportfactors[y] +
                            (tuning::climate::circulation::enableLaggedDiabaticCoupling ? 0.0f : persistentheating[x][y]) *
                                tuning::climate::moistureadvection::heatingVerticalVelocityResponse *
                                polarconvectionfactors[y];
                        phasesaturationcapacity[x][y] =
                            climatephysics::saturationColumnWaterAtPressure(
                                phasetemperature[x][y],
                                surfacepressure[x][y],
                                gravitymultiplier);
                        phaseboundarycapacity[x][y] = phasesaturationcapacity[x][y] *
                            tuning::climate::moistureadvection::boundaryLayerCapacityFraction;
                        surfacewindspeed[x][y] = std::max(
                            tuning::climate::moistureadvection::minimumSurfaceWind,
                            std::sqrt(
                                phasewindu[x][y] * phasewindu[x][y] +
                                phasewindv[x][y] * phasewindv[x][y]));
                    }
                }
            });

            };
            buildphase(startweather, iteration);
            previousphasewindu = phasewindu; previousphasewindv = phasewindv;
            previousphaseupperwindu = phaseupperwindu; previousphaseupperwindv = phaseupperwindv;
            buildphase(weatherstate, iteration + 1);
            auto climatesample = climateweather::makeState(width + 1, height + 1, 2);
            climatesample.elapsedSeconds = (startday * 86400.0) + (iteration + 0.5) * timestepseconds;
            for (int y = 0; y <= height; ++y)
                for (int x = 0; x <= width; ++x)
                {
                    const int cell = y * (width + 1) + x;
                    climatesample.layers[0].eastWindMps[cell] = 0.5f * (previousphasewindu[x][y] + phasewindu[x][y]);
                    climatesample.layers[0].southWindMps[cell] = 0.5f * (previousphasewindv[x][y] + phasewindv[x][y]);
                    climatesample.layers[1].eastWindMps[cell] = 0.5f * (previousphaseupperwindu[x][y] + phaseupperwindu[x][y]);
                    climatesample.layers[1].southWindMps[cell] = 0.5f * (previousphaseupperwindv[x][y] + phaseupperwindv[x][y]);
                }
            seasonalweather[quarter].push_back(std::move(climatesample));
            seasonalweights[quarter].push_back(timestepseconds);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        const bool sea = seacell[x][y] != 0;
                        const float localtemperature = phasetemperature[x][y];
                        const float surfacetemperature = phasesurfacetemperature[x][y];
                        const float retained = std::max(
                            0.0f,
                            boundarymoisture[x][y] + freemoisture[x][y]);
                        // The humidity deficit uses surface saturation and the
                        // actual vapour store, independently of parcel cooling.
                        const float surfacecapacity = climatephysics::saturationColumnWaterAtPressure(
                            surfacetemperature, surfacepressure[x][y], gravitymultiplier);
                        const float relativehumidity = std::clamp(retained / surfacecapacity, 0.0f, 1.0f);
                        const float airtemperature = localtemperature;
                        const float frozencover = sea ? std::clamp(icefraction[x][y], 0.0f, 1.0f)
                            : std::clamp(snowwater[x][y] /
                                tuning::climate::moistureadvection::fullSnowCoverWaterEquivalentMm, 0.0f, 1.0f);
                        const float roughnesslength = sea
                            ? tuning::climate::moistureadvection::oceanRoughnessLengthMetres
                            : tuning::climate::moistureadvection::landRoughnessLengthMetres;
                        const float icecoefficient = climatephysics::neutralDragCoefficient(
                            tuning::climate::moistureadvection::iceRoughnessLengthMetres,
                            tuning::climate::moistureadvection::exchangeReferenceHeightMetres);
                        // Exchange coefficients are area-weighted over open and
                        // frozen surface; trace snow must not switch a whole cell.
                        const float neutralcoefficient =
                            (1.0f - frozencover) * climatephysics::neutralDragCoefficient(
                                roughnesslength,
                                tuning::climate::moistureadvection::exchangeReferenceHeightMetres) +
                            frozencover * icecoefficient;
                        const float referencecoefficient =
                            climatephysics::neutralDragCoefficient(
                                tuning::climate::moistureadvection::oceanRoughnessLengthMetres,
                                tuning::climate::moistureadvection::exchangeReferenceHeightMetres);
                        const float stabilitymultiplier =
                            climatephysics::bulkRichardsonExchangeMultiplier(
                                surfacetemperature,
                                airtemperature,
                                surfacewindspeed[x][y],
                                gravitymultiplier,
                                tuning::climate::moistureadvection::exchangeReferenceHeightMetres,
                                tuning::climate::moistureadvection::minimumStabilityExchangeMultiplier,
                                tuning::climate::moistureadvection::maximumStabilityExchangeMultiplier);
                        const float transfercoefficient =
                            tuning::climate::moistureadvection::surfaceExchangeCoefficient *
                            neutralcoefficient / std::max(0.0001f, referencecoefficient) *
                            stabilitymultiplier;
                        sensibleheat[x][y] = tuning::climate::atmosphere::surfaceAirDensityKgM3 * 1004.0f *
                            transfercoefficient * surfacewindspeed[x][y] * (surfacetemperature - airtemperature);
                        const float potentialevaporation =
                            climatephysics::bulkAerodynamicEvaporationMmAtPressure(
                            surfacetemperature,
                            surfacepressure[x][y],
                            surfacewindspeed[x][y],
                            relativehumidity,
                            timestepseconds,
                            transfercoefficient);
                        float seaevaporation = 0.0f;
                        float landevaporation = 0.0f;

                        if (sea)
                        {
                            seaevaporation = potentialevaporation;

                            seaevaporation *= 1.0f - icefraction[x][y] *
                                (1.0f - tuning::climate::moistureadvection::seaIceFactor);

                            totalseaevaporation[x][y] += seaevaporation;
                        }
                        else
                        {
                            const float melt = climatehydrology::snowMeltAmount(
                                snowwater[x][y],
                                localtemperature,
                                timestepseconds,
                                tuning::climate::moistureadvection::degreeDaySnowMeltMmPerDegreeC);
                            snowwater[x][y] -= melt;
                            const float meltstorageavailable = std::max(
                                0.0f,
                                tuning::climate::moistureadvection::landSoilMoistureCapacity -
                                    soilmoisture[x][y]);
                            const float meltinfiltration = std::min(
                                meltstorageavailable,
                                melt * tuning::climate::moistureadvection::landInfiltrationFraction);
                            soilmoisture[x][y] += meltinfiltration;
                            runoff[x][y] += melt - meltinfiltration;

                            if (snowwater[x][y] > 0.0f)
                            {
                                landevaporation = std::min(
                                    snowwater[x][y],
                                    potentialevaporation *
                                        tuning::climate::moistureadvection::snowSublimationResistance);
                                snowwater[x][y] -= landevaporation;
                            }
                            else if (soilmoisture[x][y] > 0.0f)
                            {
                                const float wateravailability =
                                    climatehydrology::soilMoistureStress(
                                    soilmoisture[x][y],
                                    tuning::climate::moistureadvection::landSoilMoistureCapacity,
                                    tuning::climate::moistureadvection::soilMoistureCriticalFraction,
                                    tuning::climate::moistureadvection::soilMoistureStressExponent);
                                landevaporation = std::min(
                                    soilmoisture[x][y],
                                    potentialevaporation * wateravailability *
                                        tuning::climate::moistureadvection::landSurfaceResistance);
                                soilmoisture[x][y] -= landevaporation;
                            }

                            totallandevaporation[x][y] += landevaporation;
                        }

                        surfaceevaporation[x][y] = seaevaporation + landevaporation;
                        boundarymoisture[x][y] = std::max(0.0f, boundarymoisture[x][y]) +
                            surfaceevaporation[x][y];
                        freemoisture[x][y] = std::max(0.0f, freemoisture[x][y]);
                    }
                }
            });

            advectlayer(
                boundarymoisture,
                phasewindu,
                phasewindv,
                nextboundarymoisture,
                timestepseconds);
            const int uppertimestepcount = iteration %
                tuning::climate::moistureadvection::freeTroposphereTransportCadence + 1;
            const bool updateuppertransport = uppertimestepcount ==
                    tuning::climate::moistureadvection::freeTroposphereTransportCadence ||
                iteration == iterations - 1;

            if (updateuppertransport)
            {
                advectlayer(
                    freemoisture,
                    phaseupperwindu,
                    phaseupperwindv,
                    nextfreemoisture,
                    timestepseconds * static_cast<float>(uppertimestepcount));
                applyeddyexchange(
                    width,
                    height,
                    rowareaweights,
                    nextfreemoisture,
                    zonaleddyexchange,
                    meridionaleddyexchange);
            }
            else
            {
                nextfreemoisture = freemoisture;
            }

            damphighwavenumbers(
                nextboundarymoisture,
                highpasslaplacian,
                highpassscratch,
                rowareaweights,
                polartransportfactors,
                width,
                height,
                tuning::climate::moistureadvection::vapourHighWavenumberDamping,
                tuning::climate::moistureadvection::highWavenumberDampingPasses,
                true);
            damphighwavenumbers(
                nextfreemoisture,
                highpasslaplacian,
                highpassscratch,
                rowareaweights,
                polartransportfactors,
                width,
                height,
                tuning::climate::moistureadvection::vapourHighWavenumberDamping,
                tuning::climate::moistureadvection::highWavenumberDampingPasses,
                true);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {

                        boundaryfluxtendency[x][y] =
                            nextboundarymoisture[x][y] - boundarymoisture[x][y];
                        freefluxtendency[x][y] =
                            nextfreemoisture[x][y] - freemoisture[x][y];
                        totalconvergence[x][y] +=
                            (boundaryfluxtendency[x][y] + freefluxtendency[x][y]) /
                            std::max(0.05f, phasesaturationcapacity[x][y]);
                    }
                }
            });

            smoothconvergencefootprint(
                boundaryfluxtendency,
                convectiveconvergence,
                convergencescratch,
                rowareaweights,
                width,
                height,
                tuning::climate::moistureadvection::convectiveConvergenceMixingFraction,
                tuning::climate::moistureadvection::convectiveConvergenceSmoothingPasses);
            damphighwavenumbers(
                convectiveconvergence,
                highpasslaplacian,
                highpassscratch,
                rowareaweights,
                polartransportfactors,
                width,
                height,
                tuning::climate::moistureadvection::convergenceHighWavenumberDamping,
                tuning::climate::moistureadvection::highWavenumberDampingPasses,
                false);

            parallelforrows(0, height, [&](int startrow, int endrow)
            {
                for (int y = startrow; y <= endrow; y++)
                {
                    for (int x = 0; x <= width; x++)
                    {
                        const bool sea = seacell[x][y] != 0;
                        const float localtemperature = phasetemperature[x][y];
                        const float dynamiccooling = std::clamp(
                            std::max(0.0f, phasedynamicvertical[x][y]) *
                                tuning::climate::moistureadvection::dynamicVerticalCooling,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float dynamicwarming = std::clamp(
                            std::max(0.0f, -phasedynamicvertical[x][y]) *
                                tuning::climate::moistureadvection::dynamicSubsidenceWarming,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float nonorographicparceltemperature =
                            localtemperature - dynamiccooling + dynamicwarming;
                        const float subgridrelief = std::max(
                            0.0f,
                            peakelevation[x][y] - meanelevation[x][y]);
                        const float terraincooling = std::clamp(
                            (uplift[x][y] * polarconvectionfactors[y] +
                                subgridrelief *
                                    tuning::climate::moistureadvection::subgridOrographicExtremeRetention *
                                    polarconvectionfactors[y]) *
                                tuning::climate::moistureadvection::topographicUpliftCoolingCPerMetre,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float terrainwarming = std::clamp(
                            descent[x][y] *
                                tuning::climate::moistureadvection::topographicDescentWarmingCPerMetre,
                            0.0f,
                            tuning::climate::moistureadvection::maximumParcelTemperatureAdjustment);
                        const float terrainadjustedparceltemperature =
                            nonorographicparceltemperature - terraincooling + terrainwarming;
                        const float nonorographiccapacity =
                            climatephysics::saturationColumnWaterAtPressure(
                                nonorographicparceltemperature,
                                surfacepressure[x][y],
                                gravitymultiplier) *
                            tuning::climate::moistureadvection::freeTroposphereCapacityFraction;
                        const float terrainadjustedcapacity =
                            climatephysics::saturationColumnWaterAtPressure(
                                terrainadjustedparceltemperature,
                                surfacepressure[x][y],
                                gravitymultiplier) *
                            tuning::climate::moistureadvection::freeTroposphereCapacityFraction;
                        const float stratiformrelativehumidity = std::clamp(
                            tuning::climate::moistureadvection::stratiformCriticalRelativeHumidity,
                            0.0f,
                            1.0f);
                        const float stratiformnonorographiccapacity =
                            nonorographiccapacity * stratiformrelativehumidity;
                        const float stratiformterrainadjustedcapacity =
                            terrainadjustedcapacity * stratiformrelativehumidity;
                        const float boundaryrelativehumidity = std::clamp(
                            nextboundarymoisture[x][y] /
                                std::max(0.05f, phaseboundarycapacity[x][y]),
                            0.0f,
                            1.0f);
                        const float freerelativehumidity = std::clamp(
                            nextfreemoisture[x][y] /
                                std::max(0.05f, nonorographiccapacity),
                            0.0f,
                            1.0f);
                        const float cloudfraction = std::max(
                            climatehydrology::diagnosticCloudFraction(
                                boundaryrelativehumidity,
                                tuning::climate::moistureadvection::cloudOnsetRelativeHumidity),
                            climatehydrology::diagnosticCloudFraction(
                                freerelativehumidity,
                                tuning::climate::moistureadvection::cloudOnsetRelativeHumidity));
                        const float dayfraction = timestepseconds /
                            tuning::climate::circulation::secondsPerDay;
                        const float backgroundexchange = -std::expm1(
                            -dayfraction /
                            tuning::climate::moistureadvection::backgroundVerticalExchangeTimeDays);
                        const float normalizedvertical = std::clamp(
                            phasedynamicvertical[x][y] /
                                tuning::climate::circulation::maximumVerticalVelocity,
                            -1.0f,
                            1.0f);
                        const float convectivesupplyfraction = std::clamp(
                            std::max(
                                0.0f,
                                convectiveconvergence[x][y] + surfaceevaporation[x][y] +
                                    std::max(0.0f, freefluxtendency[x][y]) *
                                        tuning::climate::moistureadvection::elevatedMoistureAccessionFraction) /
                                std::max(0.05f, phaseboundarycapacity[x][y]),
                            0.0f,
                            1.0f);
                        const float environmentalfreetemperature =
                            temperature[x][y] -
                            tuning::climate::moistureadvection::freeTroposphereEnvironmentalLapseC *
                                polarconvectionfactors[y] +
                            dynamicwarming;
                        const float parcelbuoyancy =
                            localtemperature - environmentalfreetemperature;
                        const float verticalwindshear = std::sqrt(
                            std::pow(phaseupperwindu[x][y] - phasewindu[x][y], 2.0f) +
                            std::pow(phaseupperwindv[x][y] - phasewindv[x][y], 2.0f));
                        const float shallowexchange =
                            climatehydrology::shallowConvectionExchangeFraction(
                                boundaryrelativehumidity,
                                freerelativehumidity,
                                parcelbuoyancy,
                                verticalwindshear,
                                timestepseconds,
                                tuning::climate::moistureadvection::shallowConvectionMixingTimeDays,
                                tuning::climate::moistureadvection::shallowConvectionHumidityOnset,
                                tuning::climate::moistureadvection::shallowConvectionFullHumidity,
                                tuning::climate::moistureadvection::shallowConvectionFullShearMps,
                                tuning::climate::moistureadvection::maximumShallowExchangeFraction);
                        const float dryexchange =
                            climatehydrology::dryConvectionExchangeFraction(
                                parcelbuoyancy,
                                timestepseconds,
                                tuning::climate::moistureadvection::dryConvectionMixingTimeDays,
                                tuning::climate::moistureadvection::dryConvectionActivationBuoyancyC,
                                tuning::climate::moistureadvection::dryConvectionFullStrengthBuoyancyC,
                                tuning::climate::moistureadvection::maximumDryExchangeFraction);
                        const float upwardfraction = std::clamp(
                            backgroundexchange +
                                std::max(0.0f, normalizedvertical) *
                                    tuning::climate::moistureadvection::ascentVerticalExchangeFraction +
                                cloudfraction * convectivesupplyfraction *
                                    tuning::climate::moistureadvection::convectiveVerticalExchangeFraction +
                                shallowexchange + dryexchange,
                            0.0f,
                            tuning::climate::moistureadvection::maximumVerticalExchangeFraction);
                        const float downwardfraction = std::clamp(
                            backgroundexchange +
                                std::max(0.0f, -normalizedvertical) *
                                    tuning::climate::moistureadvection::subsidenceVerticalExchangeFraction +
                                0.5f * shallowexchange + dryexchange,
                            0.0f,
                            tuning::climate::moistureadvection::maximumVerticalExchangeFraction);
                        const climatehydrology::MoistureLayerExchange layerexchange =
                            climatehydrology::exchangeMoistureLayers(
                                nextboundarymoisture[x][y],
                                nextfreemoisture[x][y],
                                upwardfraction,
                                downwardfraction);
                        const float availablemoisture =
                            layerexchange.boundaryLayerMm + layerexchange.freeTroposphereMm;
                        const float condensablewater = std::max(
                            0.0f,
                            layerexchange.freeTroposphereMm -
                                stratiformterrainadjustedcapacity);
                        const climatehydrology::PrecipitationPartition partition =
                            climatehydrology::partitionTwoLayerPrecipitation(
                            layerexchange.boundaryLayerMm,
                            layerexchange.freeTroposphereMm,
                            phaseboundarycapacity[x][y],
                            stratiformnonorographiccapacity,
                            stratiformterrainadjustedcapacity,
                            convectiveconvergence[x][y],
                            surfaceevaporation[x][y],
                            localtemperature,
                            environmentalfreetemperature,
                            updateuppertransport
                                ? timestepseconds * static_cast<float>(uppertimestepcount)
                                : 0.0f,
                            condensationconversiontimeseconds,
                            tuning::climate::moistureadvection::kuoCriticalRelativeHumidity,
                            tuning::climate::moistureadvection::convectiveConversionEfficiency,
                            tuning::climate::moistureadvection::convectiveActivationBuoyancyC,
                            tuning::climate::moistureadvection::convectiveFullStrengthBuoyancyC,
                            tuning::climate::moistureadvection::moistAdjustmentIterations,
                            tuning::climate::moistureadvection::latentHeatingCPerMillimetre,
                            tuning::climate::moistureadvection::saturationCapacityTemperatureSensitivityPerC,
                            freefluxtendency[x][y],
                            tuning::climate::moistureadvection::elevatedMoistureAccessionFraction,
                            tuning::climate::moistureadvection::kuoHumidityExponent);
                        const float condensate = partition.totalMm();
                        const float boundaryaftercondensation = std::max(
                            0.0f,
                            layerexchange.boundaryLayerMm - partition.convectiveMm);
                        const float freeaftercondensation = std::max(
                            0.0f,
                            layerexchange.freeTroposphereMm -
                                partition.stratiformMm - partition.orographicMm);
                        const float falloutrelativehumidity = std::clamp(
                            boundaryaftercondensation /
                                std::max(0.05f, phaseboundarycapacity[x][y]),
                            0.0f,
                            1.0f);
                        const climatehydrology::FallingPrecipitation falling =
                            climatehydrology::processFallingPrecipitation(
                                condensate,
                                localtemperature,
                                falloutrelativehumidity,
                                tuning::climate::moistureadvection::maximumFallingPrecipitationReevaporationFraction,
                                std::max(
                                    0.0f,
                                    phaseboundarycapacity[x][y] - boundaryaftercondensation),
                                tuning::climate::moistureadvection::allSnowTemperatureC,
                                tuning::climate::moistureadvection::allRainTemperatureC);
                        const float precipitation = falling.surfaceTotalMm();
                        const int cell = y * (width + 1) + x;
                        climateatmosphere::ColumnHeatingInput heatinput;
                        heatinput.incomingSolarWm2 = climateenergy::dailyMeanInsolationWm2(climatelatitudes[y], solardeclination, solardistance);
                        const float frozencover = seacell[x][y] ? std::clamp(icefraction[x][y], 0.0f, 1.0f)
                            : std::clamp(snowwater[x][y] /
                                tuning::climate::moistureadvection::fullSnowCoverWaterEquivalentMm, 0.0f, 1.0f);
                        const double openalbedo = seacell[x][y] ? tuning::climate::energybalance::oceanAlbedo
                            : tuning::climate::energybalance::landAlbedo;
                        const double frozenalbedo = seacell[x][y] ? tuning::climate::energybalance::seaIceAlbedo
                            : tuning::climate::energybalance::snowAlbedo;
                        heatinput.surfaceAlbedo = openalbedo + frozencover * (frozenalbedo - openalbedo);
                        heatinput.surfaceTemperatureK = phasesurfacetemperature[x][y] + 273.15;
                        heatinput.airTemperatureK = {localtemperature + 273.15,
                            environmentalfreetemperature + 273.15};
                        heatinput.longwaveOpticalDepth = {tuning::climate::circulation::lowerLongwaveOpticalDepth,
                            tuning::climate::circulation::upperLongwaveOpticalDepth};
                        heatinput.shortwaveOpticalDepth = {tuning::climate::circulation::lowerShortwaveOpticalDepth,
                            tuning::climate::circulation::upperShortwaveOpticalDepth};
                        heatinput.sensibleHeatingWm2 = sensibleheat[x][y];
                        heatinput.condensationMm = {partition.convectiveMm, partition.stratiformMm + partition.orographicMm};
                        heatinput.reevaporationMm = falling.reevaporatedMm;
                        heatinput.surfaceEvaporationMm = surfaceevaporation[x][y];
                        heatinput.accumulationSeconds = timestepseconds;
                        const auto heat = climateatmosphere::diagnoseColumnHeating(heatinput);
                        for (int layer = 0; layer < 2; ++layer)
                        {
                            exchanges.radiativeHeatingWm2[layer][cell] += static_cast<float>(heat.radiativeWm2[layer] * timestepseconds);
                            exchanges.latentHeatingWm2[layer][cell] += static_cast<float>(heat.latentWm2[layer] * timestepseconds);
                        }
                        exchanges.sensibleHeatingWm2[cell] += sensibleheat[x][y] * timestepseconds;
                        exchanges.surfaceNetHeatingWm2[cell] += static_cast<float>(heat.surfaceNetHeatingWm2 * timestepseconds);
                        exchanges.columnEnergyResidualWm2[cell] += static_cast<float>(heat.closureResidualWm2 * timestepseconds);
                        exchanges.columnWaterMm[cell] += (boundaryaftercondensation + falling.reevaporatedMm + freeaftercondensation) * timestepseconds;
                        exchanges.ascentHpaPerDay[cell] += phasedynamicvertical[x][y] * timestepseconds;
                        boundarymoisture[x][y] =
                            boundaryaftercondensation + falling.reevaporatedMm;
                        freemoisture[x][y] = freeaftercondensation;
                        const float heatingdecay = std::exp(
                            -timestepseconds /
                            (tuning::climate::moistureadvection::persistentHeatingDecayTimeDays *
                                tuning::climate::circulation::secondsPerDay));
                        persistentheating[x][y] = std::clamp(
                            persistentheating[x][y] * heatingdecay +
                                precipitation *
                                    tuning::climate::moistureadvection::persistentHeatingRetentionPerMillimetreC,
                            0.0f,
                            tuning::climate::moistureadvection::maximumPersistentHeatingC);
                        const float clouddecay = std::exp(
                            -timestepseconds /
                            (tuning::climate::moistureadvection::cloudMemoryTimeDays *
                                tuning::climate::circulation::secondsPerDay));
                        cloudmemory[x][y] = std::clamp(
                            cloudmemory[x][y] * clouddecay +
                                cloudfraction * (1.0f - clouddecay),
                            0.0f,
                            1.0f);
                        const double areaweight = rowareaweights[y];
                        climatephysics::CondensationActivityDiagnostics& activity =
                            rowcondensation[y];
                        climatephysics::PrecipitationProcessDiagnostics& processes =
                            rowprocesses[y];
                        activity.cellStepAreaWeight += areaweight;
                        activity.atmosphericWaterAreaWeighted +=
                            areaweight * static_cast<double>(availablemoisture);

                        if (precipitation > 0.0f)
                        {
                            activity.activeCellStepAreaWeight += areaweight;
                            activity.activeAtmosphericWaterAreaWeighted +=
                                areaweight * static_cast<double>(availablemoisture);
                            activity.excessWaterAreaWeighted +=
                                areaweight * static_cast<double>(
                                    std::max(condensablewater, precipitation));
                        }

                        activity.precipitationAreaWeighted +=
                            areaweight * static_cast<double>(precipitation);
                        processes.stratiformPrecipitation +=
                            areaweight * static_cast<double>(
                                condensate > 0.0f
                                    ? partition.stratiformMm * precipitation / condensate
                                    : 0.0f);
                        processes.orographicPrecipitation +=
                            areaweight * static_cast<double>(
                                condensate > 0.0f
                                    ? partition.orographicMm * precipitation / condensate
                                    : 0.0f);
                        processes.convectivePrecipitation +=
                            areaweight * static_cast<double>(
                                condensate > 0.0f
                                    ? partition.convectiveMm * precipitation / condensate
                                    : 0.0f);
                        processes.reevaporatedPrecipitation +=
                            areaweight * static_cast<double>(falling.reevaporatedMm);
                        processes.snowfall +=
                            areaweight * static_cast<double>(falling.snowMm);
                        processes.upwardMoistureTransfer +=
                            areaweight * static_cast<double>(layerexchange.upwardTransferMm);
                        processes.downwardMoistureTransfer +=
                            areaweight * static_cast<double>(layerexchange.downwardTransferMm);
                        processes.cloudFractionAreaWeighted += areaweight * cloudfraction;

                        const float combinedfluxtendency =
                            boundaryfluxtendency[x][y] + freefluxtendency[x][y];

                        if (combinedfluxtendency >= 0.0f)
                        {
                            processes.positiveMoistureFluxConvergence +=
                                areaweight * static_cast<double>(combinedfluxtendency);
                        }
                        else
                        {
                            processes.negativeMoistureFluxConvergence +=
                                areaweight * static_cast<double>(combinedfluxtendency);
                        }
                        totalrain[x][y] += precipitation;

                        if (sea == false)
                        {
                            const climatehydrology::SnowAccumulation snowaccumulation =
                                climatehydrology::accumulateSnowfall(
                                    snowwater[x][y],
                                    falling.snowMm,
                                    tuning::climate::moistureadvection::maximumSnowStorageMm);
                            snowwater[x][y] = snowaccumulation.storageMm;
                            const float storageavailable = std::max(
                                0.0f,
                                tuning::climate::moistureadvection::landSoilMoistureCapacity - soilmoisture[x][y]);
                            const float infiltration = std::min(
                                storageavailable,
                                falling.rainMm *
                                    tuning::climate::moistureadvection::landInfiltrationFraction);

                            soilmoisture[x][y] += infiltration;
                            runoff[x][y] += falling.rainMm - infiltration +
                                snowaccumulation.overflowMm;
                        }
                    }
                }
            });
        }

        climatephysics::CondensationActivityDiagnostics condensationactivity;
        climatephysics::PrecipitationProcessDiagnostics processdiagnostics;

        for (const climatephysics::CondensationActivityDiagnostics& row : rowcondensation)
        {
            condensationactivity.cellStepAreaWeight += row.cellStepAreaWeight;
            condensationactivity.activeCellStepAreaWeight += row.activeCellStepAreaWeight;
            condensationactivity.atmosphericWaterAreaWeighted +=
                row.atmosphericWaterAreaWeighted;
            condensationactivity.activeAtmosphericWaterAreaWeighted +=
                row.activeAtmosphericWaterAreaWeighted;
            condensationactivity.excessWaterAreaWeighted += row.excessWaterAreaWeighted;
            condensationactivity.precipitationAreaWeighted +=
                row.precipitationAreaWeighted;
        }

        for (const climatephysics::PrecipitationProcessDiagnostics& row : rowprocesses)
        {
            processdiagnostics.stratiformPrecipitation +=
                row.stratiformPrecipitation;
            processdiagnostics.orographicPrecipitation +=
                row.orographicPrecipitation;
            processdiagnostics.convectivePrecipitation +=
                row.convectivePrecipitation;
            processdiagnostics.reevaporatedPrecipitation +=
                row.reevaporatedPrecipitation;
            processdiagnostics.snowfall += row.snowfall;
            processdiagnostics.upwardMoistureTransfer +=
                row.upwardMoistureTransfer;
            processdiagnostics.downwardMoistureTransfer +=
                row.downwardMoistureTransfer;
            processdiagnostics.cloudFractionAreaWeighted +=
                row.cloudFractionAreaWeighted;
            processdiagnostics.positiveMoistureFluxConvergence +=
                row.positiveMoistureFluxConvergence;
            processdiagnostics.negativeMoistureFluxConvergence +=
                row.negativeMoistureFluxConvergence;
        }

        parallelforrows(0, height, [&](int startrow, int endrow)
        {
            for (int y = startrow; y <= endrow; y++)
            {
                for (int x = 0; x <= width; x++)
                {
                    const float rawrain = totalrain[x][y] *
                        tuning::climate::moistureadvection::rainfallScale;
                    rawannualrain[x][y] += rawrain;
                    quarterrain[x][y] += rawrain;
                    quarterconvergence[x][y] +=
                        totalconvergence[x][y] / static_cast<float>(iterations);

                    if (rawrain > rawmaximummonthlyrain[x][y])
                    {
                        rawmaximummonthlyrain[x][y] = rawrain;
                        rawmaximumrainmonth[x][y] = static_cast<signed char>(month);
                    }

                    if (month % 3 == 2)
                    {
                        const float physicalrain = std::max(
                            0.0f,
                            quarterrain[x][y] / 3.0f);
                        storedseasonalrain[quarter][x][y] = physicalrain;
                        storedseasonalmoisture[quarter][x][y] =
                            seasonalprocess[quarter].columnWaterMm[y * (width + 1) + x] /
                            static_cast<float>(seasonalprocess[quarter].durationSeconds);
                        storedseasonalconvergence[quarter][x][y] =
                            quarterconvergence[x][y] / 3.0f *
                            tuning::climate::moistureadvection::convergenceStorageScale;
                    }
                }
            }
        });

        for (int y = 0; y <= height; y++)
        {
            const double areaweight = rowareaweights[y];

            for (int x = 0; x <= width; x++)
            {
                if (seacell[x][y] != 0)
                {
                    budget.oceanEvaporation += totalseaevaporation[x][y];
                    budget.oceanPrecipitation += totalrain[x][y];
                    areaweightedbudget.oceanEvaporation += areaweight * totalseaevaporation[x][y];
                    areaweightedbudget.oceanPrecipitation += areaweight * totalrain[x][y];
                }
                else
                {
                    budget.landEvaporation += totallandevaporation[x][y];
                    budget.landPrecipitation += totalrain[x][y];
                    budget.runoff += runoff[x][y];
                    budget.soilStorage += soilmoisture[x][y];
                    budget.snowStorage += snowwater[x][y];
                    areaweightedbudget.landEvaporation += areaweight * totallandevaporation[x][y];
                    areaweightedbudget.landPrecipitation += areaweight * totalrain[x][y];
                    areaweightedbudget.runoff += areaweight * runoff[x][y];
                    areaweightedbudget.soilStorage += areaweight * soilmoisture[x][y];
                    areaweightedbudget.snowStorage += areaweight * snowwater[x][y];
                }

                const float atmosphericwater =
                    boundarymoisture[x][y] + freemoisture[x][y];
                budget.atmosphericStorage += atmosphericwater;
                areaweightedbudget.atmosphericStorage += areaweight * atmosphericwater;
            }
        }

        continuedweather = weatherstate;
        return monthlysolverresult{
            budget,
            areaweightedbudget,
            condensationactivity,
            processdiagnostics,
            climateweather::serializeState(weatherstate)
        };
    };

    climatephysics::HydrologySpinupDiagnostics diagnostics;

    for (int cycle = 0; cycle < tuning::climate::moistureadvection::maximumSpinupCycles; cycle++)
    {
        // Replay the same evolving year while hydrological storage spins up.
        // Months within a year share one dynamical state, not independent noise.
        continuedweather = {};
        for (auto& samples : seasonalweather) samples.clear();
        for (auto& weights : seasonalweights) weights.clear();
        for (auto& fields : seasonalprocess)
        {
            fields = {};
            fields.columns = width + 1; fields.rows = height + 1;
            const int cells = fields.columns * fields.rows;
            for (int layer = 0; layer < 2; ++layer)
            {
                fields.eastIntegratedFlux[layer].assign(cells, 0.0);
                fields.southIntegratedFlux[layer].assign(cells, 0.0);
                fields.radiativeHeatingWm2[layer].assign(cells, 0.0f);
                fields.latentHeatingWm2[layer].assign(cells, 0.0f);
            }
            fields.sensibleHeatingWm2.assign(cells, 0.0f);
            fields.surfaceNetHeatingWm2.assign(cells, 0.0f);
            fields.columnWaterMm.assign(cells, 0.0f);
            fields.ascentHpaPerDay.assign(cells, 0.0f);
            fields.columnEnergyResidualWm2.assign(cells, 0.0f);
        }
        const auto cyclestart = std::chrono::steady_clock::now();
        const hydrologyfloatgrid initialboundarymoisture = boundarymoisture;
        const hydrologyfloatgrid initialfreemoisture = freemoisture;
        const hydrologyfloatgrid initialsoilmoisture = soilmoisture;
        const hydrologyfloatgrid initialsnowwater = snowwater;

        rawannualrain.fill(0.0f);
        rawmaximummonthlyrain.fill(0.0f);
        rawmaximumrainmonth.fill(static_cast<signed char>(-1));

        finalbudgets = {};
        finalareaweightedbudgets = {};
        finalcondensationactivity = {};
        finalprocessdiagnostics = {};

        auto appendbudget = [](
            climatephysics::WaterBudget& quarterbudget,
            const climatephysics::WaterBudget& monthbudget,
            bool firstmonth)
        {
            if (firstmonth)
            {
                quarterbudget.initialAtmosphericStorage =
                    monthbudget.initialAtmosphericStorage;
                quarterbudget.initialSoilStorage = monthbudget.initialSoilStorage;
                quarterbudget.initialSnowStorage = monthbudget.initialSnowStorage;
            }

            quarterbudget.oceanEvaporation += monthbudget.oceanEvaporation;
            quarterbudget.landEvaporation += monthbudget.landEvaporation;
            quarterbudget.oceanPrecipitation += monthbudget.oceanPrecipitation;
            quarterbudget.landPrecipitation += monthbudget.landPrecipitation;
            quarterbudget.runoff += monthbudget.runoff;
            quarterbudget.atmosphericStorage = monthbudget.atmosphericStorage;
            quarterbudget.soilStorage = monthbudget.soilStorage;
            quarterbudget.snowStorage = monthbudget.snowStorage;
        };

        for (int month = 0; month < climatehydrology::monthCount; month++)
        {
            const monthlysolverresult result = runsolver(month);
            if (month == climatehydrology::monthCount - 1)
                finalweatherstate = result.weatherstate;
            const int quarter = month / 3;
            appendbudget(finalbudgets[quarter], result.budget, month % 3 == 0);
            appendbudget(
                finalareaweightedbudgets[quarter],
                result.areaweightedbudget,
                month % 3 == 0);

            auto& activity = finalcondensationactivity[quarter];
            activity.cellStepAreaWeight += result.condensation.cellStepAreaWeight;
            activity.activeCellStepAreaWeight +=
                result.condensation.activeCellStepAreaWeight;
            activity.atmosphericWaterAreaWeighted +=
                result.condensation.atmosphericWaterAreaWeighted;
            activity.activeAtmosphericWaterAreaWeighted +=
                result.condensation.activeAtmosphericWaterAreaWeighted;
            activity.excessWaterAreaWeighted +=
                result.condensation.excessWaterAreaWeighted;
            activity.precipitationAreaWeighted +=
                result.condensation.precipitationAreaWeighted;

            auto& processes = finalprocessdiagnostics[quarter];
            processes.stratiformPrecipitation +=
                result.processes.stratiformPrecipitation;
            processes.orographicPrecipitation +=
                result.processes.orographicPrecipitation;
            processes.convectivePrecipitation +=
                result.processes.convectivePrecipitation;
            processes.reevaporatedPrecipitation +=
                result.processes.reevaporatedPrecipitation;
            processes.snowfall += result.processes.snowfall;
            processes.upwardMoistureTransfer +=
                result.processes.upwardMoistureTransfer;
            processes.downwardMoistureTransfer +=
                result.processes.downwardMoistureTransfer;
            processes.cloudFractionAreaWeighted +=
                result.processes.cloudFractionAreaWeighted;
            processes.positiveMoistureFluxConvergence +=
                result.processes.positiveMoistureFluxConvergence;
            processes.negativeMoistureFluxConvergence +=
                result.processes.negativeMoistureFluxConvergence;
        }

        double atmosphericabsolutechange = 0.0;
        double atmosphericreferencestorage = 0.0;
        double soilabsolutechange = 0.0;
        double soilreferencestorage = 0.0;
        double snowabsolutechange = 0.0;
        double snowreferencestorage = 0.0;
        double snowcoverchangedarea = 0.0;
        double landarea = 0.0;
        double atmosphericstorage = 0.0;
        double soilstorage = 0.0;
        double snowstorage = 0.0;

        for (int y = 0; y <= height; y++)
        {
            const double areaweight = rowareaweights[y];

            for (int x = 0; x <= width; x++)
            {
                const double atmosphericwater =
                    boundarymoisture[x][y] + freemoisture[x][y];
                const double initialatmosphericwater =
                    initialboundarymoisture[x][y] + initialfreemoisture[x][y];
                atmosphericabsolutechange += areaweight * std::abs(
                    atmosphericwater - initialatmosphericwater);
                atmosphericreferencestorage += areaweight * 0.5 *
                    (atmosphericwater + initialatmosphericwater);
                atmosphericstorage += areaweight * atmosphericwater;

                if (seacell[x][y] == 0)
                {
                    landarea += areaweight;
                    soilabsolutechange += areaweight * std::abs(static_cast<double>(
                        soilmoisture[x][y] - initialsoilmoisture[x][y]));
                    soilreferencestorage += areaweight * 0.5 * static_cast<double>(
                        soilmoisture[x][y] + initialsoilmoisture[x][y]);
                    soilstorage += areaweight * soilmoisture[x][y];
                    snowabsolutechange += areaweight * std::abs(static_cast<double>(
                        snowwater[x][y] - initialsnowwater[x][y]));
                    snowreferencestorage += areaweight * 0.5 * static_cast<double>(
                        snowwater[x][y] + initialsnowwater[x][y]);
                    snowstorage += areaweight * snowwater[x][y];
                    const bool hassnow = snowwater[x][y] >=
                        tuning::climate::moistureadvection::snowCoverConvergenceThresholdMm;
                    const bool initiallyhadsnow = initialsnowwater[x][y] >=
                        tuning::climate::moistureadvection::snowCoverConvergenceThresholdMm;

                    if (hassnow != initiallyhadsnow)
                        snowcoverchangedarea += areaweight;
                }
            }
        }

        diagnostics.cyclesCompleted = cycle + 1;
        diagnostics.relativeAtmosphericStorageChange = atmosphericabsolutechange /
            std::max(1.0, atmosphericreferencestorage);
        diagnostics.relativeSoilStorageChange = soilabsolutechange /
            std::max(1.0, soilreferencestorage);
        diagnostics.relativeSnowStorageChange = snowabsolutechange /
            std::max(1.0, snowreferencestorage);
        diagnostics.relativeSnowCoverChange = snowcoverchangedarea /
            std::max(1.0, landarea);
        diagnostics.relativeStorageChange = std::max(
            diagnostics.relativeAtmosphericStorageChange,
            std::max(
                diagnostics.relativeSoilStorageChange,
                std::max(diagnostics.relativeSnowCoverChange, diagnostics.relativeSnowStorageChange)));
        diagnostics.atmosphericStorage = atmosphericstorage;
        diagnostics.soilStorage = soilstorage;
        diagnostics.snowStorage = snowstorage;
        diagnostics.converged = diagnostics.cyclesCompleted >=
                tuning::climate::moistureadvection::minimumSpinupCycles &&
            diagnostics.relativeStorageChange <=
                tuning::climate::moistureadvection::spinupRelativeStorageTolerance;
        const double cycleelapsedseconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - cyclestart).count();

        std::cout
            << "Climate hydrology spinup cycle=" << diagnostics.cyclesCompleted
            << " relative_storage_change=" << diagnostics.relativeStorageChange
            << " relative_atmospheric_change="
            << diagnostics.relativeAtmosphericStorageChange
            << " relative_soil_change=" << diagnostics.relativeSoilStorageChange
            << " relative_snow_change=" << diagnostics.relativeSnowStorageChange
            << " relative_snow_cover_change=" << diagnostics.relativeSnowCoverChange
            << " atmospheric_storage=" << diagnostics.atmosphericStorage
            << " soil_storage=" << diagnostics.soilStorage
            << " snow_storage=" << diagnostics.snowStorage
            << " elapsed_seconds=" << cycleelapsedseconds
            << " converged=" << (diagnostics.converged ? 1 : 0)
            << '\n';

        if (diagnostics.converged)
            break;
    }

    hydrologyfloatgrid outputrawannualrain(
        outputwidth + 1,
        outputheight + 1,
        0.0f);
    hydrologyfloatgrid outputrawmaximummonthlyrain(
        outputwidth + 1,
        outputheight + 1,
        0.0f);
    hydrologybytegrid outputrawmaximumrainmonth(
        outputwidth + 1,
        outputheight + 1,
        static_cast<signed char>(-1));
    vector<float> outputrowareaweights(outputheight + 1, 0.0f);

    for (int y = 0; y <= outputheight; y++)
        outputrowareaweights[y] = static_cast<float>(rowareaweight(y, outputheight));

    parallelforrows(0, outputheight, [&](int startrow, int endrow)
    {
        for (int y = startrow; y <= endrow; y++)
        {
            const int nearestsourcey = outputheight > 0
                ? std::clamp(
                    static_cast<int>(std::round(
                        static_cast<float>(y * height) /
                        static_cast<float>(outputheight))),
                    0,
                    height)
                : 0;

            for (int x = 0; x <= outputwidth; x++)
            {
                outputrawannualrain[x][y] = samplehydrologygrid(
                    rawannualrain,
                    width,
                    height,
                    x,
                    y,
                    outputwidth,
                    outputheight);
                outputrawmaximummonthlyrain[x][y] = samplehydrologygrid(
                    rawmaximummonthlyrain,
                    width,
                    height,
                    x,
                    y,
                    outputwidth,
                    outputheight);
                const float sourcex = (static_cast<float>(x) + 0.5f) *
                    static_cast<float>(width + 1) /
                    static_cast<float>(outputwidth + 1) - 0.5f;
                const int nearestsourcex = wrapx(
                    static_cast<int>(std::round(sourcex)),
                    width);
                outputrawmaximumrainmonth[x][y] =
                    rawmaximumrainmonth[nearestsourcex][nearestsourcey];

                for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                {
                    const float seasonalrain = std::max(
                        0.0f,
                        samplehydrologygrid(
                            storedseasonalrain[season],
                            width,
                            height,
                            x,
                            y,
                            outputwidth,
                            outputheight));
                    const float seasonalmoisture = samplehydrologygrid(
                        storedseasonalmoisture[season],
                        width,
                        height,
                        x,
                        y,
                        outputwidth,
                        outputheight);
                    const float seasonalconvergence = samplehydrologygrid(
                        storedseasonalconvergence[season],
                        width,
                        height,
                        x,
                        y,
                        outputwidth,
                        outputheight);
                    world.setseasonalrainfloat(season, x, y, seasonalrain);
                    world.setseasonalmoisture(
                        season,
                        x,
                        y,
                        std::clamp(
                            static_cast<int>(std::round(seasonalmoisture)),
                            static_cast<int>(std::numeric_limits<short>::min()),
                            static_cast<int>(std::numeric_limits<short>::max())));
                    world.setseasonalconvergence(
                        season,
                        x,
                        y,
                        std::clamp(
                            static_cast<int>(std::round(seasonalconvergence)),
                            static_cast<int>(std::numeric_limits<short>::min()),
                            static_cast<int>(std::numeric_limits<short>::max())));

                    const int roundedrain = std::clamp(
                        static_cast<int>(std::round(seasonalrain)),
                        0,
                        static_cast<int>(std::numeric_limits<short>::max()));
                    if (season == seasonjanuary)
                        world.setjanrain(x, y, roundedrain);
                    else if (season == seasonjuly)
                        world.setjulrain(x, y, roundedrain);
                }
            }
        }
    });

    struct distributionaccumulator
    {
        climatephysics::PrecipitationDistributionScope metrics;
        long long rawzeros = 0;
        long long storedzeros = 0;
        long long belowone = 0;
        double areaweightedrawzeros = 0.0;
        double areaweightedstoredzeros = 0.0;
        double areaweightedbelowone = 0.0;
        double precipitationtotal = 0.0;
        double areaweightedprecipitationtotal = 0.0;
        vector<float> precipitation;
    };

    distributionaccumulator landdistribution;
    distributionaccumulator oceandistribution;

    for (int y = 0; y <= outputheight; y++)
    {
        const double areaweight = outputrowareaweights[y];

        for (int x = 0; x <= outputwidth; x++)
        {
            distributionaccumulator& distribution = world.sea(x, y) == 1
                ? oceandistribution
                : landdistribution;
            const float meanmonthlyrain = outputrawannualrain[x][y] /
                static_cast<float>(climatehydrology::monthCount);
            int storedseasonalrain = 0;

            for (int season = 0; season < CLIMATESEASONCOUNT; season++)
                storedseasonalrain += world.seasonalrain(season, x, y);

            distribution.metrics.cells++;
            distribution.metrics.areaWeight += areaweight;
            distribution.precipitationtotal += meanmonthlyrain;
            distribution.areaweightedprecipitationtotal += areaweight * meanmonthlyrain;
            if (outputrawmaximummonthlyrain[x][y] >
                distribution.metrics.maximumMonthlyPrecipitationMm)
            {
                const int maximummonth = static_cast<int>(outputrawmaximumrainmonth[x][y]);
                distribution.metrics.maximumMonthlyPrecipitationMm =
                    static_cast<double>(outputrawmaximummonthlyrain[x][y]);
                distribution.metrics.maximumMonthlyPrecipitationX = x;
                distribution.metrics.maximumMonthlyPrecipitationY = y;
                distribution.metrics.maximumMonthlyPrecipitationMonth = maximummonth;
                distribution.metrics.roundedPrecipitationAtMaximumMm =
                    maximummonth >= 0
                        ? static_cast<int>(std::round(outputrawmaximummonthlyrain[x][y]))
                        : 0;
            }
            distribution.precipitation.push_back(meanmonthlyrain);

            if (outputrawannualrain[x][y] <= 1.0e-6f)
            {
                distribution.rawzeros++;
                distribution.areaweightedrawzeros += areaweight;
            }

            if (storedseasonalrain == 0)
            {
                distribution.storedzeros++;
                distribution.areaweightedstoredzeros += areaweight;
            }

            if (meanmonthlyrain < 1.0f)
            {
                distribution.belowone++;
                distribution.areaweightedbelowone += areaweight;
            }
        }
    }

    auto finalisedistribution = [](distributionaccumulator& distribution)
    {
        auto& metrics = distribution.metrics;
        const double inversecellcount = metrics.cells > 0
            ? 1.0 / static_cast<double>(metrics.cells)
            : 0.0;
        const double inverseareaweight = metrics.areaWeight > 0.0
            ? 1.0 / metrics.areaWeight
            : 0.0;

        metrics.rawZeroFraction = distribution.rawzeros * inversecellcount;
        metrics.storedZeroFraction = distribution.storedzeros * inversecellcount;
        metrics.belowOneMillimetreFraction = distribution.belowone * inversecellcount;
        metrics.areaWeightedRawZeroFraction = distribution.areaweightedrawzeros * inverseareaweight;
        metrics.areaWeightedStoredZeroFraction = distribution.areaweightedstoredzeros * inverseareaweight;
        metrics.areaWeightedBelowOneMillimetreFraction = distribution.areaweightedbelowone * inverseareaweight;
        metrics.meanMonthlyPrecipitationMm = distribution.precipitationtotal * inversecellcount;
        metrics.areaWeightedMeanMonthlyPrecipitationMm =
            distribution.areaweightedprecipitationtotal * inverseareaweight;

        std::sort(distribution.precipitation.begin(), distribution.precipitation.end());
        const size_t wettestcount = static_cast<size_t>(std::ceil(
            static_cast<double>(distribution.precipitation.size()) * 0.10));
        double wettesttotal = 0.0;

        for (size_t index = distribution.precipitation.size() - wettestcount;
            index < distribution.precipitation.size(); index++)
        {
            wettesttotal += distribution.precipitation[index];
        }

        if (distribution.precipitationtotal > 0.0)
            metrics.wettestTenPercentShare = wettesttotal / distribution.precipitationtotal;
    };

    finalisedistribution(landdistribution);
    finalisedistribution(oceandistribution);

    climatephysics::PrecipitationDistributionDiagnostics precipitationdiagnostics;
    precipitationdiagnostics.land = landdistribution.metrics;
    precipitationdiagnostics.ocean = oceandistribution.metrics;
    climatephysics::setLastPrecipitationDistributionDiagnostics(precipitationdiagnostics);
    climatephysics::setLastHydrologySpinupDiagnostics(diagnostics);

    std::cout
        << "Climate precipitation distribution land_raw_zero="
        << precipitationdiagnostics.land.areaWeightedRawZeroFraction
        << " land_stored_zero="
        << precipitationdiagnostics.land.areaWeightedStoredZeroFraction
        << " land_below_1mm="
        << precipitationdiagnostics.land.areaWeightedBelowOneMillimetreFraction
        << " land_wettest_10_percent_share="
        << precipitationdiagnostics.land.wettestTenPercentShare
        << '\n';

    for (int season = 0; season < CLIMATESEASONCOUNT; season++)
    {
        const climatephysics::WaterBudget& budget = finalbudgets[season];
        const climatephysics::WaterBudget& areaweightedbudget = finalareaweightedbudgets[season];
        climatephysics::setLastWaterBudget(season, budget);
        climatephysics::setLastAreaWeightedWaterBudget(season, areaweightedbudget);
        climatephysics::setLastCondensationActivityDiagnostics(
            season,
            finalcondensationactivity[season]);
        climatephysics::setLastPrecipitationProcessDiagnostics(
            season,
            finalprocessdiagnostics[season]);
        std::cout
            << "Climate area-weighted water budget season=" << season
            << " initial_atmospheric_storage=" << areaweightedbudget.initialAtmosphericStorage
            << " initial_soil_storage=" << areaweightedbudget.initialSoilStorage
            << " initial_snow_storage=" << areaweightedbudget.initialSnowStorage
            << " ocean_evaporation=" << areaweightedbudget.oceanEvaporation
            << " land_evaporation=" << areaweightedbudget.landEvaporation
            << " ocean_precipitation=" << areaweightedbudget.oceanPrecipitation
            << " land_precipitation=" << areaweightedbudget.landPrecipitation
            << " runoff=" << areaweightedbudget.runoff
            << " atmospheric_storage=" << areaweightedbudget.atmosphericStorage
            << " soil_storage=" << areaweightedbudget.soilStorage
            << " snow_storage=" << areaweightedbudget.snowStorage
            << " residual=" << areaweightedbudget.residual()
            << " relative_residual=" << areaweightedbudget.relativeResidual()
            << '\n';
    }

    // Climate sampling does not initialize/advance the deferred gameplay state.
    climatestorage = { &world, width + 1, height + 1, boundarymoisture, freemoisture, soilmoisture, snowwater };
    processfieldworld = &world;
    for (int season = 0; season < CLIMATESEASONCOUNT; ++season)
    {
        auto& fields = seasonalprocess[season];
        const auto normalize = [&](vector<float>& values) {
            for (float& value : values) value /= static_cast<float>(std::max(1.0, fields.durationSeconds)); };
        for (int layer = 0; layer < 2; ++layer)
        {
            normalize(fields.radiativeHeatingWm2[layer]); normalize(fields.latentHeatingWm2[layer]);
        }
        normalize(fields.sensibleHeatingWm2); normalize(fields.surfaceNetHeatingWm2);
        normalize(fields.columnWaterMm); normalize(fields.ascentHpaPerDay); normalize(fields.columnEnergyResidualWm2);
        climatevalidation::captureclimateprocessfields(world, season, fields);
        processfields[season] = std::move(fields);
        for (int layer = 0; layer < 2; ++layer)
            climatevalidation::captureweatherstatistics(world, season, width + 1, height + 1,
                climateweather::calculateStatistics(seasonalweather[season], layer, seasonalweights[season]), layer);
    }
}
