#include "simulations/geology/mineral_systems/gold_formation.hpp"

#include "simulations/geology/mineral_systems/mineral_random.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace mineralsystems
{
namespace
{
constexpr double kPi = 3.14159265358979323846;
constexpr double kKernelFloor = 1.0e-3; // Gaussian kernels are truncated below this.
constexpr double kLognormalTailSigmas = 4.0;
constexpr int kInheritedPlacementAttempts = 4;

double smoothstep(double edge0, double edge1, double value)
{
    if (edge1 <= edge0)
        return value >= edge1 ? 1.0 : 0.0;
    const double t = std::clamp((value - edge0) / (edge1 - edge0), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

bool hostcell(const TectonicFrame& frame, std::size_t cell, double minimumElevationM)
{
    return frame.elevationM != nullptr && frame.elevationM[cell] >= minimumElevationM;
}

struct PlateEdge
{
    std::uint32_t cells = 0;
    double continentalCrustM = 0.0;
    double crustAgeMyr = 0.0;

    bool buoyant(double thresholdM) const { return cells > 0 && continentalCrustM >= thresholdM * cells; }
};

// Crust within one cell of each convergent segment, on its left and right plates.
std::vector<std::array<PlateEdge, 2>> convergentedges(const TectonicFrame& frame)
{
    std::vector<std::array<PlateEdge, 2>> edges(frame.boundaries.size());
    if (frame.continentalCrustM == nullptr)
        return edges;
    for (std::size_t cell = 0; cell < frame.cellCount(); ++cell)
    {
        const std::uint32_t id = frame.nearestBoundaryId[cell];
        if (frame.boundaryDistance[cell] > 1 || id >= frame.boundaries.size())
            continue;
        const BoundarySummary& boundary = frame.boundaries[id];
        if (boundary.boundaryType != boundarycode::convergent)
            continue;
        const std::uint32_t plate = frame.plateIndex[cell];
        const int side = plate == boundary.leftPlate ? 0 : plate == boundary.rightPlate ? 1 : -1;
        if (side < 0)
            continue;
        PlateEdge& edge = edges[id][side];
        ++edge.cells;
        edge.continentalCrustM += frame.continentalCrustM[cell];
        edge.crustAgeMyr += frame.crustAgeMyr != nullptr ? frame.crustAgeMyr[cell] : 0.0;
    }
    return edges;
}

// The plate that overrides a subduction margin, or kNoPlate when the margin
// is a collision or its polarity is unknown.
std::uint32_t overridingplate(const BoundarySummary& boundary, const PlateEdge& left, const PlateEdge& right,
    double buoyantCrustM)
{
    const bool leftbuoyant = left.buoyant(buoyantCrustM);
    const bool rightbuoyant = right.buoyant(buoyantCrustM);
    if (leftbuoyant != rightbuoyant)
        return leftbuoyant ? boundary.leftPlate : boundary.rightPlate;
    if (leftbuoyant || left.cells == 0 || right.cells == 0)
        return kNoPlate;
    // Older, colder ocean floor subducts.
    return left.crustAgeMyr / left.cells <= right.crustAgeMyr / right.cells ? boundary.leftPlate : boundary.rightPlate;
}

// Boundary cells border the plate edge, so a cell at distance d has its centre
// (d + 0.5) cells from the edge; the midpoint keeps kernel integrals independent
// of cell size.
std::vector<float> gaussiankernel(double cellKm, double centreKm, double widthKm)
{
    std::vector<float> kernel;
    if (cellKm <= 0.0 || widthKm <= 0.0)
        return kernel;
    const double cutoffKm = centreKm + widthKm * std::sqrt(-std::log(kKernelFloor));
    kernel.resize(static_cast<std::size_t>(std::ceil(cutoffKm / cellKm)) + 1);
    for (std::size_t distance = 0; distance < kernel.size(); ++distance)
    {
        const double z = ((static_cast<double>(distance) + 0.5) * cellKm - centreKm) / widthKm;
        const double value = std::exp(-z * z);
        kernel[distance] = value < kKernelFloor ? 0.0f : static_cast<float>(value);
    }
    return kernel;
}

float kernelat(const std::vector<float>& kernel, std::uint16_t distance)
{
    return distance < kernel.size() ? kernel[distance] : 0.0f;
}

double lognormal(detail::Random& random, const LognormalRange& range)
{
    if (range.median <= 0.0)
        return 0.0;
    const double z = std::clamp(random.normal(), -kLognormalTailSigmas, kLognormalTailSigmas);
    return range.median * std::exp(range.sigmaLn * z);
}

double wrapcoordinate(double value, std::uint32_t period)
{
    const double size = static_cast<double>(period);
    double wrapped = std::fmod(value, size);
    if (wrapped < 0.0)
        wrapped += size;
    return wrapped >= size ? 0.0 : wrapped;
}

// Cells weighted for area-proportional sampling.
class WeightedCells
{
public:
    void add(std::size_t cell, double weight)
    {
        cells_.push_back(cell);
        cumulative_.push_back(total() + weight);
    }
    double total() const { return cumulative_.empty() ? 0.0 : cumulative_.back(); }
    std::size_t pick(double unit) const
    {
        const auto found = std::upper_bound(cumulative_.begin(), cumulative_.end(), unit * total());
        const auto position = static_cast<std::size_t>(found - cumulative_.begin());
        return cells_[std::min(position, cells_.size() - 1)];
    }

private:
    std::vector<std::size_t> cells_;
    std::vector<double> cumulative_;
};
}

SettingField::SettingField(const TectonicFrame& frame, const FormationParameters& parameters)
    : frame_(frame), minimumHostElevationM_(parameters.minimumHostElevationM)
{
    const double cellKm = frame.cellKm();
    orogenicKernel_ = gaussiankernel(cellKm, 0.0, parameters.orogenicWidthKm);
    arcKernel_ = gaussiankernel(cellKm, parameters.arcPeakKm, parameters.arcWidthKm);
    boundaries_.assign(frame.boundaries.size(), BoundarySetting{});
    const auto edges = convergentedges(frame);

    const double full = std::max(1.0, parameters.fullActivityScore);
    for (std::size_t id = 0; id < frame.boundaries.size(); ++id)
    {
        const BoundarySummary& boundary = frame.boundaries[id];
        const double convergence = std::clamp(boundary.convergenceScore / full, 0.0, 1.0);
        const double shear = std::clamp(boundary.shearScore / full, 0.0, 1.0);
        const double maturity = smoothstep(parameters.orogenicOnsetMyr, parameters.orogenicMatureMyr, boundary.ageMyr);
        BoundarySetting& setting = boundaries_[id];
        if (boundary.boundaryType == boundarycode::convergent)
        {
            const PlateEdge& left = edges[id][0];
            const PlateEdge& right = edges[id][1];
            const std::uint32_t overriding = overridingplate(boundary, left, right, parameters.buoyantCrustM);
            if (left.buoyant(parameters.buoyantCrustM) && right.buoyant(parameters.buoyantCrustM))
            {
                setting.regime = regimecode::continent_collision;
                setting.orogenic = static_cast<float>(parameters.collisionWeight * convergence * maturity);
            }
            else if (overriding != kNoPlate)
            {
                setting.regime = regimecode::convergent_arc;
                setting.plate = overriding;
                setting.orogenic = static_cast<float>(parameters.accretionaryWeight * convergence * maturity);
                setting.arc = static_cast<float>(convergence);
            }
        }
        else if (boundary.boundaryType == boundarycode::transform)
        {
            setting.regime = regimecode::transform;
            setting.orogenic = static_cast<float>(parameters.transformWeight * shear * maturity);
        }
    }
}

SettingIntensity SettingField::at(std::size_t cell) const
{
    if (!hostcell(frame_, cell, minimumHostElevationM_))
        return {};
    const std::uint32_t id = frame_.nearestBoundaryId[cell];
    if (id >= boundaries_.size())
        return {};
    const BoundarySetting& setting = boundaries_[id];
    if (setting.plate != kNoPlate && frame_.plateIndex[cell] != setting.plate)
        return {};
    const std::uint16_t distance = frame_.boundaryDistance[cell];
    return { setting.orogenic * kernelat(orogenicKernel_, distance), setting.arc * kernelat(arcKernel_, distance) };
}

std::uint8_t SettingField::regime(std::uint32_t boundaryId) const
{
    return boundaryId < boundaries_.size() ? boundaries_[boundaryId].regime : regimecode::stable;
}

void sampledepositbody(Deposit& deposit, const SystemParameters& parameters, std::uint64_t seed)
{
    detail::Random random(seed);
    deposit.oreMt = lognormal(random, parameters.oreMt);
    deposit.goldGradeGt = lognormal(random, parameters.goldGradeGt);
    deposit.copperGradePct = lognormal(random, parameters.copperGradePct);
    deposit.formationDepthM = random.uniform(parameters.formationDepthM.minimum, parameters.formationDepthM.maximum);
    deposit.verticalExtentM = random.uniform(parameters.verticalExtentM.minimum, parameters.verticalExtentM.maximum);
}

FormationSample sampleformation(const TectonicFrame& frame, const FormationParameters& parameters, double dtMyr,
    std::uint64_t seed)
{
    FormationSample sample;
    if (dtMyr <= 0.0 || frame.cellCount() == 0)
        return sample;

    const SettingField field(frame, parameters);
    WeightedCells orogenic;
    WeightedCells arc;
    for (std::uint32_t y = 0; y < frame.height; ++y)
    {
        const double area = frame.cellAreaKm2();
        for (std::uint32_t x = 0; x < frame.width; ++x)
        {
            const std::size_t cell = frame.cellIndex(x, y);
            const SettingIntensity setting = field.at(cell);
            if (setting.orogenic > 0.0f)
                orogenic.add(cell, area * setting.orogenic);
            if (setting.arc > 0.0f)
                arc.add(cell, area * setting.arc);
        }
    }
    sample.orogenicAreaKm2 = orogenic.total();
    sample.arcAreaKm2 = arc.total();

    detail::Random random(seed);
    const double cellKm = frame.cellKm();
    const auto form = [&](DepositSystem system, const SystemParameters& systemparameters, const WeightedCells& cells)
    {
        const double expected = systemparameters.ratePerMkm2Myr * 1.0e-6 * cells.total() * dtMyr;
        const std::uint64_t count = random.poisson(expected);
        for (std::uint64_t k = 0; k < count; ++k)
        {
            const std::size_t cell = cells.pick(random.uniform());
            NewDeposit item;
            item.x = static_cast<double>(cell % frame.width) + random.uniform();
            item.y = static_cast<double>(cell / frame.width) + random.uniform();

            Deposit& deposit = item.deposit;
            deposit.system = system;
            deposit.origin = DepositOrigin::simulated;
            sampledepositbody(deposit, systemparameters, random.next());

            const std::uint32_t boundaryid = frame.nearestBoundaryId[cell];
            const BoundarySummary boundary = boundaryid < frame.boundaries.size() ? frame.boundaries[boundaryid] : BoundarySummary{};
            const SettingIntensity setting = field.at(cell);
            FormationRecord& formation = deposit.formation;
            formation.timeMyr = frame.timeMyr - dtMyr * random.uniform();
            formation.geologicRegime = field.regime(boundaryid);
            formation.boundaryType = boundary.boundaryType;
            formation.hostCrustClass = frame.crustClass != nullptr ? frame.crustClass[cell] : crustcode::none;
            formation.boundarySegmentId = boundaryid;
            formation.boundaryDistanceKm = static_cast<float>((frame.boundaryDistance[cell] + 0.5) * cellKm);
            formation.boundaryAgeMyr = boundary.ageMyr;
            formation.settingIntensity = system == DepositSystem::orogenic_gold ? setting.orogenic : setting.arc;
            sample.deposits.push_back(item);
        }
    };
    form(DepositSystem::orogenic_gold, parameters.orogenic, orogenic);
    form(DepositSystem::epithermal_gold, parameters.epithermal, arc);
    form(DepositSystem::porphyry_copper_gold, parameters.porphyry, arc);
    return sample;
}

std::vector<NewDeposit> sampleinheritedprovinces(const TectonicFrame& frame, const FormationParameters& parameters,
    std::uint64_t seed, std::uint32_t& provinceCount)
{
    std::vector<NewDeposit> deposits;
    provinceCount = 0;
    WeightedCells hosts;
    for (std::uint32_t y = 0; y < frame.height; ++y)
    {
        const double area = frame.cellAreaKm2();
        for (std::uint32_t x = 0; x < frame.width; ++x)
        {
            const std::size_t cell = frame.cellIndex(x, y);
            if (hostcell(frame, cell, parameters.minimumHostElevationM))
                hosts.add(cell, area);
        }
    }
    const double cellKm = frame.cellKm();
    if (hosts.total() <= 0.0 || cellKm <= 0.0)
        return deposits;

    detail::Random random(seed);
    const std::uint64_t provinces = random.poisson(parameters.provincesPerMkm2 * 1.0e-6 * hosts.total());
    for (std::uint64_t province = 0; province < provinces; ++province)
    {
        const std::size_t centre = hosts.pick(random.uniform());
        const double centrex = static_cast<double>(centre % frame.width) + random.uniform();
        const double centrey = static_cast<double>(centre / frame.width) + random.uniform();
        const double angle = random.uniform(0.0, kPi);
        const double halflength = random.uniform(parameters.provinceHalfLengthKm.minimum, parameters.provinceHalfLengthKm.maximum) / cellKm;
        const double halfwidth = random.uniform(parameters.provinceHalfWidthKm.minimum, parameters.provinceHalfWidthKm.maximum) / cellKm;
        const std::uint64_t count = random.poisson(parameters.depositsPerProvince);
        ++provinceCount;

        for (std::uint64_t k = 0; k < count; ++k)
        {
            for (int attempt = 0; attempt < kInheritedPlacementAttempts; ++attempt)
            {
                const double radius = std::sqrt(random.uniform());
                const double phi = 2.0 * kPi * random.uniform();
                const double along = halflength * radius * std::cos(phi);
                const double across = halfwidth * radius * std::sin(phi);
                const double x = wrapcoordinate(centrex + along * std::cos(angle) - across * std::sin(angle), frame.width);
                const double y = wrapcoordinate(centrey + along * std::sin(angle) + across * std::cos(angle), frame.height);
                const std::size_t cell = frame.cellIndex(std::min(frame.width - 1, static_cast<std::uint32_t>(x)),
                    std::min(frame.height - 1, static_cast<std::uint32_t>(y)));
                if (!hostcell(frame, cell, parameters.minimumHostElevationM))
                    continue;

                NewDeposit item;
                item.x = x;
                item.y = y;
                Deposit& deposit = item.deposit;
                deposit.system = DepositSystem::orogenic_gold;
                deposit.origin = DepositOrigin::inherited;
                sampledepositbody(deposit, parameters.orogenic, random.next());
                const double topdepth = random.uniform(parameters.inheritedTopDepthM.minimum, parameters.inheritedTopDepthM.maximum);
                deposit.formationDepthM = std::max(deposit.formationDepthM, topdepth);
                deposit.exhumationM = deposit.formationDepthM - topdepth;
                deposit.formation.timeMyr = frame.timeMyr - random.uniform(parameters.inheritedAgeMyr.minimum, parameters.inheritedAgeMyr.maximum);
                deposit.formation.hostCrustClass = frame.crustClass != nullptr ? frame.crustClass[cell] : crustcode::none;
                deposit.formation.provinceId = provinceCount;
                deposits.push_back(item);
                break;
            }
        }
    }
    return deposits;
}
}
