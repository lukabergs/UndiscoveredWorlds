#pragma once

#include "simulations/geology/mineral_systems/mineral_deposit.hpp"
#include "simulations/geology/mineral_systems/tectonic_frame.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

// Formation of primary gold deposits from the tectonic setting of each native
// cell. Settings are resolution independent: they depend on the distance in km
// to the nearest boundary segment and on that segment's type, activity and age,
// and rates are per unit area. Host crust stands near sea level or above
// (continents, arcs and accreted terranes); deep ocean floor hosts no deposits
// in this model.
namespace mineralsystems
{
struct LognormalRange
{
    double median = 0.0; // Zero disables the quantity.
    double sigmaLn = 0.0;
};

struct UniformRange
{
    double minimum = 0.0;
    double maximum = 0.0;
};

// Size and depth assumptions for one mineral system. The grade-tonnage medians
// are order-of-magnitude choices after USGS descriptive and grade-tonnage models
// (Cox and Singer 1986; Singer, Berger and Moring 2008). They are calibration
// inputs, not validated deposit statistics.
struct SystemParameters
{
    double ratePerMkm2Myr = 0.0; // Expected deposits per 10^6 km^2 per Myr at unit intensity.
    LognormalRange oreMt;
    LognormalRange goldGradeGt;
    LognormalRange copperGradePct;
    UniformRange formationDepthM; // Depth of the ore-body top.
    UniformRange verticalExtentM;
};

struct FormationParameters
{
    // Mesozonal lode gold: large vertical continuity, formed at mid-crustal depth.
    SystemParameters orogenic{ 2.0, { 2.0, 1.6 }, { 5.0, 0.6 }, {}, { 3000.0, 15000.0 }, { 800.0, 2500.0 } };
    // Low-sulphidation veins in the upper 1 km above arc magmatism.
    SystemParameters epithermal{ 3.0, { 1.0, 1.5 }, { 5.0, 0.8 }, {}, { 100.0, 1000.0 }, { 200.0, 700.0 } };
    // Bulk-tonnage copper with by-product gold, emplaced at 1.5-4 km.
    SystemParameters porphyry{ 0.5, { 200.0, 1.2 }, { 0.3, 0.7 }, { 0.45, 0.35 }, { 1500.0, 4000.0 }, { 800.0, 2000.0 } };

    double orogenicWidthKm = 300.0; // Gaussian width of the orogenic belt from the boundary.
    double arcPeakKm = 150.0;       // Arc-trench gap.
    double arcWidthKm = 100.0;
    double collisionWeight = 1.0;
    double accretionaryWeight = 0.4; // Overriding plate of a subduction margin.
    double transformWeight = 0.25;
    // Orogenic fluids follow sustained deformation; boundary age sets maturity.
    double orogenicOnsetMyr = 5.0;
    double orogenicMatureMyr = 20.0;
    // Segment mean convergence (shear for transforms) score. The engine scores
    // 0..100 over 0..5 cm/yr, so full activity is 2.5 cm/yr.
    double fullActivityScore = 50.0;
    // A plate edge whose mean continental crust within one cell of a convergent
    // boundary reaches this is too buoyant to subduct; two such edges collide.
    double buoyantCrustM = 20000.0;
    double minimumHostElevationM = -1000.0;

    // Provinces seeded in the initial host crust (continents and shelves) to
    // represent history before the run, e.g. greenstone belts in cratons.
    double provincesPerMkm2 = 0.4;
    double depositsPerProvince = 25.0;
    UniformRange provinceHalfLengthKm{ 150.0, 400.0 };
    UniformRange provinceHalfWidthKm{ 30.0, 80.0 };
    UniformRange inheritedAgeMyr{ 300.0, 2500.0 };
    // Old provinces are deeply eroded: ore formed at mid-crustal depth now lies
    // within a few km of the surface.
    UniformRange inheritedTopDepthM{ 0.0, 3000.0 };
};

struct SettingIntensity
{
    float orogenic = 0.0f;
    float arc = 0.0f;
};

// Per-frame setting evaluator; kernels are tabulated by integer cell distance.
// Convergent segments are classified from the crust each plate brings to the
// boundary, as in the column solver: plates collide when both edges carry
// buoyant continental crust. Otherwise the continental plate, or else the one
// with younger ocean floor, overrides the subducting plate and alone hosts the
// arc and accretionary settings. The engine's segment regimes rest on coarse
// crust classes and are not used.
class SettingField
{
public:
    SettingField(const TectonicFrame& frame, const FormationParameters& parameters);
    SettingIntensity at(std::size_t cell) const;
    // Contract regime code of the segment's setting: continent collision,
    // convergent arc (any subduction margin), transform, or stable.
    std::uint8_t regime(std::uint32_t boundaryId) const;

private:
    struct BoundarySetting
    {
        float orogenic = 0.0f;
        float arc = 0.0f;
        std::uint32_t plate = kNoPlate; // Sole host plate; kNoPlate for both sides.
        std::uint8_t regime = regimecode::stable;
    };

    const TectonicFrame& frame_;
    double minimumHostElevationM_;
    std::vector<float> orogenicKernel_;
    std::vector<float> arcKernel_;
    std::vector<BoundarySetting> boundaries_;
};

struct NewDeposit
{
    Deposit deposit;
    mesh::Vec3 position{ 1.0, 0.0, 0.0 }; // Unit vector, inside an engine cell.
};

struct FormationSample
{
    std::vector<NewDeposit> deposits;
    // Effective active areas, sum of cell area times intensity.
    double orogenicAreaKm2 = 0.0;
    double arcAreaKm2 = 0.0;
};

// Deposits formed during dtMyr under the frame's settings. Ids are unset.
FormationSample sampleformation(const TectonicFrame& frame, const FormationParameters& parameters, double dtMyr,
    std::uint64_t seed);

// Inherited orogenic provinces in the frame's host crust. Ids are unset.
std::vector<NewDeposit> sampleinheritedprovinces(const TectonicFrame& frame, const FormationParameters& parameters,
    std::uint64_t seed, std::uint32_t& provinceCount);

// Samples size and depth for one deposit of the given system.
void sampledepositbody(Deposit& deposit, const SystemParameters& parameters, std::uint64_t seed);
}
