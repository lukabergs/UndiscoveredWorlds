#include "simulations/geology/mineral_systems/mineral_deposit.hpp"

#include "simulations/geology/mineral_systems/tectonic_frame.hpp"

#include <algorithm>
#include <cstdio>
#include <ostream>
#include <sstream>
#include <string>

namespace mineralsystems
{
namespace
{
const char* crustlabel(std::uint8_t crust)
{
    switch (crust)
    {
    case crustcode::oceanic: return "oceanic crust";
    case crustcode::transitional: return "transitional (arc or thinned) crust";
    case crustcode::continental: return "continental crust";
    default: return "unclassified crust";
    }
}

const char* regimelabel(std::uint8_t regime)
{
    switch (regime)
    {
    case regimecode::convergent_arc: return "subduction-arc";
    case regimecode::continent_collision: return "continent-collision";
    case regimecode::divergent_rift: return "continental-rift";
    case regimecode::transform: return "transform";
    case regimecode::passive_margin: return "passive-margin";
    case regimecode::mid_ocean_ridge: return "mid-ocean-ridge";
    case regimecode::trench_adjacent: return "trench";
    default: return "stable";
    }
}

std::string format(const char* pattern, double a = 0.0, double b = 0.0, double c = 0.0, double d = 0.0)
{
    char buffer[256];
    std::snprintf(buffer, sizeof(buffer), pattern, a, b, c, d);
    return buffer;
}

const Deposit* findsource(const std::vector<Deposit>* deposits, std::uint32_t id)
{
    if (deposits == nullptr)
        return nullptr;
    for (const Deposit& deposit : *deposits)
    {
        if (deposit.id == id)
            return &deposit;
    }
    return nullptr;
}
}

const char* depositsystemlabel(DepositSystem system)
{
    switch (system)
    {
    case DepositSystem::orogenic_gold: return "Orogenic gold";
    case DepositSystem::epithermal_gold: return "Epithermal gold";
    case DepositSystem::porphyry_copper_gold: return "Porphyry copper-gold";
    case DepositSystem::placer_gold: return "Placer gold";
    }
    return "Unknown deposit";
}

const char* depositoriginlabel(DepositOrigin origin)
{
    switch (origin)
    {
    case DepositOrigin::simulated: return "formed during the tectonic run";
    case DepositOrigin::inherited: return "inherited from before the run";
    case DepositOrigin::redistributed: return "concentrated by rivers";
    }
    return "unknown origin";
}

const char* depositfatelabel(DepositFate fate)
{
    switch (fate)
    {
    case DepositFate::preserved: return "preserved";
    case DepositFate::eroded: return "eroded away";
    case DepositFate::subducted: return "subducted";
    }
    return "unknown";
}

const char* exposurelabel(Exposure exposure)
{
    switch (exposure)
    {
    case Exposure::outcropping: return "outcropping";
    case Exposure::shallow: return "shallow";
    case Exposure::buried: return "buried";
    case Exposure::deep: return "deep";
    }
    return "unknown";
}

Exposure exposureclass(double topDepthM)
{
    if (topDepthM <= 0.0)
        return Exposure::outcropping;
    if (topDepthM <= 100.0)
        return Exposure::shallow;
    if (topDepthM <= 1000.0)
        return Exposure::buried;
    return Exposure::deep;
}

std::string describedeposit(const Deposit& deposit, const std::vector<Deposit>* deposits)
{
    std::string text = std::string(depositsystemlabel(deposit.system)) + " #" + std::to_string(deposit.id) + ", "
        + depositoriginlabel(deposit.origin) + ".\n";

    if (deposit.system == DepositSystem::placer_gold)
    {
        text += format("%.2f t gold in %.2f Mt of gravel at %.2f g/t, under about %.0f m of overburden.\n",
            deposit.containedGoldT(), deposit.oreMt, deposit.goldGradeGt, deposit.formationDepthM);
        text += "Gold liberated from eroding lodes and trapped here by the present rivers:";
        double listed = 0.0;
        for (const SourceShare& share : deposit.sources)
        {
            const Deposit* source = findsource(deposits, share.depositId);
            text += format(" %.0f%% from #", 100.0 * share.goldTonnes / std::max(1.0e-12, deposit.containedGoldT()))
                + std::to_string(share.depositId);
            if (source != nullptr)
                text += std::string(" (") + depositsystemlabel(source->system) + ")";
            text += ";";
            listed += share.goldTonnes;
        }
        if (deposit.containedGoldT() - listed > 1.0e-9)
            text += format(" %.0f%% from minor sources.", 100.0 * (deposit.containedGoldT() - listed) / deposit.containedGoldT());
        text += "\n";
        return text;
    }

    text += format("%.1f t gold in %.2f Mt ore at %.2f g/t", deposit.containedGoldT(), deposit.oreMt, deposit.goldGradeGt);
    if (deposit.copperGradePct > 0.0)
        text += format(", with %.0f kt copper at %.2f%%", deposit.containedCopperT() / 1000.0, deposit.copperGradePct);
    text += format("; %.0f%% of the ore remains.\n", 100.0 * (1.0 - deposit.erodedFraction));

    if (deposit.origin == DepositOrigin::inherited)
    {
        text += format("Part of inherited province %.0f, dated %.0f Myr before the run began, in ",
            deposit.formation.provinceId, -deposit.formation.timeMyr);
        text += std::string(crustlabel(deposit.formation.hostCrustClass)) + ".\n";
    }
    else
    {
        text += format("Formed %.1f Myr into the run in ", deposit.formation.timeMyr);
        text += std::string(crustlabel(deposit.formation.hostCrustClass));
        text += format(", %.0f km from a %.0f Myr old ", deposit.formation.boundaryDistanceKm, deposit.formation.boundaryAgeMyr);
        text += std::string(regimelabel(deposit.formation.geologicRegime)) + " boundary (segment "
            + std::to_string(deposit.formation.boundarySegmentId) + ")"
            + format("; setting intensity %.2f.\n", deposit.formation.settingIntensity);
    }

    text += format("Formed with its top at %.0f m depth; %.0f m of rock has since been removed above it.\n",
        deposit.formationDepthM, deposit.exhumationM);
    if (deposit.fate == DepositFate::preserved)
    {
        const double top = deposit.topDepthM();
        text += std::string("Now ") + exposurelabel(exposureclass(top))
            + (top > 0.0 ? format(": top at %.0f m depth.\n", top) : std::string(": erosion is cutting into the ore body.\n"));
    }
    else
        text += std::string("Fate: ") + depositfatelabel(deposit.fate) + format(" at %.1f Myr.\n", deposit.fateTimeMyr);

    if (deposit.recentLiberatedGoldT > 0.0)
        text += format("%.2f t of recently liberated gold is available to rivers.\n", deposit.recentLiberatedGoldT);
    if (deposit.surface.evaluated)
    {
        text += deposit.surface.submerged ? std::string("Present surface: submerged")
            : format("Present surface: land at %.0f m, local relief %.0f m, %.0f mm/yr rain",
                deposit.surface.elevationM, deposit.surface.localReliefM, deposit.surface.annualRainfallMm);
        text += deposit.surface.iceCovered ? ", ice covered.\n" : ".\n";
    }
    return text;
}
void writedepositscsv(std::ostream& out, const std::vector<Deposit>& deposits, bool cells)
{
    const std::streamsize precision = out.precision(10);
    out << "id,system,origin,fate,fate_myr,u,v,formed_myr,regime,boundary_type,host_crust,boundary_segment,"
           "province,boundary_km,boundary_age_myr,intensity,ore_mt,gold_g_t,copper_pct,gold_t,remaining_gold_t,"
           "formation_depth_m,extent_m,exhumation_m,top_depth_m,eroded_fraction,recent_liberated_t"
        << (cells ? ",cell\n" : "\n");
    for (const Deposit& d : deposits)
    {
        const FormationRecord& f = d.formation;
        out << d.id << ',' << static_cast<int>(d.system) << ',' << static_cast<int>(d.origin) << ','
            << static_cast<int>(d.fate) << ',' << d.fateTimeMyr << ',' << d.u << ',' << d.v << ',' << f.timeMyr << ','
            << static_cast<int>(f.geologicRegime) << ',' << static_cast<int>(f.boundaryType) << ','
            << static_cast<int>(f.hostCrustClass) << ',' << f.boundarySegmentId << ',' << f.provinceId << ','
            << f.boundaryDistanceKm << ',' << f.boundaryAgeMyr << ',' << f.settingIntensity << ',' << d.oreMt << ','
            << d.goldGradeGt << ',' << d.copperGradePct << ',' << d.containedGoldT() << ',' << d.remainingGoldT() << ','
            << d.formationDepthM << ',' << d.verticalExtentM << ',' << d.exhumationM << ',' << d.topDepthM() << ','
            << d.erodedFraction << ',' << d.recentLiberatedGoldT;
        if (cells)
            out << ',' << d.surface.cell;
        out << '\n';
    }
    out.precision(precision);
}

std::string mineralbudgetjson(const MineralBudget& budget)
{
    std::ostringstream out;
    out.precision(12);
    out << "{\n  \"schema\": \"mineral-budget/v1\",\n  \"formed\": {";
    for (int system = 0; system < kDepositSystemCount; ++system)
    {
        out << (system == 0 ? "" : ",") << "\n    \"" << depositsystemlabel(static_cast<DepositSystem>(system))
            << "\": {\"count\": " << budget.formedCount[system] << ", \"gold_t\": " << budget.formedGoldT[system] << "}";
    }
    out << "\n  },\n"
        << "  \"inherited\": {\"count\": " << budget.inheritedCount << ", \"gold_t\": " << budget.inheritedGoldT << "},\n"
        << "  \"eroded\": {\"count\": " << budget.erodedCount << ", \"gold_t\": " << budget.erodedGoldT << "},\n"
        << "  \"subducted\": {\"count\": " << budget.subductedCount << ", \"gold_t\": " << budget.subductedGoldT << "},\n"
        << "  \"placer\": {\"injected_t\": " << budget.placerInjectedGoldT << ", \"deposit_t\": " << budget.placerDepositGoldT
        << ", \"dispersed_t\": " << budget.placerDispersedGoldT << ", \"lake_t\": " << budget.placerLakeGoldT
        << ", \"exported_t\": " << budget.placerExportedGoldT << "}\n}\n";
    return out.str();
}
}
