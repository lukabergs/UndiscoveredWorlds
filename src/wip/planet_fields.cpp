#include "planet_fields.hpp"

#include <array>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#include "climate_world_fields.hpp"
#include "planet.hpp"

using worldmodel::ElementType;
using worldmodel::FieldId;
using worldmodel::TableId;

std::uint64_t fnv1a64(const void* data, std::size_t size, std::uint64_t hash)
{
    const auto* bytes = static_cast<const unsigned char*>(data);
    for (std::size_t index = 0; index < size; index++)
    {
        hash ^= bytes[index];
        hash *= 1099511628211ULL;
    }
    return hash;
}

namespace
{
std::string csvtext(const std::string& value)
{
    if (value.find_first_of(",\"\n\r") == std::string::npos)
        return value;
    std::string quoted = "\"";
    for (char character : value)
    {
        if (character == '"')
            quoted += '"';
        quoted += character;
    }
    return quoted + "\"";
}
}

// Friend of planet: maps each catalogue id onto the current storage. Fields
// owned by the climate module are read through climate_world_fields.
struct planetfieldaccess
{
    static PlanetFieldData read(const planet& world, FieldId id)
    {
        namespace encode = worldmodel::encoding;
        PlanetFieldData data;
        if (readclimateworldfield(world.climate_, id, data, false))
            return data;
        // Absent diagnostics read as empty data: bundles skip them, hashes see them appear.
        if (isclimatediagnosticfield(id))
        {
            data.id = id;
            return data;
        }

        const std::size_t columns = static_cast<std::size_t>(world.body_.itswidth) + 1;
        const std::size_t rows = static_cast<std::size_t>(world.body_.itsheight) + 1;
        switch (id)
        {
        case FieldId::seed: return encode::scalar(id, world.body_.itsseed);
        case FieldId::size_class: return encode::scalar(id, world.body_.itssize);
        case FieldId::width_index: return encode::scalar(id, world.body_.itswidth);
        case FieldId::height_index: return encode::scalar(id, world.body_.itsheight);
        case FieldId::prograde_rotation: return encode::scalar(id, world.body_.itsrotation);
        case FieldId::axial_tilt: return encode::scalar(id, world.body_.itstilt);
        case FieldId::eccentricity: return encode::scalar(id, world.body_.itseccentricity);
        case FieldId::perihelion: return encode::scalar(id, world.body_.itsperihelion);
        case FieldId::gravity: return encode::scalar(id, world.body_.itsgravity);
        case FieldId::lunar_strength: return encode::scalar(id, world.body_.itslunar);
        case FieldId::lapse_rate: return encode::scalar(id, world.atmosphere_.itstempdecrease);
        case FieldId::north_polar_adjust: return encode::scalar(id, world.atmosphere_.itsnorthpolaradjust);
        case FieldId::south_polar_adjust: return encode::scalar(id, world.atmosphere_.itssouthpolaradjust);
        case FieldId::average_temperature: return encode::scalar(id, world.atmosphere_.itsaveragetemp);
        case FieldId::water_pickup: return encode::scalar(id, world.atmosphere_.itswaterpickup);
        case FieldId::river_factor: return encode::scalar(id, world.hydrosphere_.itsriverfactor);
        case FieldId::river_land_reduce: return encode::scalar(id, world.hydrosphere_.itsriverlandreduce);
        case FieldId::estuary_limit: return encode::scalar(id, world.hydrosphere_.itsestuarylimit);
        case FieldId::glacial_temperature: return encode::scalar(id, world.cryosphere_.itsglacialtemp);
        case FieldId::glacier_temperature: return encode::scalar(id, world.cryosphere_.itsglaciertemp);
        case FieldId::mountain_reduce: return encode::scalar(id, world.geosphere_.itsmountainreduce);
        case FieldId::climate_type_count: return encode::scalar(id, world.atmosphere_.itsclimateno);
        case FieldId::max_elevation: return encode::scalar(id, world.geosphere_.itsmaxheight);
        case FieldId::sea_level: return encode::scalar(id, world.geosphere_.itssealevel);
        case FieldId::max_river_flow: return encode::scalar(id, world.hydrosphere_.itsmaxriverflow);
        case FieldId::land_cells: return encode::scalar(id, world.geosphere_.itslandtotal);
        case FieldId::sea_cells: return encode::scalar(id, world.geosphere_.itsseatotal);
        case FieldId::crater_count: return encode::scalar(id, world.geosphere_.itscraterno);
        case FieldId::tectonic_time_origin_step: return encode::scalar(id, world.geosphere_.itstectonictimeoriginstep);
        case FieldId::tectonic_time: return encode::scalar(id, world.geosphere_.itstectonictimemyr);
        case FieldId::tectonic_time_step: return encode::scalar(id, world.geosphere_.itstectonicdeltatimemyr);
        case FieldId::tectonic_plate_count: return encode::scalar(id, world.geosphere_.itstectonicplatecount);
        case FieldId::tectonic_sea_level: return encode::scalar(id, world.geosphere_.itstectonicsealevelm);
        case FieldId::base_elevation: return encode::grid(id, columns, rows, world.geosphere_.mapnom);
        case FieldId::continental_shelf: return encode::grid(id, columns, rows, world.geosphere_.shelfmap);
        case FieldId::mountain_height: return encode::grid(id, columns, rows, world.geosphere_.mountainheights);
        case FieldId::mountain_ridges: return encode::grid(id, columns, rows, world.geosphere_.mountainridges);
        case FieldId::extra_elevation: return encode::grid(id, columns, rows, world.geosphere_.extraelevmap);
        case FieldId::crater_rims: return encode::grid(id, columns, rows, world.geosphere_.craterrims);
        case FieldId::crater_centres: return encode::grid(id, columns, rows, world.geosphere_.cratercentres);
        case FieldId::volcano: return encode::grid(id, columns, rows, world.geosphere_.volcanomap);
        case FieldId::stratovolcano: return encode::grid(id, columns, rows, world.geosphere_.stratomap);
        case FieldId::island: return encode::grid(id, columns, rows, world.geosphere_.islandmap);
        case FieldId::mountain_island: return encode::grid(id, columns, rows, world.geosphere_.mountainislandmap);
        case FieldId::ocean_ridge: return encode::grid(id, columns, rows, world.geosphere_.oceanridgemap);
        case FieldId::ocean_ridge_height: return encode::grid(id, columns, rows, world.geosphere_.oceanridgeheightmap);
        case FieldId::ocean_rift: return encode::grid(id, columns, rows, world.geosphere_.oceanriftmap);
        case FieldId::ocean_ridge_offset: return encode::grid(id, columns, rows, world.geosphere_.oceanridgeoffsetmap);
        case FieldId::ocean_ridge_angle: return encode::grid(id, columns, rows, world.geosphere_.oceanridgeanglemap);
        case FieldId::roughness: return encode::grid(id, columns, rows, world.geosphere_.roughnessmap);
        case FieldId::detail_noise: return encode::grid(id, columns, rows, world.geosphere_.noisemap);
        case FieldId::geologic_regime: return encode::grid(id, columns, rows, world.geosphere_.geologicregimemap);
        case FieldId::tectonic_convergence: return encode::grid(id, columns, rows, world.geosphere_.tectonicconvergencemap);
        case FieldId::tectonic_divergence: return encode::grid(id, columns, rows, world.geosphere_.tectonicdivergencemap);
        case FieldId::tectonic_shear: return encode::grid(id, columns, rows, world.geosphere_.tectonicshearmap);
        case FieldId::crust_age: return encode::grid(id, columns, rows, world.geosphere_.tectoniccrustagemyrmap);
        case FieldId::crust_thickness: return encode::grid(id, columns, rows, world.geosphere_.tectoniccrustthicknessmap);
        case FieldId::crust_class: return encode::grid(id, columns, rows, world.geosphere_.tectoniccrustclassmap);
        case FieldId::uplift_tendency: return encode::grid(id, columns, rows, world.geosphere_.tectonicuplifttendencymap);
        case FieldId::subsidence_tendency: return encode::grid(id, columns, rows, world.geosphere_.tectonicsubsidencetendencymap);
        case FieldId::accumulated_strain: return encode::grid(id, columns, rows, world.geosphere_.tectonicaccumulatedstrainmap);
        case FieldId::boundary_type: return encode::grid(id, columns, rows, world.geosphere_.tectonicboundarytypemap);
        case FieldId::boundary_distance: return encode::grid(id, columns, rows, world.geosphere_.tectonicboundarydistancemap);
        case FieldId::boundary_segment: return encode::grid(id, columns, rows, world.geosphere_.tectonicboundarysegmentidmap);
        case FieldId::nearest_boundary: return encode::grid(id, columns, rows, world.geosphere_.tectonicnearestboundaryidmap);
        case FieldId::boundary_history: return encode::grid(id, columns, rows, world.geosphere_.tectonicboundaryhistorymap);
        case FieldId::deforming_region: return encode::grid(id, columns, rows, world.geosphere_.tectonicdeformingregionidmap);
        case FieldId::deforming_region_type: return encode::grid(id, columns, rows, world.geosphere_.tectonicdeformingregiontypemap);
        case FieldId::deformation_rate: return encode::grid(id, columns, rows, world.geosphere_.tectonicdeformationratemap);
        case FieldId::deformation_velocity_x: return encode::grid(id, columns, rows, world.geosphere_.tectonicdeformationvelocityxmap);
        case FieldId::deformation_velocity_y: return encode::grid(id, columns, rows, world.geosphere_.tectonicdeformationvelocityymap);
        case FieldId::erosion_potential: return encode::grid(id, columns, rows, world.geosphere_.erosionpotentialmap);
        case FieldId::deposition_potential: return encode::grid(id, columns, rows, world.geosphere_.depositionpotentialmap);
        case FieldId::floodplain_fertility: return encode::grid(id, columns, rows, world.geosphere_.floodplainfertilitymap);
        case FieldId::metal_ore_potential: return encode::grid(id, columns, rows, world.geosphere_.metalorereservemap);
        case FieldId::placer_potential: return encode::grid(id, columns, rows, world.geosphere_.placerreservemap);
        case FieldId::evaporite_potential: return encode::grid(id, columns, rows, world.geosphere_.evaporitereservemap);
        case FieldId::volcanic_potential: return encode::grid(id, columns, rows, world.geosphere_.volcanicreservemap);
        case FieldId::river_direction: return encode::grid(id, columns, rows, world.hydrosphere_.rivermapdir);
        case FieldId::river_flow_january: return encode::grid(id, columns, rows, world.hydrosphere_.rivermapjan);
        case FieldId::river_flow_july: return encode::grid(id, columns, rows, world.hydrosphere_.rivermapjul);
        case FieldId::lake_surface: return encode::grid(id, columns, rows, world.hydrosphere_.lakemap);
        case FieldId::rift_lake_surface: return encode::grid(id, columns, rows, world.hydrosphere_.riftlakemapsurface);
        case FieldId::rift_lake_bed: return encode::grid(id, columns, rows, world.hydrosphere_.riftlakemapbed);
        case FieldId::lake_start: return encode::grid(id, columns, rows, world.hydrosphere_.lakestartmap);
        case FieldId::special_features: return encode::grid(id, columns, rows, world.hydrosphere_.specials);
        case FieldId::tide: return encode::grid(id, columns, rows, world.hydrosphere_.tidalmap);
        case FieldId::basin_class: return encode::grid(id, columns, rows, world.hydrosphere_.basinclassmap);
        case FieldId::delta_direction: return encode::grid(id, columns, rows, world.hydrosphere_.deltamapdir);
        case FieldId::delta_flow_january: return encode::grid(id, columns, rows, world.hydrosphere_.deltamapjan);
        case FieldId::delta_flow_july: return encode::grid(id, columns, rows, world.hydrosphere_.deltamapjul);
        case FieldId::fishery_potential: return encode::grid(id, columns, rows, world.biosphere_.fisheryreservemap);
        case FieldId::settlement_suitability: return encode::grid(id, columns, rows, world.anthroposphere_.settlement_suitability);
        case FieldId::infrastructure: return encode::grid(id, columns, rows, world.anthroposphere_.social_infrastructure);
        case FieldId::agricultural_capacity: return encode::grid(id, columns, rows, world.anthroposphere_.agricultural_capacity);
        case FieldId::route_traffic: return encode::grid(id, columns, rows, world.anthroposphere_.route_traffic);
        case FieldId::river_access: return encode::grid(id, columns, rows, world.anthroposphere_.river_access);
        case FieldId::harbor_score: return encode::grid(id, columns, rows, world.anthroposphere_.harbor_score);
        case FieldId::owner_settlement: return encode::grid(id, columns, rows, world.anthroposphere_.owner_settlement_id);
        case FieldId::owner_polity: return encode::grid(id, columns, rows, world.anthroposphere_.owner_polity_id);
        case FieldId::no_shade: return encode::grid(id, columns, rows, world.presentation_.noshademap);
        case FieldId::debug_test: return encode::grid(id, columns, rows, world.presentation_.testmap);
        default:
            throw std::logic_error("World field has no storage: " + worldmodel::qualifiedname(id));
        }
    }

    static std::string table(const planet& world, TableId id, std::string& extension)
    {
        std::ostringstream out;
        out << std::setprecision(9);
        extension = ".csv";
        switch (id)
        {
        case TableId::tectonic_boundary_segments:
            out << "id,left_plate,right_plate,cells,persistence_steps,centroid_x,centroid_y,length_cells,normal_motion,shear_motion,"
                   "convergence_score,divergence_score,shear_score,boundary_type,geologic_regime,age_myr\n";
            for (const TectonicBoundarySegment& s : world.geosphere_.tectonicboundarysegmentlist)
                out << s.id << ',' << s.leftPlateId << ',' << s.rightPlateId << ',' << s.cellCount << ',' << s.persistenceSteps << ','
                    << s.centroidX << ',' << s.centroidY << ',' << s.lengthCells << ',' << s.averageNormalMotion << ','
                    << s.averageShearMotion << ',' << s.averageConvergenceScore << ',' << s.averageDivergenceScore << ','
                    << s.averageShearScore << ',' << static_cast<int>(s.boundaryType) << ',' << static_cast<int>(s.geologicRegime) << ','
                    << s.ageMyr << '\n';
            break;
        case TableId::tectonic_deforming_regions:
            out << "id,boundary_segment,primary_plate,secondary_plate,cells,persistence_steps,centroid_x,centroid_y,deformation_rate,"
                   "velocity_x,velocity_y,normal_motion,shear_motion,type,age_myr\n";
            for (const TectonicDeformingRegion& r : world.geosphere_.tectonicdeformingregionlist)
                out << r.id << ',' << r.boundarySegmentId << ',' << r.primaryPlateId << ',' << r.secondaryPlateId << ',' << r.cellCount << ','
                    << r.persistenceSteps << ',' << r.centroidX << ',' << r.centroidY << ',' << r.averageDeformationRate << ','
                    << r.averageInterpolatedVelocityX << ',' << r.averageInterpolatedVelocityY << ',' << r.averageNormalMotion << ','
                    << r.averageShearMotion << ',' << static_cast<int>(r.type) << ',' << r.ageMyr << '\n';
            break;
        case TableId::mineral_deposits:
            mineralsystems::writedepositscsv(out, world.geosphere_.mineraldepositlist);
            break;
        case TableId::mineral_budget:
            extension = ".json";
            out << mineralsystems::mineralbudgetjson(world.geosphere_.mineralbudgetstate);
            break;
        case TableId::craters:
            out << "w,x,y,z\n";
            for (int index = 0; index < world.geosphere_.itscraterno && index < MAXCRATERS; index++)
            {
                const fourshorts& crater = world.geosphere_.cratercentreslist[index];
                out << crater.w << ',' << crater.x << ',' << crater.y << ',' << crater.z << '\n';
            }
            break;
        case TableId::settlements:
            out << "id,x,y,name,urban_population,rural_population,carrying_capacity,infrastructure,market_strength,harbor,river_access,polity\n";
            for (const Settlement& s : world.anthroposphere_.settlementlist)
                out << s.id << ',' << s.x << ',' << s.y << ',' << csvtext(s.name) << ',' << s.urbanPopulation << ',' << s.ruralPopulation << ','
                    << s.carryingCapacity << ',' << s.infrastructure << ',' << s.marketStrength << ',' << s.harbor << ','
                    << s.riverAccess << ',' << s.polityId << '\n';
            break;
        case TableId::polities:
            out << "id,name,capital,cohesion,population,infrastructure,wealth,military_pressure\n";
            for (const Polity& p : world.anthroposphere_.politylist)
                out << p.id << ',' << csvtext(p.name) << ',' << p.capitalSettlementId << ',' << p.cohesion << ',' << p.population << ','
                    << p.infrastructure << ',' << p.wealth << ',' << p.militaryPressure << '\n';
            break;
        case TableId::route_edges:
            out << "from,to,mode,cost,capacity,traffic\n";
            for (const RouteEdge& e : world.anthroposphere_.routeedgelist)
                out << e.fromSettlementId << ',' << e.toSettlementId << ',' << static_cast<int>(e.mode) << ',' << e.cost << ','
                    << e.capacity << ',' << e.traffic << '\n';
            break;
        case TableId::history_events:
            out << "year,type,primary_settlement,primary_polity,secondary_settlement,secondary_polity,summary\n";
            for (const HistoryEvent& e : world.anthroposphere_.historyeventlist)
                out << e.year << ',' << csvtext(e.type) << ',' << e.primarySettlementId << ',' << e.primaryPolityId << ','
                    << e.secondarySettlementId << ',' << e.secondaryPolityId << ',' << csvtext(e.summary) << '\n';
            break;
        case TableId::weather_state:
        case TableId::atmosphere_levels:
        {
            std::string text;
            readclimateworldtable(world.climate_, id, text, extension);
            out << text;
            break;
        }
        }
        return out.str();
    }
};

PlanetFieldData readplanetfield(const planet& world, FieldId id)
{
    return planetfieldaccess::read(world, id);
}

std::uint64_t hashplanetfield(const planet& world, FieldId id)
{
    const PlanetFieldData data = readplanetfield(world, id);
    std::uint64_t hash = fnv1a64(data.shape.data(), data.shape.size() * sizeof(std::size_t));
    return fnv1a64(data.bytes.data(), data.bytes.size(), hash);
}

std::string planettabletext(const planet& world, TableId id, std::string* extension)
{
    std::string suffix;
    std::string text = planetfieldaccess::table(world, id, suffix);
    if (extension != nullptr)
        *extension = suffix;
    return text;
}

std::uint64_t hashplanettable(const planet& world, TableId id)
{
    const std::string text = planettabletext(world, id);
    return fnv1a64(text.data(), text.size());
}
