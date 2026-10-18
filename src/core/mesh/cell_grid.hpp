#pragma once

// The table-driven cell grid of the spherical world mesh
// (docs/architecture/sphere-mesh.md, sections 2 and 3). One class, two kinds:
// the equiangular cubed sphere (orientation uw-cube-1, ids f N^2 + j N + i) and
// the cell-centred W x W/2 latitude-longitude grid (row-major, row 0 north,
// longitudes east from -180), whose tables reproduce landscape::SphereMetric
// exactly. Neighbour slots are the D8 codes 1 N, 2 NE, 3 E, 4 SE, 5 S, 6 SW,
// 7 W, 8 NW in the cell's own frame (geographic on lat-lon; toward v and u on a
// cube face, also across face edges); -1 where no neighbour exists. A subgrid
// (a window) is a grid of the same kind over a subset of cells, local ids in
// ascending global order, neighbours outside the subset -1.

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace mesh
{
using CellId = std::uint32_t;
constexpr CellId nocell = 0xffffffffu;  // locate() outside a subgrid
constexpr int slotcount = 8;            // neighbour slots 1..8; 0 is none

struct Vec3
{
    double x = 0.0, y = 0.0, z = 0.0;
};

inline Vec3 operator+(const Vec3& p, const Vec3& q) { return { p.x + q.x, p.y + q.y, p.z + q.z }; }
inline Vec3 operator-(const Vec3& p, const Vec3& q) { return { p.x - q.x, p.y - q.y, p.z - q.z }; }
inline Vec3 operator*(double s, const Vec3& p) { return { s * p.x, s * p.y, s * p.z }; }
inline double dot(const Vec3& p, const Vec3& q) { return p.x * q.x + p.y * q.y + p.z * q.z; }
inline Vec3 cross(const Vec3& p, const Vec3& q)
{
    return { p.y * q.z - p.z * q.y, p.z * q.x - p.x * q.z, p.x * q.y - p.y * q.x };
}
inline double length(const Vec3& p) { return std::sqrt(dot(p, p)); }
inline Vec3 normalise(const Vec3& p) { return (1.0 / length(p)) * p; }

// The unit vector of a latitude and longitude (radians): (cos lat cos lon, cos lat sin lon, sin lat).
Vec3 unitvector(double latitude_rad, double longitude_rad);
double latitudeof(const Vec3& p);   // atan2(z, hypot(x, y))
double longitudeof(const Vec3& p);  // atan2(y, x), 0 at a pole
// Local unit tangents at a latitude and longitude: east (-sin lon, cos lon, 0),
// north (-sin lat cos lon, -sin lat sin lon, cos lat). At a pole (longitude 0)
// east is +Y and north is the limit along longitude 0 (-X at the north pole).
Vec3 eastof(double latitude_rad, double longitude_rad);
Vec3 northof(double latitude_rad, double longitude_rad);

enum class GridKind
{
    cube,
    latlon,
};

class CellGrid
{
public:
    GridKind kind() const { return kind_; }
    int facecells() const { return facecells_; }  // cube N (0 on lat-lon)
    int columns() const { return columns_; }      // lat-lon W (0 on the cube)
    int rows() const { return rows_; }            // lat-lon W / 2 (0 on the cube)
    std::size_t cellcount() const { return cellcount_; }
    double radius_m() const { return radius_m_; }

    Vec3 centre(CellId c) const;
    double latitude_rad(CellId c) const;
    double longitude_rad(CellId c) const;
    Vec3 east(CellId c) const { return eastof(latitude_rad(c), longitude_rad(c)); }
    Vec3 north(CellId c) const { return northof(latitude_rad(c), longitude_rad(c)); }

    double area_m2(CellId c) const;
    double totalarea_m2() const { return totalarea_m2_; }

    // Slot k in 1..8 (anything else: -1, 0, 0, 0, 0).
    std::int32_t neighbour(CellId c, int k) const;
    int reverse(CellId c, int k) const;  // the slot by which neighbour(c, k) reaches c; 0 when none
    double distance_m(CellId c, int k) const;
    double facelength_m(CellId c, int k) const;
    double bearing_rad(CellId c, int k) const;  // clockwise from geographic north, [0, 2 pi)

    // The cell containing a unit vector; on a subgrid nocell when the cell is outside it.
    CellId locate(const Vec3& p) const;

    // Lat-lon only (row-major, row 0 north). On a subgrid x, y, column and row
    // are those of the full grid; cellof returns nocell outside the subset.
    CellId cellof(int x, int y) const;
    int column(CellId c) const;
    int row(CellId c) const;

    // Windows: the global id of a local id (the id itself on a full grid).
    bool issubgrid() const { return subgrid_; }
    CellId globalid(CellId local) const { return subgrid_ ? globalids_[local] : local; }
    // The local id of a global id, nocell when outside (the id itself on a full grid).
    CellId localid(CellId global) const;

    // The smallest centre spacing over every cell and cardinal slot (cellsalongarc's step).
    double minimumspacing_m() const { return minimumspacing_m_; }

private:
    friend std::shared_ptr<const CellGrid> makecubegrid(int, double);
    friend std::shared_ptr<const CellGrid> makelatlongrid(int, double);
    friend std::shared_ptr<const CellGrid> subgrid(const CellGrid&, const std::vector<CellId>&);

    CellId locateglobal(const Vec3& p) const;

    GridKind kind_ = GridKind::latlon;
    int facecells_ = 0;
    int columns_ = 0;
    int rows_ = 0;
    std::size_t cellcount_ = 0;
    double radius_m_ = 0.0;
    double totalarea_m2_ = 0.0;
    double minimumspacing_m_ = 0.0;

    // Stored tables (a cube grid, any subgrid), per cell or per cell and slot (8 per cell).
    bool stored_ = false;
    std::vector<Vec3> centre_;
    std::vector<double> latitude_;
    std::vector<double> longitude_;
    std::vector<double> area_;
    std::vector<std::int32_t> neighbour_;
    std::vector<std::int8_t> reverse_;
    std::vector<double> distance_;
    std::vector<double> facelength_;
    std::vector<double> bearing_;

    // Analytic lat-lon tables (a full lat-lon grid), per row and per column.
    double dy_m_ = 0.0;
    std::vector<double> rowarea_;
    std::vector<double> rowdx_;
    std::vector<double> rowlatitude_;
    std::vector<double> columnlongitude_;

    // A subgrid: the ascending global ids of its cells.
    bool subgrid_ = false;
    std::vector<CellId> globalids_;
};

// The equiangular cubed sphere of N cells per face edge (N >= 1).
std::shared_ptr<const CellGrid> makecubegrid(int facecells, double radius_m);
// The cell-centred W x W/2 grid (W even, >= 2).
std::shared_ptr<const CellGrid> makelatlongrid(int columns, double radius_m);
// A grid over a subset of a full grid's cells (sorted and deduplicated);
// throws std::invalid_argument for a subgrid of a subgrid or an id out of range.
std::shared_ptr<const CellGrid> subgrid(const CellGrid& grid, const std::vector<CellId>& globalids);
}
