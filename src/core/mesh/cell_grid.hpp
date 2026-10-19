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
// ascending global order, neighbours outside the subset -1. Windows of a
// refinement level (section 8) are built from 64-bit level keys without the
// whole level: cubewindowgrid (a key set, ascending key order) and
// latlonwindowgrid (a box of a level's lat-lon grid, row-major box order).

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
// A cell's id at a refinement level (sphere-mesh.md section 8): on the cube
// f N_L^2 + j N_L + i with N_L up to 2^30; on lat-lon J W_L + I (row-major, row 0 north).
using LevelKey = std::uint64_t;
constexpr LevelKey nokey = ~LevelKey{ 0 };

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

struct LatLonWindowBox;

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

    // Windows: the global id of a local id (the id itself on a full grid; on a
    // level window the key when it fits 32 bits, else nocell: use globalkey).
    bool issubgrid() const { return subgrid_; }
    CellId globalid(CellId local) const;
    // The local id of a global id, nocell when outside (the id itself on a full grid).
    CellId localid(CellId global) const;
    // The level key of a local id (the global id on a full grid or a subgrid) and
    // its reverse, nocell when outside: a binary search on a cube window, closed
    // form on a lat-lon window.
    LevelKey globalkey(CellId local) const;
    CellId localofkey(LevelKey key) const;
    // A cube window's tables computed on demand (sphere-mesh.md section 8).
    bool iscompact() const { return compact_; }

    // The smallest centre spacing over every cell and cardinal slot (cellsalongarc's step).
    double minimumspacing_m() const { return minimumspacing_m_; }

private:
    friend std::shared_ptr<const CellGrid> makecubegrid(int, double);
    friend std::shared_ptr<const CellGrid> makelatlongrid(int, double);
    friend std::shared_ptr<const CellGrid> subgrid(const CellGrid&, const std::vector<CellId>&);
    friend std::shared_ptr<const CellGrid> cubewindowgrid(int, double, const std::vector<LevelKey>&, bool);
    friend std::shared_ptr<const CellGrid> latlonwindowgrid(int, double, const LatLonWindowBox&);
    friend std::shared_ptr<const CellGrid> makeplanarwindow(int, int, double);

    LevelKey locatekey(const Vec3& p) const;  // the level key of the full grid's cell
    int analyticwidth() const { return latlonwindow_ ? windowcolumns_ : columns_; }

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

    // A cube window: the ascending level keys of its cells; compact keeps only
    // keys, centres, areas, neighbours and reverses and computes the rest.
    bool keywindow_ = false;
    bool compact_ = false;
    std::vector<LevelKey> keys_;

    // A lat-lon window: a box of the level grid (columns_ x rows_), its first
    // level column and row and its size; the analytic tables are per box row
    // and box column.
    bool latlonwindow_ = false;
    std::int64_t windowcolumn0_ = 0;
    int windowrow0_ = 0;
    int windowcolumns_ = 0;
    int windowrows_ = 0;
};

// A box of a level's lat-lon grid in level cells: columns column0 ..
// column0 + columns - 1, each modulo the level's columns (the box may cross the
// date line), rows row0 .. row0 + rows - 1 (no wrap over a pole).
struct LatLonWindowBox
{
    std::int64_t column0 = 0;
    int columns = 0;
    int row0 = 0;
    int rows = 0;
};

// The equiangular cubed sphere of N cells per face edge (N >= 1).
std::shared_ptr<const CellGrid> makecubegrid(int facecells, double radius_m);
// The cell-centred W x W/2 grid (W even, >= 2).
std::shared_ptr<const CellGrid> makelatlongrid(int columns, double radius_m);
// A grid over a subset of a full grid's cells (sorted and deduplicated);
// throws std::invalid_argument for a subgrid of a subgrid or an id out of range.
std::shared_ptr<const CellGrid> subgrid(const CellGrid& grid, const std::vector<CellId>& globalids);

// A window of the cube level with facecells_L cells per face edge (1 .. 2^30)
// over a set of level keys (sorted and deduplicated here), built in closed form
// from the keys: local ids in ascending key order, neighbours outside the set -1,
// distance, face length and bearing those of the full level (the cube's own
// formulas), areas exact. Compact (about 80 bytes per cell) computes latitude,
// longitude, distance, face length and bearing on demand, bit-identical to the
// stored mode. Throws std::invalid_argument for N_L out of range or a key >= 6 N_L^2.
std::shared_ptr<const CellGrid> cubewindowgrid(int facecells, double radius_m, const std::vector<LevelKey>& keys,
                                               bool compact = false);
// A box of the cell-centred levelcolumns x levelcolumns/2 lat-lon grid, in
// row-major box order (local id y columns + x), nothing wrapping inside the box:
// the tables of landscape::finepatchmetric (latitude, row spacing, area, the
// diagonal with the cell's own row beyond the box edge), longitudes as
// makelatlongrid's. Throws std::invalid_argument for an empty or oversized box.
std::shared_ptr<const CellGrid> latlonwindowgrid(int levelcolumns, double radius_m, const LatLonWindowBox& box);
// A planar patch of columns x rows square cells of side cell_m for hand-built benches
// (the analytical relief's tests): a lat-lon window's tables (local id y columns + x,
// nothing wrapping, -1 beyond the edge, globalkey the local id) with every row's
// spacing cell_m, area cell_m^2 and latitude 0 (diagonals hypot(cell_m, cell_m)).
// Throws std::invalid_argument for an empty patch or a cell side that is not positive.
std::shared_ptr<const CellGrid> makeplanarwindow(int columns, int rows, double cell_m);
}
