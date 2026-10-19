#pragma once

// The lat-lon export and import of a world (docs/architecture/sphere-mesh.md
// section 6). Every field of a spatial layout on the world grid moves by its
// catalogue remap rule (world_fields.hpp, RemapRule): mean, total or nearest;
// direction fields are never remapped but recomputed on the export by a
// callback (the lat-lon flow routing, supplied by hydrology) or exported as 0;
// land-only and sea-only codes (landawarefield) by class against the destination's
// own sea mask. Scalars, zonal profiles, tables and the climate solver's own grids are copied.

#include <cstdint>
#include <functional>
#include <memory>
#include <vector>

#include "core/mesh/cell_grid.hpp"
#include "core/mesh/mesh_remap.hpp"
#include "core/world/field_encoding.hpp"

class planet;

// Recomputes the direction fields (rule direction) of a lat-lon export from its
// exported fields (the filled elevation).
using LatlonDirectionRecompute = std::function<void(planet& latlon)>;

// One field's data moved between two full grids by its catalogue rule. The data
// is in the planet's encoding (field_encoding.hpp): [layer,] rows, columns on
// lat-lon, [layer,] cells on the cube. The weights are built on first use; the
// grids must outlive the remap.
class FieldRemap
{
public:
    FieldRemap(const mesh::CellGrid& source, const mesh::CellGrid& destination, mesh::RemapOptions options = {});
    ~FieldRemap();

    const mesh::CellGrid& source() const { return source_; }
    const mesh::CellGrid& destination() const { return destination_; }
    // Equal full grids (the same kind and size): every field is copied as it is.
    bool identity() const { return identity_; }

    // True when the data's trailing axes are the source grid's cells: a field on
    // the world grid (not a scalar, a zonal profile or a solver's own grid).
    bool ongrid(const worldmodel::FieldData& data) const;
    // The data on the destination grid by the field's rule; a copy between equal
    // grids. A direction field between other grids is an error (std::logic_error).
    worldmodel::FieldData remap(const worldmodel::FieldData& data) const;
    // A land-aware code (landawarefield) by class: each destination cell takes the
    // code with the largest overlap among the source cells of its own land/sea class
    // (1 sea: surfacebelowsea, per cell of each grid), ties to the smaller code; the
    // closest source cell of its class within two neighbour rings of the cell at its
    // centre where none overlaps, else that cell. A copy between equal grids.
    worldmodel::FieldData remapbyclass(const worldmodel::FieldData& data, const std::vector<std::uint8_t>& sourcesea,
        const std::vector<std::uint8_t>& destinationsea) const;
    // Zeros of the field's type in the destination shape.
    worldmodel::FieldData zeros(const worldmodel::FieldData& data) const;

private:
    std::size_t layers(const worldmodel::FieldData& data) const;
    std::vector<std::size_t> destinationshape(const worldmodel::FieldData& data) const;
    const mesh::RemapWeights& weights() const;

    const mesh::CellGrid& source_;
    const mesh::CellGrid& destination_;
    mesh::RemapOptions options_;
    bool identity_ = false;
    mutable std::unique_ptr<mesh::RemapWeights> weights_;
};

// A lat-lon planet of columns x columns/2 with every field of the world (a copy
// on a lat-lon world of the same columns). Direction fields are recomputed by the
// callback, or 0 without one; land-aware codes are moved after the terrain, against
// the export's own sea mask. Settlements and deposits move to the export's cell of
// their cell's centre (x, y its column and row); polities and routes have no positions. The
// planet is returned on the heap: it holds the crater list inline and has no move
// constructor.
std::unique_ptr<planet> exportlatlon(const planet& world, int columns, const LatlonDirectionRecompute& recompute = {});

// Codes and masks that hold only on land or only at sea, moved by class (FieldRemap::
// remapbyclass) on export and import: on land climate_class, biome, special_features,
// island, mountain_island, basin_class, prograded, owner_settlement, owner_polity; at
// sea sea_ice, continental_shelf. Continuous fields and other codes move by their rule.
bool landawarefield(worldmodel::FieldId id);
// The sea mask of a world by cell id (1 where surfacebelowsea).
std::vector<std::uint8_t> seamask(const planet& world);

// The reverse for the named fields (climate's outputs): each is remapped from
// the lat-lon planet onto the world's grid by its rule (copied when the grids are
// equal, so the import is the identity on a lat-lon world of the same columns), a
// land-aware code by class against the world's own sea mask. A direction field onto
// another grid is an error; scalars and zonal profiles are copied.
void importfields(planet& world, const planet& latlon, const std::vector<worldmodel::FieldId>& fields);
