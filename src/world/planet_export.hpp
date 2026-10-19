#pragma once

// The lat-lon export and import of a world (docs/architecture/sphere-mesh.md
// section 6). Every field of a spatial layout on the world grid moves by its
// catalogue remap rule (world_fields.hpp, RemapRule): mean, total or nearest;
// direction fields are never remapped but recomputed on the export by a
// callback (the lat-lon flow routing, supplied by hydrology) or exported as 0.
// Scalars, zonal profiles, tables and the climate solver's own grids are copied.

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
// callback, or 0 without one. Settlements and deposits move to the export's cell of
// their cell's centre (x, y its column and row); polities and routes have no positions. The
// planet is returned on the heap: it holds the crater list inline and has no move
// constructor.
std::unique_ptr<planet> exportlatlon(const planet& world, int columns, const LatlonDirectionRecompute& recompute = {});

// The reverse for the named fields (climate's outputs): each is remapped from
// the lat-lon planet onto the world's grid by its rule (copied when the grids are
// equal, so the import is the identity on a lat-lon world of the same columns). A
// direction field onto another grid is an error; scalars and zonal profiles are copied.
void importfields(planet& world, const planet& latlon, const std::vector<worldmodel::FieldId>& fields);
