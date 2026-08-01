#pragma once


#include "boundary.h"
#include "prism/field/ifield.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
class IGradient : public prism::boundary::BHManagerProvider<boundary::IGradSchemeBoundaryHandler> {
  public:
    IGradient();
    IGradient(const IGradient&) = default;
    IGradient(IGradient&&) = default;
    auto operator=(const IGradient&) -> IGradient& = default;
    auto operator=(IGradient&&) -> IGradient& = default;
    virtual ~IGradient() = default;

    virtual auto gradAtFace(const mesh::Face& f, const field::IScalar& field) -> Vector3d;

    /** @brief Returns the gradient of the given field at a cell.
     *
     * Gradients are cached per field: the first call after the field's cell values changed
     * recomputes the gradient over the whole mesh (one pass) and stamps the field's event
     * number; subsequent calls return the stored value in O(1). Callers never need to
     * explicitly refresh the gradient cache. */
    virtual auto gradAtCell(const mesh::Cell& c, const field::IScalar& field) -> Vector3d = 0;

    /** @brief Forces a whole-mesh recomputation of the cell gradients of the given field.
     * After this call, gradAtCell() returns the freshly computed gradients in O(1) until the
     * field's cell values change again. */
    virtual void computeAllCellGradients(const field::IScalar& field) = 0;

  protected:
    /** @brief Computes and stores the gradient at a single cell (no staleness check).
     * Used by computeAllCellGradients() and the self-healing gradAtCell() implementation. */
    virtual auto computeGradAtCell(const mesh::Cell& c, const field::IScalar& field)
        -> Vector3d = 0;

    virtual auto gradAtBoundaryFace(const mesh::Face& f, const field::IScalar& field) -> Vector3d;
};

} // namespace prism::gradient
