#pragma once


#include "boundary.h"
#include "prism/field/ifield.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
class IGradient
    : public prism::boundary::BHManagerProvider<boundary::IGradSchemeBoundaryHandler> {
  public:
    IGradient();
    IGradient(const IGradient&) = default;
    IGradient(IGradient&&) = default;
    auto operator=(const IGradient&) -> IGradient& = default;
    auto operator=(IGradient&&) -> IGradient& = default;
    virtual ~IGradient() = default;

    virtual auto gradAtCell(const mesh::Cell& c, const field::IScalar& field) -> Vector3d = 0;

    /// TODO: gradAtCellStored should not require the field as an argument. Remove it.
    virtual auto gradAtCellStored(const mesh::Cell& c,
                                  const field::IScalar& field) -> Vector3d = 0;
    virtual auto gradAtFace(const mesh::Face& f, const field::IScalar& field) -> Vector3d;

  protected:
    virtual auto gradAtBoundaryFace(const mesh::Face& f, const field::IScalar& field) -> Vector3d;
};

} // namespace prism::gradient
