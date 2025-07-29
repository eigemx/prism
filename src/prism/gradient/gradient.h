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
    IGradient() = delete;
    IGradient(const IGradient&) = default;
    IGradient(IGradient&&) = default;
    auto operator=(const IGradient&) -> IGradient& = default;
    auto operator=(IGradient&&) -> IGradient& = default;
    virtual ~IGradient() = default;

    IGradient(field::IScalar* field);

    virtual auto gradAtCell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradAtCellStored(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradAtFace(const mesh::Face& f) -> Vector3d;

  protected:
    auto field() -> field::IScalar* { return _field; }
    virtual auto gradAtBoundaryFace(const mesh::Face& f) -> Vector3d;

  private:
    field::IScalar* _field;
};

} // namespace prism::gradient
