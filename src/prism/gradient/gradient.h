#pragma once

#include <fmt/format.h>

#include <vector>

#include "boundary.h"
#include "prism/field/ifield.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::gradient {
// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class IGradient
    : public prism::boundary::BHManagersProvider<boundary::IGradSchemeBoundaryHandler> {
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

class GreenGauss : public IGradient {
  public:
    GreenGauss(field::IScalar* field);
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) -> Vector3d override;

  private:
    auto correctSkewness(const mesh::Face& face,
                         const mesh::Cell& cell,
                         const mesh::Cell& nei) const -> double;
    auto gradAtCell_(const mesh::Cell& cell, bool correct_skewness = true) -> Vector3d;
    auto boundaryFaceIntegral(const mesh::Face& f) -> Vector3d;

    std::vector<Vector3d> _cell_gradients;
};

class LeastSquares : public IGradient {
  public:
    LeastSquares(field::IScalar* field);
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) -> Vector3d override;

  private:
    void setPseudoInvMatrices();

    std::vector<Vector3d> _cell_gradients;
    std::vector<Matrix3d> _pinv_matrices; // pseudo-inverse matrices
};

} // namespace prism::gradient