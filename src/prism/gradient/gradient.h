#pragma once

#include <fmt/core.h>

#include <vector>

#include "prism/boundary.h"
#include "prism/field/scalar.h"
#include "prism/field/vector.h"
#include "prism/gradient/boundary.h"
#include "prism/mesh/face.h"
#include "prism/types.h"


namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class IGradient {
  public:
    IGradient() = delete;
    IGradient(const IGradient&) = default;
    IGradient(IGradient&&) = default;
    auto operator=(const IGradient&) -> IGradient& = default;
    auto operator=(IGradient&&) -> IGradient& = default;
    virtual ~IGradient() = default;

    IGradient(field::Scalar field);

    virtual auto gradAtCell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradAtFace(const mesh::Face& f) -> Vector3d;
    virtual auto gradField() -> field::Vector;

    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<IGradient, boundary::GradSchemeBoundaryHandler>;
    auto boundaryHandlersManager() -> BoundaryHandlersManager& { return _bh_manager; }

  protected:
    virtual auto gradAtBoundaryFace(const mesh::Face& f) -> Vector3d;

  private:
    field::Scalar _field;
    BoundaryHandlersManager _bh_manager;
};

class GreenGauss : public IGradient {
  public:
    GreenGauss(field::Scalar field);
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;

  private:
    auto correctSkewness(const mesh::Face& face,
                         const mesh::Cell& cell,
                         const mesh::Cell& nei) const -> double;
    auto gradAtCell_(const mesh::Cell& cell, bool correct_skewness = true) -> Vector3d;
    auto boundaryFaceIntegral(const mesh::Face& f) -> Vector3d;

    field::Scalar _field;
    std::vector<Vector3d> _cell_gradients;
};

class LeastSquares : public IGradient {
  public:
    LeastSquares(field::Scalar field);
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;

  private:
    void setPseudoInvMatrices();

    field::Scalar _field;
    std::vector<MatrixX3d> _pinv_matrices; // pseudo-inverse matrices
};

} // namespace prism::gradient