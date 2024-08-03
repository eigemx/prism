#pragma once

#include <fmt/core.h>

#include <optional>
#include <vector>

#include "prism/boundary.h"
#include "prism/field/field.h"
#include "prism/gradient/boundary.h"
#include "prism/mesh/face.h"
#include "prism/types.h"


namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class IGradient {
  public:
    IGradient() = delete;
    IGradient(const field::Scalar& field);
    IGradient(const IGradient&) = default;
    IGradient(IGradient&&) = default;
    auto operator=(const IGradient&) -> IGradient& = delete;
    auto operator=(IGradient&&) -> IGradient& = delete;
    virtual ~IGradient() = default;

    virtual auto gradient_at_cell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradient_at_face(const mesh::Face& f) -> Vector3d;
    virtual auto gradient_field() -> field::Vector;

    using BoundaryHandlersManager =
        prism::boundary::BoundaryHandlersManager<IGradient, boundary::GradSchemeBoundaryHandler>;
    auto bh_manager() -> BoundaryHandlersManager& { return _bh_manager; }

  protected:
    virtual auto gradient_at_boundary_face(const mesh::Face& f) -> Vector3d;

  private:
    field::Scalar _field;
    BoundaryHandlersManager _bh_manager;
};

class GreenGauss : public IGradient {
  public:
    GreenGauss(const field::Scalar& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;

  private:
    auto skewness_correction(const mesh::Face& face,
                             const mesh::Cell& cell,
                             const mesh::Cell& nei) const -> double;
    auto _gradient_at_cell(const mesh::Cell& cell, bool correct_skewness = true) -> Vector3d;
    auto green_gauss_face_integral(const mesh::Face& f) -> Vector3d;

    field::Scalar _field;
    std::vector<Vector3d> _cell_gradients;
};

class LeastSquares : public IGradient {
  public:
    LeastSquares(const field::Scalar& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;

  private:
    void set_pseudo_inv_matrices();
    auto boundary_face_phi(const mesh::Face& face) -> std::optional<double>;

    field::Scalar _field;
    std::vector<MatrixX3d> _pinv_matrices; // pseudo-inverse matrices
};

} // namespace prism::gradient