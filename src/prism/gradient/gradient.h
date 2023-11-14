#pragma once

#include "prism/constants.h"
#include "prism/field.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"
#include "prism/types.h"


namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class AbstractGradient {
  public:
    AbstractGradient() = delete;
    AbstractGradient(const ScalarField& field) : _field(field) {} // NOLINT
    AbstractGradient(const AbstractGradient&) = default;
    AbstractGradient(AbstractGradient&&) = default;
    auto operator=(const AbstractGradient&) -> AbstractGradient& = delete;
    auto operator=(AbstractGradient&&) -> AbstractGradient& = delete;

    virtual auto gradient_at_cell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradient_at_face(const mesh::Face& f) -> Vector3d;
    virtual auto gradient_field() -> VectorField;

    static auto gradient_at_boundary_face(const mesh::Face& f, const ScalarField& field)
        -> Vector3d;

    virtual ~AbstractGradient() = default;

  private:
    const ScalarField _field;
};

class GreenGauss : public AbstractGradient {
  public:
    GreenGauss(const ScalarField& field)
        : _field(field),
          _cell_gradients(MatrixX3d::Zero(field.mesh().n_cells(), 3)),
          AbstractGradient(field) {}

    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;

  private:
    const ScalarField& _field;
    MatrixX3d _cell_gradients;
};

class LeastSquares : public AbstractGradient {
  public:
    LeastSquares(const ScalarField& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;

  private:
    void set_lsq_matrices();
    auto boundary_face_phi(const mesh::Face& face) -> double;

    const ScalarField& _field;
    std::vector<MatrixX3d> _pinv_matrices; // pseudo-inverse matrices
};

auto inline AbstractGradient::gradient_at_boundary_face(const mesh::Face& face,
                                                        const ScalarField& field) -> Vector3d {
    const auto& boundary_patch = field.mesh().boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(field.name());
    auto bc_type = boundary_condition.bc_type();

    switch (bc_type) {
        case mesh::BoundaryConditionType::Empty:
        case mesh::BoundaryConditionType::Outlet:
        case mesh::BoundaryConditionType::Symmetry: {
            return Vector3d {0., 0., 0.};
        }

        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            const auto& phi_name = field.name();
            auto phi = boundary_patch.get_scalar_bc(phi_name);
            return phi * face.area_vector();
        }

        case mesh::BoundaryConditionType::FixedGradient: {
            const auto& phi_name = field.name();
            auto flux = boundary_patch.get_scalar_bc(phi_name + "-flux");
            return flux * face.area_vector();
        }

        default:
            throw std::runtime_error(
                fmt::format("GradientSchemeBase::gradient_at_boundary_face(): "
                            "Non-implemented boundary type for boundary patch: '{}'",
                            boundary_patch.name()));
    }
}

auto inline AbstractGradient::gradient_at_face(const mesh::Face& face) -> Vector3d {
    // interpolate gradient at surrounding cells to the face center
    // interior face
    if (face.has_neighbor()) {
        const auto& mesh = _field.mesh();
        const auto& owner_cell = mesh.cell(face.owner());
        auto owner_grad = gradient_at_cell(owner_cell);

        const auto& neighbor_cell = mesh.cell(face.neighbor().value());
        auto neighbor_grad = gradient_at_cell(neighbor_cell);

        auto gc = mesh::geo_weight(owner_cell, neighbor_cell, face);

        // Equation 9.33 without the correction part, a simple linear interpolation.
        Vector3d grad = gc * owner_grad + (1. - gc) * neighbor_grad;

        // correct the interpolation
        auto phi_C = _field[owner_cell.id()];
        auto phi_F = _field[neighbor_cell.id()];
        auto d_CF = neighbor_cell.center() - owner_cell.center();
        auto d_CF_norm = d_CF.norm();
        Vector3d e_CF = d_CF / d_CF_norm;

        auto correction = (phi_F - phi_C) / (d_CF_norm + EPSILON);
        correction -= grad.dot(e_CF);

        return grad + (e_CF * correction);
    }

    // boundary face
    return gradient_at_boundary_face(face, _field);
}

auto inline AbstractGradient::gradient_field() -> VectorField {
    auto grad_field_name = _field.name() + "_grad";
    VectorField grad_field {grad_field_name, _field.mesh()};

    for (const auto& cell : _field.mesh().cells()) {
        grad_field.data().row(cell.id()) = gradient_at_cell(cell);
    }

    return grad_field;
}

} // namespace prism::gradient