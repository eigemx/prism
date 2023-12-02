#pragma once

#include <fmt/core.h>

#include <optional>
#include <stdexcept>
#include <vector>

#include "prism/constants.h"
#include "prism/exceptions.h"
#include "prism/field.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/face.h"
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

    virtual ~AbstractGradient() = default;

  protected:
    virtual auto gradient_at_boundary_face(const mesh::Face& f) -> Vector3d;

  private:
    const ScalarField _field;
};

class GreenGauss : public AbstractGradient {
  public:
    GreenGauss(const ScalarField& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;

  private:
    auto skewness_correction(const mesh::Face& face,
                             const mesh::Cell& cell,
                             const mesh::Cell& nei) const -> double;
    auto _gradient_at_cell(const mesh::Cell& cell, bool correct_skewness = true) -> Vector3d;
    auto green_gauss_face_integral(const mesh::Face& f) -> Vector3d;

    const ScalarField _field;
    std::vector<Vector3d> _cell_gradients;
};

class LeastSquares : public AbstractGradient {
  public:
    LeastSquares(const ScalarField& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;

  private:
    void set_pseudo_inv_matrices();
    auto boundary_face_phi(const mesh::Face& face) -> std::optional<double>;

    const ScalarField _field;
    std::vector<MatrixX3d> _pinv_matrices; // pseudo-inverse matrices
};

auto inline AbstractGradient::gradient_at_face(const mesh::Face& face) -> Vector3d {
    // interpolate gradient at surrounding cells to the face center
    if (face.is_interior()) {
        // interior face
        const auto& mesh = _field.mesh();
        const auto& owner_cell = mesh.cell(face.owner());
        auto owner_grad = gradient_at_cell(owner_cell);

        const auto& neighbor_cell = mesh.cell(face.neighbor().value());
        auto neighbor_grad = gradient_at_cell(neighbor_cell);

        auto gc = mesh::geo_weight(owner_cell, neighbor_cell, face);

        // Equation 9.33 without the correction part, a simple linear interpolation.
        Vector3d grad = (gc * owner_grad) + ((1. - gc) * neighbor_grad);

        // correct the interpolation
        auto phi_C = _field.value_at_cell(owner_cell.id());
        auto phi_F = _field.value_at_cell(neighbor_cell.id());

        auto d_CF = neighbor_cell.center() - owner_cell.center();
        auto d_CF_norm = d_CF.norm();
        Vector3d e_CF = d_CF / d_CF_norm;

        auto correction = (phi_F - phi_C) / (d_CF_norm + EPSILON);
        correction -= grad.dot(e_CF);

        return grad + (e_CF * correction);
    }

    // boundary face
    return gradient_at_boundary_face(face);
}

auto inline AbstractGradient::gradient_at_boundary_face(const mesh::Face& face) -> Vector3d {
    const auto& boundary_patch = _field.mesh().boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_field.name());
    auto bc_type = boundary_condition.kind();

    switch (bc_type) {
        case mesh::BoundaryConditionKind::Empty:
        case mesh::BoundaryConditionKind::Symmetry:
        case mesh::BoundaryConditionKind::Outlet: {
            return {0.0, 0.0, 0.0};
        }

        case mesh::BoundaryConditionKind::Inlet:
        case mesh::BoundaryConditionKind::Fixed: {
            const auto& owner = _field.mesh().cell(face.owner());
            Vector3d e = face.center() - owner.center();
            double d_Cf = e.norm();
            e = e / e.norm();

            double delta_phi = _field.value_at_face(face) - _field.value_at_cell(owner);
            return (delta_phi / d_Cf) * e;
        }

        case mesh::BoundaryConditionKind::FixedGradient: {
            return boundary_patch.get_vector_bc(_field.name());
        }

        default: {
            throw error::NonImplementedBoundaryCondition(
                "AbstractGradient::gradient_at_boundary_face()",
                boundary_patch.name(),
                boundary_condition.kind_string());
        }
    }

    // We should never reach this
    return {0.0, 0.0, 0.0};
}

auto inline AbstractGradient::gradient_field() -> VectorField {
    // TODO: This function is VERY expensive
    auto grad_field_name = fmt::format("grad({})", _field.name());
    const auto& mesh = _field.mesh();

    auto n_cells = mesh.n_cells();
    auto n_faces = mesh.n_faces();

    VectorXd grad_x = VectorXd::Zero(n_cells);
    VectorXd grad_x_face_data = VectorXd::Zero(n_faces);

    VectorXd grad_y = VectorXd::Zero(n_cells);
    VectorXd grad_y_face_data = VectorXd::Zero(n_faces);

    VectorXd grad_z = VectorXd::Zero(n_cells);
    VectorXd grad_z_face_data = VectorXd::Zero(n_faces);

    for (std::size_t i = 0; i < n_cells; ++i) {
        const auto& cell_grad = gradient_at_cell(mesh.cell(i));
        grad_x[i] = cell_grad[0];
        grad_y[i] = cell_grad[1];
        grad_z[i] = cell_grad[2];
    }

    for (std::size_t j = 0; j < n_faces; ++j) {
        const auto& face_grad = gradient_at_face(mesh.face(j));
        grad_x_face_data[j] = face_grad[0];
        grad_y_face_data[j] = face_grad[1];
        grad_z_face_data[j] = face_grad[2];
    }

    auto components_fields = std::array<ScalarField, 3> {
        ScalarField(grad_field_name + "_x", mesh, std::move(grad_x), std::move(grad_x_face_data)),
        ScalarField(grad_field_name + "_y", mesh, std::move(grad_y), std::move(grad_y_face_data)),
        ScalarField(grad_field_name + "_z", mesh, std::move(grad_z), std::move(grad_z_face_data)),
    };

    return {grad_field_name, mesh, components_fields};
}

} // namespace prism::gradient