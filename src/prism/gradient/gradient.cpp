#include "gradient.h"

#include "prism/exceptions.h"
#include "prism/mesh/utilities.h"

namespace prism::gradient {
IGradient::IGradient(const field::Scalar& field) : _field(field) { // NOLINT
    _bh_manager.add_handler<boundary::Fixed>();
    _bh_manager.add_handler<boundary::FixedGradient>();
    _bh_manager.add_handler<boundary::Empty>();
    _bh_manager.add_handler<boundary::FixedGradient>();
    _bh_manager.add_handler<boundary::Outlet>();
    _bh_manager.add_handler<boundary::Symmetry>();
    _bh_manager.add_handler<boundary::VelocityInlet>();
}

auto IGradient::gradient_at_face(const mesh::Face& face) -> Vector3d {
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
        return (gc * owner_grad) + ((1. - gc) * neighbor_grad);
    }

    // boundary face
    return gradient_at_boundary_face(face);
}

auto IGradient::gradient_at_boundary_face(const mesh::Face& face) -> Vector3d {
    const auto& boundary_patch = _field.mesh().boundary_patch(face);
    const auto& boundary_condition = boundary_patch.get_bc(_field.name());

    auto handler = _bh_manager.get_handler(boundary_condition.kind_string());

    if (handler == nullptr) {
        throw prism::error::NonImplementedBoundaryCondition(
            "IGradient::gradient_at_boundary_face",
            boundary_patch.name(),
            boundary_condition.kind_string());
    }

    return handler->get(_field, face);
}

auto IGradient::gradient_field() -> field::Vector {
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

    auto components_fields = std::array<field::Scalar, 3> {
        field::Scalar(
            grad_field_name + "_x", mesh, std::move(grad_x), std::move(grad_x_face_data)),
        field::Scalar(
            grad_field_name + "_y", mesh, std::move(grad_y), std::move(grad_y_face_data)),
        field::Scalar(
            grad_field_name + "_z", mesh, std::move(grad_z), std::move(grad_z_face_data)),
    };

    return {grad_field_name, mesh, components_fields};
}
} // namespace prism::gradient