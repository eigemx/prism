#pragma once

#include <fmt/format.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "prism/boundary.h"
#include "prism/constants.h"
#include "prism/field/scalar.h"
#include "prism/field/vector.h"
#include "prism/gradient/boundary.h"
#include "prism/log.h"
#include "prism/mesh/face.h"
#include "prism/mesh/utilities.h"
#include "prism/types.h"

namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
template <field::IScalarBased Field>
class IGradient
    : public prism::boundary::BHManagersProvider<boundary::IGradSchemeBoundaryHandler> {
  public:
    IGradient() = delete;
    IGradient(const IGradient&) = default;
    IGradient(IGradient&&) = default;
    auto operator=(const IGradient&) -> IGradient& = default;
    auto operator=(IGradient&&) -> IGradient& = default;
    virtual ~IGradient() = default;

    IGradient(Field field);

    virtual auto gradAtCell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradAtFace(const mesh::Face& f) -> Vector3d;
    virtual auto gradField() -> field::Vector;
    auto field() -> Field { return _field; }

  protected:
    virtual auto gradAtBoundaryFace(const mesh::Face& f) -> Vector3d;

  private:
    Field _field;
};

template <field::IScalarBased Field>
class GreenGauss : public IGradient<Field> {
  public:
    GreenGauss(Field field);
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;

  private:
    auto correctSkewness(const mesh::Face& face,
                         const mesh::Cell& cell,
                         const mesh::Cell& nei) const -> double;
    auto gradAtCell_(const mesh::Cell& cell, bool correct_skewness = true) -> Vector3d;
    auto boundaryFaceIntegral(const mesh::Face& f) -> Vector3d;

    std::vector<Vector3d> _cell_gradients;
};

template <field::IScalarBased Field>
class LeastSquares : public IGradient<Field> {
  public:
    LeastSquares(Field field);
    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;

  private:
    void setPseudoInvMatrices();

    std::vector<MatrixX3d> _pinv_matrices; // pseudo-inverse matrices
};


template <field::IScalarBased Field, typename G>
class GradientProvider {
  public:
    GradientProvider(Field field) : _gradient_scheme(field) {}
    using GradSchemeType = G;
    auto gradScheme() -> GradSchemeType& { return _gradient_scheme; }

  private:
    GradSchemeType _gradient_scheme;
};

template <field::IScalarBased Field>
IGradient<Field>::IGradient(Field field) : _field(field) { // NOLINT
    log::debug("prism::gradient::IGradient() adding default boundary handlers for IGradient");
    this->boundaryHandlersManager().template addHandler<boundary::Fixed>();
    this->boundaryHandlersManager().template addHandler<boundary::FixedGradient>();
    this->boundaryHandlersManager().template addHandler<boundary::Empty>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient>();
    this->boundaryHandlersManager().template addHandler<boundary::VelocityInlet>();
}

template <field::IScalarBased Field>
auto IGradient<Field>::gradAtFace(const mesh::Face& face) -> Vector3d {
    // interpolate gradient at surrounding cells to the face center
    if (face.isInterior()) {
        // interior face
        const auto& mesh = _field.mesh();
        const auto& owner_cell = mesh.cell(face.owner());
        auto owner_grad = gradAtCell(owner_cell);

        const auto& neighbor_cell = mesh.cell(face.neighbor().value());
        auto neighbor_grad = gradAtCell(neighbor_cell);

        auto gc = mesh::geometricWeight(owner_cell, neighbor_cell, face);

        // Equation 9.33 without the correction part, a simple linear interpolation.
        return (gc * owner_grad) + ((1. - gc) * neighbor_grad);
    }

    // boundary face
    return gradAtBoundaryFace(face);
}

template <field::IScalarBased Field>
auto IGradient<Field>::gradAtBoundaryFace(const mesh::Face& face) -> Vector3d {
    const auto& boundary_patch = _field.mesh().boundaryPatch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(_field.name());

    auto handler = this->boundaryHandlersManager().getHandler(boundary_condition.kindString());

    if (handler == nullptr) {
        throw prism::error::NonImplementedBoundaryCondition(
            "prism::gradient::IGradient::gradAtBoundaryFace",
            boundary_patch.name(),
            boundary_condition.kindString());
    }

    return handler->get(_field, face);
}

template <field::IScalarBased Field>
auto IGradient<Field>::gradField() -> field::Vector {
    // TODO: This function is VERY expensive
    auto grad_field_name = fmt::format("grad({})", _field.name());
    const auto& mesh = _field.mesh();

    auto n_cells = mesh.nCells();
    auto n_faces = mesh.nFaces();

    VectorXd grad_x = VectorXd::Zero(n_cells);
    VectorXd grad_x_face_data = VectorXd::Zero(n_faces);

    VectorXd grad_y = VectorXd::Zero(n_cells);
    VectorXd grad_y_face_data = VectorXd::Zero(n_faces);

    VectorXd grad_z = VectorXd::Zero(n_cells);
    VectorXd grad_z_face_data = VectorXd::Zero(n_faces);

    for (std::size_t i = 0; i < n_cells; ++i) {
        const auto& cell_grad = gradAtCell(mesh.cell(i));
        grad_x[i] = cell_grad[0];
        grad_y[i] = cell_grad[1];
        grad_z[i] = cell_grad[2];
    }

    for (std::size_t j = 0; j < n_faces; ++j) {
        const auto& face_grad = gradAtFace(mesh.face(j));
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

template <field::IScalarBased Field>
auto GreenGauss<Field>::correctSkewness(const mesh::Face& face,
                                        const mesh::Cell& cell,
                                        const mesh::Cell& nei) const -> double {
    auto grad_sum = _cell_gradients[cell.id()] + _cell_gradients[nei.id()];
    auto vec = face.center() - (0.5 * (cell.center() + nei.center()));

    return 0.5 * grad_sum.dot(vec);
}

template <field::IScalarBased Field>
GreenGauss<Field>::GreenGauss(Field field) : IGradient<Field>(field) { // NOLINT
    // We need to perform a first run for calculating gradient at cells,
    // to make the cell gradient vector _cell_gradients available if the user desires to call
    // gradient_at_face() which requires a first run of gradient calculations, to perform
    // correction for faces with skewness.
    const std::size_t n_cells = this->field().mesh().nCells();


    _cell_gradients.reserve(n_cells);
    for (const auto& cell : this->field().mesh().cells()) {
        // caclulate the gradient without skewness correction
        _cell_gradients.emplace_back(gradAtCell_(cell, false));
    }
}

template <field::IScalarBased Field>
auto GreenGauss<Field>::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    return gradAtCell_(cell, true);
}

template <field::IScalarBased Field>
auto GreenGauss<Field>::gradAtCell_(const mesh::Cell& cell, bool correct_skewness) -> Vector3d {
    Vector3d grad {0., 0., 0.};
    const auto& mesh = this->field().mesh();

    for (auto face_id : cell.facesIds()) {
        const auto& face = mesh.face(face_id);

        // This is a boundary face
        if (face.is_boundary()) {
            grad += boundaryFaceIntegral(face);
            continue;
        }

        // This is an internal face
        // Area normal vector, poitning out of the cell
        auto Sf = mesh::outwardAreaVector(face, cell);
        const auto& nei = this->field().mesh().otherSharingCell(cell, face);
        auto face_phi = 0.5 * (this->field()[cell.id()] + this->field()[nei.id()]);

        if (correct_skewness) {
            face_phi += correctSkewness(face, cell, nei);
        }
        grad += Sf * face_phi;
    }
    grad /= cell.volume();

    // store the gradient to use it in next iterations, for skewness correction
    _cell_gradients[cell.id()] = grad;

    return grad;
}

template <field::IScalarBased Field>
auto GreenGauss<Field>::boundaryFaceIntegral(const mesh::Face& face) -> Vector3d {
    // returns the Green-Gauss face integral at boundary face `face`
    const auto& boundary_patch = this->field().mesh().boundary_patch(face);
    const auto& boundary_condition = boundary_patch.getBoundaryCondition(this->field().name());

    if (boundary_patch.isEmpty()) {
        return {0.0, 0.0, 0.0};
    }

    auto phi = this->field().valueAtFace(face);
    return phi * face.area_vector();
}

template <field::IScalarBased Field>
LeastSquares<Field>::LeastSquares(Field field) : IGradient<Field>(field) {
    setPseudoInvMatrices();
}

template <field::IScalarBased Field>
void LeastSquares<Field>::setPseudoInvMatrices() {
    // This function is based on section 9.3 'Least-Square Gradient'
    const auto& mesh = this->field().mesh();

    // resize the pseudo-inverse matrices vector
    _pinv_matrices.resize(mesh.nCells());

    for (const auto& cell : mesh.cells()) {
        // A 3x3 matrix of the left hand side of equation (9.27)
        // for the k-th cell, we calculate this distance matrix D
        // and push the pseudo-inverse [(D * D^T)^{-1} * D^T] to _pinv_matrix vector
        Matrix3d d_matrix = Matrix3d::Zero();

        for (auto face_id : cell.facesIds()) {
            const auto& face = mesh.face(face_id);

            // This will hold the distance vector from neighbor cell center to k-th cell
            // center, or in case we have a boundary face, r_CF will be the distance vector
            // from boundary face center to the k-th cell center.
            // check equation (9.22)
            Vector3d r_CF = {.0, .0, .0};

            // difference of field values between the k-th cell and its i-th neigbor
            // or its boundary face field value
            double delta_phi = 0.0;

            Matrix3d d_matrix_k = Matrix3d::Zero();

            if (face.isInterior()) {
                // interior face
                const auto neighbor = mesh.otherSharingCell(cell, face);
                r_CF = neighbor.center() - cell.center();
                delta_phi = this->field().valueAtCell(neighbor) - this->field().valueAtCell(cell);
            } else {
                // boundary face
                r_CF = face.center() - cell.center();
                delta_phi = this->field().valueAtFace(face) - this->field().valueAtCell(cell);
            }

            // weight factor defined in equation (9.28)
            const double wk = 1 / (r_CF.norm() + EPSILON);
            const double dx = r_CF.x(); // equation (9.24)
            const double dy = r_CF.y(); // equation (9.24)
            const double dz = r_CF.z(); // equation (9.24)

            // clang-format off
            // left hand side of equation (9.27) for the k-th cell, before summing and before
            // multiplying the weight factor wk
            d_matrix_k << (dx * dx), (dx * dy), (dx * dz), 
                          (dy * dx), (dy * dy), (dy * dz),
                          (dz * dx), (dz * dy), (dz * dz);
            // clang-format on

            d_matrix += d_matrix_k * wk;
        }

        const auto& d_matrix_t = d_matrix.transpose();
        _pinv_matrices[cell.id()] = (d_matrix_t * d_matrix).inverse() * d_matrix_t;
    }
}

template <field::IScalarBased Field>
auto LeastSquares<Field>::gradAtCell(const mesh::Cell& cell) -> Vector3d {
    const auto& mesh = this->field().mesh();

    // right hand side of equation (9.27)
    Vector3d b {0.0, 0.0, 0.0};

    for (auto face_id : cell.facesIds()) {
        const auto& face = mesh.face(face_id);

        double delta_phi = 0.0;
        auto phi_cell = this->field().valueAtCell(cell);
        Vector3d r_CF = {.0, .0, .0};

        if (face.isInterior()) {
            // interior face
            const auto neighbor = mesh.otherSharingCell(cell, face);
            r_CF = neighbor.center() - cell.center();
            auto nei_phi = this->field().valueAtCell(neighbor);
            delta_phi = nei_phi - phi_cell;

        } else {
            // boundary face
            auto bface_phi = this->field().valueAtFace(face);
            r_CF = face.center() - cell.center();
            delta_phi = bface_phi - phi_cell;
        }
        const double wk = 1 / (r_CF.norm() + EPSILON);
        const double dx = r_CF.x();
        const double dy = r_CF.y();
        const double dz = r_CF.z();

        // update the right hand side of equation (9.27)
        b += Vector3d {dx * delta_phi, dy * delta_phi, dz * delta_phi} * wk;
    }
    return _pinv_matrices[cell.id()] * b;
}
} // namespace prism::gradient