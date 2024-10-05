#pragma once

#include <fmt/format.h>

#include "prism/field/scalar.h"
#include "prism/field/vector.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"
#include "prism/types.h"

namespace prism::ops {

// Calculates the divergence of a vector field U: ∇.U
template <typename Vector>
auto div(const Vector& U) -> field::Scalar;

// Calculates the laplacian of a scalar field ϕ: ∇.∇ϕ
auto laplacian(const field::Scalar& phi, bool return_face_data = true) -> field::Scalar;

// Calculates the magnitude of a vector field U
template <typename Vector>
auto mag(const Vector& U, bool return_face_data = true) -> field::Scalar;

// Calculates the curl of a vector field U
template <typename Vector>
auto curl(const Vector& U, bool return_face_data = true) -> field::Vector;

template <field::IScalarBased Field>
auto grad(const Field& field, Coord coord) -> field::Scalar;

template <field::IScalarBased Field>
auto grad(const Field& field) -> field::Vector;

// face mass flow rate
auto inline faceFlowRate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}

namespace detail {
auto inline coordToIndex(Coord coord) -> std::uint8_t {
    switch (coord) {
        case Coord::X: return 0;
        case Coord::Y: return 1;
        case Coord::Z: return 2;
    }
}

template <typename Vector>
auto fluxSumAtCell(const mesh::PMesh& mesh, const mesh::Cell& cell, const Vector& U) -> double;
} // namespace detail

template <typename Vector>
auto div(const Vector& U) -> field::Scalar {
    // applying green-gauss theorem to the divergence of a vector field
    // ∫(∇.U) dV = ∫U.dS
    // (∇.U) V = Σ U.S
    // ∇.U = (Σ U.S)/V = ops::div(U)
    std::string name = fmt::format("div({})", U.name());
    const mesh::PMesh& mesh = U.mesh();

    VectorXd cell_data;
    cell_data.resize(mesh.cellCount());

    for (const auto& cell : mesh.cells()) {
        cell_data[cell.id()] = detail::fluxSumAtCell(mesh, cell, U) / cell.volume();
    }
    return {name, mesh, cell_data};
}

template <field::IScalarBased Field>
auto grad(const Field& field, Coord coord) -> field::Scalar {
    auto grad_field_name = fmt::format("grad({})_{}", field.name(), field::coordToStr(coord));
    const auto& mesh = field.mesh();

    const auto n_cells = mesh.cellCount();
    const auto n_faces = mesh.faceCount();

    VectorXd grad_values = VectorXd::Zero(n_cells);
    VectorXd grad_face_values = VectorXd::Zero(n_faces);
    auto i = detail::coordToIndex(coord);

    for (std::size_t cell_i = 0; cell_i < n_cells; ++cell_i) {
        grad_values[cell_i] = field.gradAtCell(mesh.cell(cell_i))[i];
    }

    for (std::size_t jface_i = 0; jface_i < n_faces; ++jface_i) {
        grad_face_values[jface_i] = field.gradAtFace(mesh.face(jface_i))[i];
    }

    return field::Scalar(grad_field_name, field.mesh(), grad_values, grad_face_values);
}

template <field::IScalarBased Field>
auto grad(const Field& field) -> field::Vector {
    std::array<field::Scalar, 3> fields {
        grad(field, Coord::X), grad(field, Coord::Y), grad(field, Coord::Z)};

    return field::Vector(fmt::format("grad({})", field.name()), field.mesh(), fields);
}

namespace detail {
template <typename Vector>
auto fluxSumAtCell(const mesh::PMesh& mesh, const mesh::Cell& cell, const Vector& U) -> double {
    double sum = 0.0;

    for (auto face_id : cell.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);

        // Skip empty faces
        if (face.isBoundary()) {
            const auto& boundary_patch = mesh.boundaryPatch(face);
            if (boundary_patch.isEmpty()) {
                continue;
            }
        }
        const Vector3d Uf = U.valueAtFace(face);
        const Vector3d& Sf = mesh::outwardAreaVector(face, cell);
        sum += Uf.dot(Sf);
    }

    return sum;
}

} // namespace detail

} // namespace prism::ops
