#pragma once

#include <fmt/format.h>

#include <cstdint>

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
auto div(const Vector& U, bool return_face_data = true) -> field::Scalar;

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
    std::uint8_t i = 0;
    switch (coord) {
        case Coord::X: {
            i = 0;
            break;
        }

        case Coord::Y: {
            i = 1;
            break;
        }

        case Coord::Z: {
            i = 2;
            break;
        }
    }
    return i;
}
template <typename Vector>
auto divAtCell(const mesh::PMesh& mesh, const mesh::Cell& cell, const Vector& U) -> double;

template <typename Vector>
auto fluxAtFace(const mesh::PMesh& mesh,
                const mesh::Cell& cell,
                const mesh::Face& face,
                const Vector& U) -> double;

template <typename Vector>
auto fluxAtBoundaryFace(const mesh::PMesh& mesh,
                        const mesh::Face& face,
                        const Vector& U) -> double;
} // namespace detail

template <typename Vector>
auto div(const Vector& U, bool return_face_data) -> field::Scalar {
    std::string name = fmt::format("div({})", U.name());
    const mesh::PMesh& mesh = U.mesh();

    VectorXd cell_data;
    cell_data.resize(mesh.nCells());

    for (const auto& cell : mesh.cells()) {
        cell_data[cell.id()] = detail::divAtCell(mesh, cell, U);
    }

    if (return_face_data) {
        VectorXd face_data;
        face_data.resize(mesh.nFaces());

        // We start with calculating the fluxes at boundary faces
        for (const auto& bface : mesh.boundaryFaces()) {
            face_data[bface.id()] = detail::fluxAtBoundaryFace(mesh, bface, U);
        }

        // for interior faces, we take the average of the divergence at the two sharing cells
        for (const auto& iface : mesh.interiorFaces()) {
            const auto& owner = mesh.cell(iface.owner());
            const auto& neighbor = mesh.cell(iface.neighbor().value());
            auto gc = mesh::geometricWeight(owner, neighbor, iface);

            auto div_f = gc * cell_data[owner.id()];
            div_f += (1 - gc) * cell_data[neighbor.id()];

            face_data[iface.id()] = div_f;
        }
        return {name, mesh, cell_data, face_data};
    }

    return {name, mesh, cell_data};
}

template <field::IScalarBased Field>
auto grad(const Field& field, Coord coord) -> field::Scalar {
    auto grad_field_name = fmt::format("grad({})_{}", field.name(), field::coordToStr(coord));
    const auto& mesh = field.mesh();

    const auto n_cells = mesh.nCells();
    const auto n_faces = mesh.nFaces();

    VectorXd grad_values = VectorXd::Zero(n_cells);
    VectorXd grad_face_values = VectorXd::Zero(n_faces);
    auto i = detail::coordToIndex(coord);

    for (std::size_t j = 0; j < n_cells; ++j) {
        grad_values[j] = field.gradAtCell(mesh.cell(j))[i];
    }

    for (std::size_t j = 0; j < n_faces; ++j) {
        grad_face_values[j] = field.gradAtFace(mesh.face(j))[i];
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
auto divAtCell(const mesh::PMesh& mesh, const mesh::Cell& cell, const Vector& U) -> double {
    double sum = 0.0;

    for (auto face_id : cell.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        sum += fluxAtFace(mesh, cell, face, U);
    }

    return sum;
}

template <typename Vector>
auto fluxAtFace(const mesh::PMesh& mesh,
                const mesh::Cell& cell,
                const mesh::Face& face,
                const Vector& U) -> double {
    if (face.isBoundary()) {
        return fluxAtBoundaryFace(mesh, face, U);
    }

    const Vector3d Uf = U.valueAtFace(face);
    auto Sf = mesh::outwardAreaVector(face, cell);
    return Uf.dot(Sf);
}

template <typename Vector>
auto fluxAtBoundaryFace(const mesh::PMesh& mesh,
                        const mesh::Face& face,
                        const Vector& U) -> double {
    const auto& Sf = face.areaVector();

    if (U.hasFaceValues()) {
        // face values of U are available, no need to manually calculate them
        const auto& Uf = U.valueAtFace(face.id());
        return Uf.dot(Sf);
    }

    const auto& boundary_patch = mesh.boundaryPatch(face);

    if (boundary_patch.isEmpty()) {
        return 0.0;
    }

    const auto& field_bc = boundary_patch.getBoundaryCondition(U.name());
    const auto& Uf = U.valueAtFace(face);
    return Uf.dot(Sf);
}
} // namespace detail

} // namespace prism::ops
