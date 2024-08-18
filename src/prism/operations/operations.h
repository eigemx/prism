#pragma once

#include <fmt/format.h>

#include "prism/field/scalar.h"
#include "prism/field/vector.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"

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

// face mass flow rate
auto inline faceFlowRate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}

namespace detail {
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
    // this is a boundary face, where normal is always pointing outside of the cell
    // no need to call mesh::outward_area_vector()
    const auto& Sf = face.area_vector();

    if (U.hasFaceValues()) {
        // face values of U are available, no need to manually calculate them
        const auto& Uf = U.valueAtFace(face.id());
        return Uf.dot(Sf);
    }

    const auto& boundary_patch = mesh.boundary_patch(face);

    if (boundary_patch.isEmpty()) {
        return 0.0;
    }

    const auto& field_bc = boundary_patch.getBoundaryCondition(U.name());
    const auto& Uf = U.valueAtFace(face);
    return Uf.dot(Sf);
}
} // namespace detail

} // namespace prism::ops
