#include <fmt/format.h>

#include "operations.h"
#include "prism/field/field.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"

namespace prism::ops {
auto div_cell(const mesh::PMesh& mesh, const mesh::Cell& cell, const field::Vector& U) -> double;

auto face_flux(const mesh::PMesh& mesh,
               const mesh::Cell& cell,
               const mesh::Face& face,
               const field::Vector& U) -> double;

auto boundary_face_flux(const mesh::PMesh& mesh, const mesh::Face& face, const field::Vector& U)
    -> double;

auto div(const field::Vector& U, bool return_face_data) -> field::Scalar {
    std::string name = fmt::format("div({})", U.name());
    const mesh::PMesh& mesh = U.mesh();

    VectorXd cell_data;
    cell_data.resize(mesh.nCells());

    for (const auto& cell : mesh.cells()) {
        cell_data[cell.id()] = div_cell(mesh, cell, U);
    }

    if (return_face_data) {
        VectorXd face_data;
        face_data.resize(mesh.nFaces());

        // We start with calculating the fluxes at boundary faces
        for (const auto& bface : mesh.boundaryFaces()) {
            face_data[bface.id()] = boundary_face_flux(mesh, bface, U);
        }

        // for interior faces, we take the average of the divergence at the two sharing cells
        for (const auto& iface : mesh.interiorFaces()) {
            const auto& owner = mesh.cell(iface.owner());
            const auto& neighbor = mesh.cell(iface.neighbor().value());
            auto gc = mesh::geo_weight(owner, neighbor, iface);

            auto div_f = gc * cell_data[owner.id()];
            div_f += (1 - gc) * cell_data[neighbor.id()];

            face_data[iface.id()] = div_f;
        }
        return {name, mesh, cell_data, face_data};
    }

    return {name, mesh, cell_data};
}

auto div_cell(const mesh::PMesh& mesh, const mesh::Cell& cell, const field::Vector& U) -> double {
    double sum = 0.0;

    for (auto face_id : cell.faces_ids()) {
        const mesh::Face& face = mesh.face(face_id);
        sum += face_flux(mesh, cell, face, U);
    }

    return sum;
}

auto face_flux(const mesh::PMesh& mesh,
               const mesh::Cell& cell,
               const mesh::Face& face,
               const field::Vector& U) -> double {
    if (face.is_boundary()) {
        return boundary_face_flux(mesh, face, U);
    }

    const Vector3d Uf = U.valueAtFace(face);
    auto Sf = mesh::outward_area_vector(face, cell);
    return Uf.dot(Sf);
}

auto boundary_face_flux(const mesh::PMesh& mesh, const mesh::Face& face, const field::Vector& U)
    -> double {
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

} // namespace prism::ops