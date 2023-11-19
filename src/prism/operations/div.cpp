#include <stdexcept>

#include "operations.h"
#include "prism/field.h"
#include "prism/mesh/boundary.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/mesh/pmesh.h"
#include "prism/mesh/utilities.h"

namespace prism::ops {
auto div_cell(const mesh::PMesh& mesh, const mesh::Cell& cell, const VectorField& U) -> double;

auto face_flux(const mesh::PMesh& mesh,
               const mesh::Cell& cell,
               const mesh::Face& face,
               const VectorField& U) -> double;

auto interior_face_flux(const mesh::PMesh& mesh,
                        const mesh::Cell& cell,
                        const mesh::Face& face,
                        const VectorField& U) -> double;

auto boundary_face_flux(const mesh::PMesh& mesh,
                        const mesh::Cell& cell,
                        const mesh::Face& face,
                        const VectorField& U) -> double;

auto div(const VectorField& U, bool return_face_data) -> ScalarField {
    // There is a room for optimization here, by first calculating the face fluxes
    // and store the result in a vector x, and provide this vector to each call of div_cell()
    // to avoid re-calling face_flux() for every cell, and avoid calling face_flux() again
    // when we update face_data vector, which in this case will be the vector x.
    std::string name = fmt::format("div({})", U.name());
    const mesh::PMesh& mesh = U.mesh();

    VectorXd data;
    data.resize(mesh.n_cells());

    for (const auto& cell : mesh.cells()) {
        data[cell.id()] = div_cell(mesh, cell, U);
    }

    if (return_face_data) {
        VectorXd face_data;
        face_data.resize(mesh.faces().size());

        for (const auto& face : mesh.faces()) {
            const auto& owner = mesh.cell(face.owner());
            face_data[face.id()] = face_flux(mesh, owner, face, U);
        }

        return {name, mesh, data, face_data};
    }

    return {name, mesh, data};
}

auto div_cell(const mesh::PMesh& mesh, const mesh::Cell& cell, const VectorField& U) -> double {
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
               const VectorField& U) -> double {
    if (face.is_boundary()) {
        return boundary_face_flux(mesh, cell, face, U);
    }
    return interior_face_flux(mesh, cell, face, U);
}

auto interior_face_flux(const mesh::PMesh& mesh,
                        const mesh::Cell& cell,
                        const mesh::Face& face,
                        const VectorField& U) -> double {
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neigbor = mesh.cell(face.neighbor().value());

    double gc = mesh::geo_weight(owner, neigbor, face);

    Vector3d Uf = gc * U[owner.id()];
    Uf += (1 - gc) * U[neigbor.id()];

    auto Sf = mesh::outward_area_vector(face, cell);

    return Uf.dot(Sf);
}

auto boundary_face_flux(const mesh::PMesh& mesh,
                        const mesh::Cell& cell,
                        const mesh::Face& face,
                        const VectorField& U) -> double {
    const auto& Sf = face.area_vector();

    if (U.has_face_data()) {
        // face values of U are available, no need to manually calculate them
        const auto& Uf = U.value_at_face(face.id());
        return Uf.dot(Sf);
    }

    const auto& boundary_patch = mesh.boundary_patch(face);
    const auto& field_boundary_type = boundary_patch.get_bc(U.name());

    switch (field_boundary_type.bc_type()) {
        case mesh::BoundaryConditionType::Empty: {
            return 0.0;
        }

        case mesh::BoundaryConditionType::Fixed:
        case mesh::BoundaryConditionType::Inlet: {
            const auto& Uf = boundary_patch.get_vector_bc(U.name());
            return Uf.dot(Sf);
        }

        case mesh::BoundaryConditionType::Outlet:
        case mesh::BoundaryConditionType::Symmetry: {
            const auto& Uf = U[cell.id()];
            return Uf.dot(Sf);
        }

        default:
            throw std::runtime_error(
                // TODO: write better error message
                "prism::ops::boundary_face_flux() was given a non-implemented boundary "
                "condition");
    }
    return 0.0;
}

} // namespace prism::ops