#pragma once

#include <fmt/format.h>

#include <algorithm>

#include "prism/field/scalar.h"
#include "prism/field/vector.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
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

namespace detail {
auto inline coordToIndex(Coord coord) -> std::uint8_t {
    switch (coord) {
        case Coord::X: return 0;
        case Coord::Y: return 1;
        case Coord::Z: return 2;
        default: break;
    }
    throw std::invalid_argument("Invalid Coord value in coordToIndex");
}

template <typename Vector>
auto fluxSumAtCell(const mesh::Cell& cell, const Vector& U) -> double;
} // namespace detail

template <typename Vector>
auto div(const Vector& U) -> field::Scalar {
    // applying green-gauss theorem to the divergence of a vector field
    // ∫(∇.U) dV = ∫U.dS
    // (∇.U) V = Σ U.S
    // ∇.U = (Σ U.S)/V = ops::div(U)
    std::string name = fmt::format("div({})", U.name());
    const auto& mesh = U.mesh();

    VectorXd cell_data;
    cell_data.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        cell_data[cell.id()] = detail::fluxSumAtCell(cell, U) / cell.volume();
    }
    return {name, mesh, cell_data};
}

template <field::IScalarBased Field>
auto grad(const Field& field, Coord coord) -> field::Scalar {
    auto grad_field_name = fmt::format("grad({})_{}", field.name(), field::coordToStr(coord));
    const auto& mesh = field.mesh();
    VectorXd grad_values = VectorXd::Zero(mesh->cellCount());
    auto coord_index = detail::coordToIndex(coord);

    std::for_each(mesh->cells().begin(), mesh->cells().end(), [&](const auto& cell) {
        grad_values[cell.id()] = field.gradAtCell(cell)[coord_index];
    });

    return field::Scalar(grad_field_name, field.mesh(), grad_values);
}

template <field::IScalarBased Field>
auto grad(const Field& field) -> field::Vector {
    std::array<field::Scalar, 3> fields {
        grad(field, Coord::X), grad(field, Coord::Y), grad(field, Coord::Z)};

    return field::Vector(fmt::format("grad({})", field.name()), field.mesh(), fields);
}

namespace detail {
template <typename Vector>
auto fluxSumAtCell(const mesh::Cell& cell, const Vector& U) -> double {
    double sum = 0.0;
    const auto& mesh = U.mesh();

    for (auto face_id : cell.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);

        // Skip empty faces
        if (face.isBoundary()) {
            const auto& boundary_patch = mesh->boundaryPatch(face);
            if (boundary_patch.isEmpty()) {
                continue;
            }
        }
        double flux = U.fluxAtFace(face);
        if (!face.isOwnedBy(cell.id())) {
            // fluxAtFace returns flux in face normal direction from owner to face center, we need
            // to reverse this value when `face` is not owned by `cell`
            flux = -flux;
        }
        sum += flux;
    }

    return sum;
}

} // namespace detail

} // namespace prism::ops
