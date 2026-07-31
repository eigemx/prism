#pragma once

#include <fmt/format.h>

#include "prism/field/scalar.h"
#include "prism/field/tensor.h"
#include "prism/field/vector.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"
#include "prism/types.h"

namespace prism::ops {

/** @brief Divergence of a vector field: \f$ \nabla \cdot \mathbf{U} \f$.
 *
 * Computed via the Green-Gauss theorem:
 * \f$ \int (\nabla \cdot \mathbf{U}) dV = \oint \mathbf{U} \cdot d\mathbf{S} \f$
 *
 * @param U Vector field (must support fluxAtFace).
 * @return Scalar field \f$ \nabla \cdot \mathbf{U} \f$. */
template <typename Vector>
auto div(const Vector& U) -> field::Scalar;

/** @brief Single component of the gradient of a scalar field.
 *
 * Extracts one coordinate (X, Y, or Z) of \f$ \nabla \phi \f$ as a Scalar field.
 *
 * @param field Scalar field whose gradient is computed.
 * @param coord Which component to extract (VectorCoord::X / Y / Z).
 * @return Scalar field containing \f$ \partial \phi / \partial x_i \f$. */
template <field::IScalarBased Field>
auto grad(Field& field, VectorCoord coord) -> field::Scalar;

/** @brief Gradient of a scalar field as a Vector field.
 *
 * Combines all three components of \f$ \nabla \phi \f$ into a Vector field.
 * Equivalent to calling grad(field, X), grad(field, Y), grad(field, Z).
 *
 * @param field Scalar field.
 * @return Vector field \f$ \nabla \phi \f$. */
template <field::IScalarBased Field>
auto grad(Field& field) -> field::Vector;

/** @brief Velocity gradient tensor as a full Tensor field.
 *
 * One 3×3 matrix per cell. Equivalent to calling velocityGradientAtCell()
 * on every cell in the mesh.
 *
 * @param U Vector field.
 * @return Tensor field \f$ \nabla \mathbf{U} \f$. */
template <typename VectorType>
auto grad(const VectorType& U) -> field::Tensor;

/** @brief Symmetric part of a tensor: \f$ \frac{1}{2}(A + A^T) \f$.
 *
 * @param A The tensor.
 * @return The symmetric part of A. */
inline auto symm(const Tensor3d& A) -> Tensor3d;

/** @brief Twice the symmetric part of a tensor: \f$ A + A^T \f$.
 *
 * @param A The tensor.
 * @return A + A^T. */
inline auto twoSymm(const Tensor3d& A) -> Tensor3d;

/** @brief Deviatoric part of a tensor: \f$ A - \frac{tr(A)}{3} I \f$.
 *
 * @param A The tensor.
 * @return The deviatoric part of A. */
inline auto dev(const Tensor3d& A) -> Tensor3d;

/** @brief Double contraction of two tensors: \f$ A : B = \sum_{ij} A_{ij} B_{ij} \f$.
 *
 * @param A First tensor.
 * @param B Second tensor.
 * @return The double contraction of A and B. */
inline auto doubleContraction(const Tensor3d& A, const Tensor3d& B) -> f64;

namespace detail {
template <typename Vector>
auto fluxSumAtCell(const mesh::Cell& cell, const Vector& U) -> f64;
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
auto grad(Field& field, VectorCoord coord) -> field::Scalar {
    const auto& mesh = field.mesh();
    VectorXd grad_values(mesh->cellCount());
    auto coord_index = field::detail::coordToIndex(coord);

    for (const auto& cell : mesh->cells()) {
        grad_values[cell.id()] = field.gradAtCell(cell)[coord_index];
    }

    return field::Scalar(fmt::format("grad({})_{}", field.name(), field::detail::coordToStr(coord)),
                         mesh,
                         std::move(grad_values));
}

template <field::IScalarBased Field>
auto grad(Field& field) -> field::Vector {
    const auto& mesh = field.mesh();
    const auto n = mesh->cellCount();
    VectorXd x_values(n);
    VectorXd y_values(n);
    VectorXd z_values(n);

    for (const auto& cell : mesh->cells()) {
        Vector3d g = field.gradAtCell(cell);
        auto i = cell.id();
        x_values[i] = g.x();
        y_values[i] = g.y();
        z_values[i] = g.z();
    }

    auto name = fmt::format("grad({})", field.name());
    std::array<field::Scalar, 3> components {
        field::Scalar(name + "_x", mesh, std::move(x_values)),
        field::Scalar(name + "_y", mesh, std::move(y_values)),
        field::Scalar(name + "_z", mesh, std::move(z_values)),
    };

    return field::Vector(name, mesh, components);
}

namespace detail {
template <typename Vector>
auto fluxSumAtCell(const mesh::Cell& cell, const Vector& U) -> f64 {
    f64 sum = 0.0;
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
        f64 flux = U.fluxAtFace(face);
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

template <typename VectorType>
auto grad(const VectorType& U) -> field::Tensor {
    const auto& mesh = U.mesh();
    std::vector<Tensor3d> data(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        data[cell.id()] = U.gradAtCell(cell);
    }

    return field::Tensor(fmt::format("grad({})", U.name()), mesh, std::move(data));
}

inline auto symm(const Tensor3d& A) -> Tensor3d {
    return 0.5 * (A + A.transpose());
}

inline auto twoSymm(const Tensor3d& A) -> Tensor3d {
    return A + A.transpose();
}

inline auto dev(const Tensor3d& A) -> Tensor3d {
    return A - Tensor3d::Identity() * (A.trace() / 3.0);
}

inline auto doubleContraction(const Tensor3d& A, const Tensor3d& B) -> f64 {
    return A.cwiseProduct(B).sum();
}

} // namespace prism::ops
