#pragma once

#include <prism/prism.h>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <string>
#include <utility>
#include <vector>

#include "prism/field/pressure.h"
#include "prism/field/velocity.h"

namespace prism::test {

using json = nlohmann::json;

/** @brief Loads a mesh and its boundary fields from a UNV mesh file.
 *
 * The boundary fields file is expected to sit next to the mesh file under the name
 * `fields.json`.
 *
 * @param mesh_file Path to the `.unv` mesh file.
 * @return The loaded mesh. */
inline auto loadMesh(const std::filesystem::path& mesh_file) -> SharedPtr<mesh::PMesh> {
    auto boundary_file = mesh_file.parent_path() / "fields.json";
    return mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();
}

/** @brief Returns the relative L2 norm of the difference between two vectors.
 *
 * Computed as ||x - x_ref|| / ||x_ref||.
 *
 * @param x The computed vector.
 * @param x_ref The reference vector.
 * @return The relative L2 norm. */
inline auto l2NormRel(const VectorXd& x, const VectorXd& x_ref) -> f64 {
    return (x - x_ref).norm() / x_ref.norm();
}

/** @brief Computes the volumetric flow rate across a named boundary patch.
 *
 * The rate is the sum of the velocity face values dotted with the face area vectors over
 * all faces belonging to the patch.
 *
 * @param U The velocity field.
 * @param mesh The mesh the field is defined on.
 * @param patch_name Name of the boundary patch to integrate over.
 * @return The total flow rate through the patch (0.0 if the patch does not exist). */
inline auto flowRate(const field::Velocity& U,
                     const mesh::PMesh& mesh,
                     const std::string& patch_name) -> f64 {
    for (const auto& patch : mesh.boundaryPatches()) {
        if (patch.name() != patch_name) {
            continue;
        }
        f64 rate = 0.0;
        for (auto fid : patch.facesIds()) {
            const auto& face = mesh.face(fid);
            rate += U.valueAtFace(face).dot(face.areaVector());
        }
        return rate;
    }
    return 0.0;
}

/** @brief Loads a reference solution from a JSON file.
 *
 * Expects a JSON object with a `p` array of pressure values and a `U` array of velocity
 * vectors, matching the format produced by OpenFOAM export.
 *
 * @param json_path Path to the JSON file.
 * @return A pair of (pressure, velocity) reference values. */
inline auto loadReference(const std::filesystem::path& json_path)
    -> std::pair<VectorXd, std::vector<Vector3d>> {
    std::ifstream file(json_path);
    json data = json::parse(file);

    auto count = data["p"].size();
    VectorXd p_ref(count);
    for (size_t i = 0; i < count; ++i) {
        p_ref(i) = data["p"][i].get<f64>();
    }

    std::vector<Vector3d> U_ref;
    for (const auto& u : data["U"]) {
        U_ref.emplace_back(u[0].get<f64>(), u[1].get<f64>(), u[2].get<f64>());
    }

    return {p_ref, U_ref};
}

/** @brief Returns the relative L2 error of a velocity component against a reference.
 *
 * Extracts the given component (X, Y or Z) from the reference velocity vectors and compares
 * it to the scalar component field.
 *
 * @param component The scalar component field (e.g. `U->x()`).
 * @param U_ref The reference velocity vectors.
 * @param coord Which velocity component to compare.
 * @return The relative L2 error. */
inline auto l2ComponentError(const SharedPtr<field::Scalar>& component,
                             const std::vector<Vector3d>& U_ref,
                             VectorCoord coord) -> f64 {
    VectorXd computed = component->values();
    VectorXd ref(component->mesh()->cellCount());
    for (size_t i = 0; i < component->mesh()->cellCount(); ++i) {
        ref(i) = U_ref[i][static_cast<size_t>(coord)];
    }
    return l2NormRel(computed, ref);
}

/** @brief Builds a scalar field by evaluating a function at every cell center.
 *
 * @tparam Fn Callable taking a cell and returning the cell value.
 * @param mesh The mesh the field is defined on.
 * @param name Name of the created field.
 * @param fn Callable evaluated for each cell to produce its value.
 * @return The constructed scalar field. */
template <typename Fn>
inline auto makeScalarField(const SharedPtr<mesh::PMesh>& mesh, const std::string& name, Fn&& fn)
    -> field::Scalar {
    auto f = std::forward<Fn>(fn);
    VectorXd sol(mesh->cellCount());
    for (const auto& cell : mesh->cells()) {
        sol[cell.id()] = f(cell);
    }
    return {name, mesh, std::move(sol)};
}

/** @brief Solves an equation with the BiCGSTAB solver using non-orthogonal correctors.
 *
 * Runs `n_ortho_correctors` outer solve cycles, each solving the equation with at most
 * `max_iters` iterations down to tolerance `tol`.
 *
 * @tparam Eqn The equation type (e.g. eqn::Transport).
 * @param eqn The equation to solve.
 * @param n_ortho_correctors Number of non-orthogonal corrector cycles.
 * @param max_iters Maximum solver iterations per cycle.
 * @param tol Solver convergence tolerance. */
template <typename Eqn>
inline void solveWithBiCGSTAB(Eqn& eqn,
                              size_t n_ortho_correctors = 5,
                              size_t max_iters = 15,
                              f64 tol = 1e-20) {
    solver::BiCGSTAB solver;
    for (size_t i = 0; i < n_ortho_correctors; ++i) {
        solver.solve(eqn, max_iters, tol);
    }
}

/** @brief Velocity, pressure and mass-flux fields for pressure-based solvers.
 *
 * `mdot` is a velocity field cloned from `U` used to carry the mass flux.
 * @see makePressureVelocityFields */
struct PressureVelocityFields {
    /** @brief The velocity field. */
    SharedPtr<field::Velocity> U;

    /** @brief The pressure field. */
    SharedPtr<field::Pressure> p;

    /** @brief The kinematic viscosity field. */
    SharedPtr<field::Scalar> nu;

    /** @brief The mass-flux velocity field. */
    SharedPtr<field::Velocity> mdot;
};

/** @brief Creates zero-initialised pressure-velocity fields for a mesh.
 *
 * @param mesh The mesh the fields are defined on.
 * @param nu_value Uniform kinematic viscosity value.
 * @return The created fields. */
inline auto makePressureVelocityFields(const SharedPtr<mesh::PMesh>& mesh, f64 nu_value = 1e-3)
    -> PressureVelocityFields {
    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {0, 0, 0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, nu_value);
    auto mdot = makeShared<field::Velocity>(U->clone());
    return {.U = U, .p = p, .nu = nu, .mdot = mdot};
}

/** @brief Steady momentum equations for the X and Y velocity components.
 *
 * Builds both component equations from the same convection, diffusion and gradient schemes,
 * sharing the fields held by @ref PressureVelocityFields. The equations must outlive any
 * span obtained from @ref eqns.
 *
 * @tparam Div The convection scheme type (e.g. scheme::convection::LinearUpwind).
 * @tparam Laplacian The diffusion scheme type (e.g. scheme::diffusion::Corrected).
 * @tparam Grad The pressure gradient scheme type (e.g. source::Gradient<Sign::Negative>).
 * @see makePressureVelocityFields */
template <typename Div, typename Laplacian, typename Grad>
struct MomentumEquations {
    /** @brief The X-component momentum equation. */
    eqn::Momentum u; // NOLINT

    /** @brief The Y-component momentum equation. */
    eqn::Momentum v; // NOLINT

    /** @brief Builds the two momentum equations from the given fields.
     * @param fields The pressure-velocity fields to build the equations from. */
    explicit MomentumEquations(const PressureVelocityFields& fields)
        : u(Div(fields.mdot, fields.U->x()),
            Laplacian(fields.nu, fields.U->x()),
            Grad(fields.p, VectorCoord::X)),

          v(Div(fields.mdot, fields.U->y()),
            Laplacian(fields.nu, fields.U->y()),
            Grad(fields.p, VectorCoord::Y)) {}

    /** @brief Returns pointers to both momentum equations.
     * @return A vector of pointers to the X and Y momentum equations. */
    auto eqns() -> std::vector<eqn::Momentum*> { return {&u, &v}; }
};

} // namespace prism::test
