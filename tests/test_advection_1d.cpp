#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>
#include <memory>

#include "prism/field/scalar.h"
#include "prism/field/velocity.h"

auto advection_analytical_solution(double u, const prism::SharedPtr<prism::mesh::PMesh>& mesh)
    -> prism::field::Scalar {
    prism::VectorXd sol;
    sol.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        double x = cell.center().x();
        double y = cell.center().y();
        double sol_cell = -(std::exp(u * x / 0.1) - 1) / (std::exp(u / 0.1) - 1);
        sol[cell.id()] = sol_cell + 1;
    }

    return {"analytical_solution", mesh, std::move(sol)};
}

TEST_CASE("solve advection equation at u = 2.5 m/s, Pe ~= 5", "[advection]") {
    using namespace prism;
    log::setLevel(log::Level::Error);

    // read mesh
    const auto* mesh_file = "tests/cases/versteeg_advection_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto T = std::make_shared<field::Scalar>("T", mesh, 0.0);

    // set up a uniform velocity field defined over the mesh
    // set the velocity of the field to be the same as the inlet value
    const auto& inlet_patch = std::find_if(
        mesh->boundaryPatches().begin(), mesh->boundaryPatches().end(), [](const auto& patch) {
            return patch.name() == "Inlet";
        });

    // Set a uniform velocity field, with value equal to inlet velocity;
    Vector3d inlet_velocity = inlet_patch->getVectorBoundaryCondition("U");

    auto U = std::make_shared<field::Velocity>("U", mesh, inlet_velocity);
    auto kappa = std::make_shared<field::Scalar>("kappa", mesh, 0.1);

    using div = scheme::convection::Upwind;
    using laplacian = scheme::diffusion::NonCorrected<field::Scalar>;

    auto eqn = eqn::Transport(div(U, T),          // ∇.(ρUT)
                              laplacian(kappa, T) // - ∇.(κ ∇T)
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    auto nOrthogonalCorrectors = 5;

    for (int i = 0; i < nOrthogonalCorrectors; ++i) {
        solver.solve(eqn, 5, 1e-20);
    }

    VectorXd diff = eqn.field()->values().array() -
                    advection_analytical_solution(inlet_velocity.x(), mesh).values().array();
    auto diff_norm = diff.norm();
    REQUIRE(diff_norm < 0.25); // it should be around 0.209

    std::vector<double> T_vec;
    for (const auto& cell : T->mesh()->cells()) {
        T_vec.push_back(T->valueAtCell(cell.id()));
    }

    REQUIRE(std::is_sorted(T_vec.rbegin(), T_vec.rend()));
}
