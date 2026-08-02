#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>
#include <memory>

#include "test_utils.h"

using namespace prism;
using namespace prism::test;

namespace {

auto advection_analytical_solution(f64 u, const SharedPtr<mesh::PMesh>& mesh) -> field::Scalar {
    return makeScalarField(mesh, "analytical_solution", [u](const auto& cell) {
        double x = cell.center().x();
        return -((std::exp(u * x / 0.1) - 1) / (std::exp(u / 0.1) - 1)) + 1;
    });
}

} // anonymous namespace

TEST_CASE("solve advection equation at u = 2.5 m/s, Pe ~= 5", "[advection]") {
    using namespace prism::scheme;
    log::setLevel(log::Level::Error);

    auto mesh = loadMesh("tests/cases/versteeg_advection_1d/mesh.unv");

    auto T = makeShared<field::Scalar>("T", mesh, 0.0);

    // set up a uniform velocity field defined over the mesh
    // set the velocity of the field to be the same as the inlet value
    const auto& inlet_patch = std::find_if(
        mesh->boundaryPatches().begin(), mesh->boundaryPatches().end(), [](const auto& patch) {
            return patch.name() == "Inlet";
        });

    // Set a uniform velocity field, with value equal to inlet velocity;
    Vector3d inlet_velocity = inlet_patch->getVectorBoundaryCondition("U");

    auto U = makeShared<field::Velocity>("U", mesh, inlet_velocity);
    auto kappa = makeShared<field::Scalar>("kappa", mesh, 0.1);

    using div = scheme::convection::Upwind;
    using laplacian = scheme::diffusion::NonCorrected;

    auto eqn = eqn::Transport(div(U, T),          // ∇.(ρUT)
                              laplacian(kappa, T) // - ∇.(κ ∇T)
    );

    solveWithBiCGSTAB(eqn, 5, 5, 1e-20);

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
