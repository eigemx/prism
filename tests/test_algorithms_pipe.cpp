#include <prism/algorithm/PISO.h>
#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "test_utils.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::scheme::convection;
using namespace prism::test;

TEST_CASE("SIMPLE algorithm mass conservation on coarse pipe", "[algo]") {
    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto mesh = loadMesh(mesh_file);

    auto fields = makePressureVelocityFields(mesh);

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    MomentumEquations<div, laplacian, grad> meq(fields);

    algo::SIMPLEControls controls;

    auto momentum_eqns = meq.eqns();
    auto nOuterIter = 200;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), fields.U, fields.mdot, fields.p);
    }

    auto inlet = flowRate(*fields.U, *mesh, "Inlet");
    auto outlet = flowRate(*fields.U, *mesh, "Outlet");
    f64 imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << ", imbalance = " << imbalance);
    REQUIRE(imbalance < 1e-6);
}

TEST_CASE("PISO auto pRef on cavity_hex (closed domain)", "[algo]") {
    auto mesh_file = std::filesystem::path("tests/cases/cavity_hex/mesh.unv");
    auto mesh = loadMesh(mesh_file);

    auto fields = makePressureVelocityFields(mesh, 0.01);

    using div = LinearUpwind;
    using lap = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    fields.U->x()->setHistorySize(2);
    fields.U->y()->setHistorySize(2);
    fields.U->x()->updatePrevTimeSteps();
    fields.U->y()->updatePrevTimeSteps();

    f64 dt = 0.005;
    auto uEqn = eqn::Momentum(temporal::BackwardEuler(fields.U->x(), dt),
                              div(fields.mdot, fields.U->x()),
                              lap(fields.nu, fields.U->x()),
                              grad(fields.p, VectorCoord::X));
    auto vEqn = eqn::Momentum(temporal::BackwardEuler(fields.U->y(), dt),
                              div(fields.mdot, fields.U->y()),
                              lap(fields.nu, fields.U->y()),
                              grad(fields.p, VectorCoord::Y));

    std::vector<eqn::Momentum*> meq {&uEqn, &vEqn};

    algo::PISOControls controls;
    controls.outer_iterations = 3;
    controls.pressure_correction_steps = 1;
    controls.non_ortho_correctors = 0;

    auto nSteps = 20;
    for (int ts = 0; ts < nSteps; ++ts) {
        algo::PISO(controls).step(std::span<eqn::Momentum*>(meq), fields.U, fields.mdot, fields.p);
        fields.U->x()->updatePrevTimeSteps();
        fields.U->y()->updatePrevTimeSteps();
    }
    REQUIRE_FALSE(std::isnan(fields.p->values().lpNorm<Eigen::Infinity>()));
    REQUIRE_FALSE(std::isinf(fields.p->values().lpNorm<Eigen::Infinity>()));
}

TEST_CASE("SIMPLE ductSIMPLE (open domain, no pRef needed)", "[algo]") {
    auto mesh_file = std::filesystem::path("tests/cases/ductSIMPLE/mesh.unv");
    auto mesh = loadMesh(mesh_file);

    auto fields = makePressureVelocityFields(mesh);

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    MomentumEquations<div, laplacian, grad> meq(fields);

    meq.u.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Transport>>();
    meq.v.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Transport>>();

    algo::SIMPLEControls controls;

    auto momentum_eqns = meq.eqns();
    auto nOuterIter = 500;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), fields.U, fields.mdot, fields.p);
    }

    auto inlet = flowRate(*fields.U, *mesh, "inlet");
    auto outlet = flowRate(*fields.U, *mesh, "outlet");
    f64 imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << "  imbalance = " << imbalance);
    REQUIRE(imbalance < 1e-6);
}

TEST_CASE("PISO no-PRIME mass conservation on coarse pipe", "[algo]") {
    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto mesh = loadMesh(mesh_file);

    auto fields = makePressureVelocityFields(mesh);

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    MomentumEquations<div, laplacian, grad> meq(fields);

    algo::PISOControls controls;
    controls.outer_iterations = 1;
    controls.momentum_implicit_steps = 1;
    controls.pressure_correction_steps = 1;
    controls.non_ortho_correctors = 2;
    // Steady pipe (no time term): steady pressure-based solvers need under-relaxation to converge.
    controls.momentum_urf = 0.7;
    controls.pressure_urf = 0.3;
    controls.momentum_residual = 1e-7;
    controls.pressure_residual = 1e-7;

    auto momentum_eqns = meq.eqns();
    auto nOuterIter = 200;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::PISO(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), fields.U, fields.mdot, fields.p);
    }

    auto inlet = flowRate(*fields.U, *mesh, "Inlet");
    auto outlet = flowRate(*fields.U, *mesh, "Outlet");
    f64 imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << ", imbalance = " << imbalance);
    REQUIRE(imbalance < 1e-6);
}

TEST_CASE("PISO with PRIME mass conservation on coarse pipe", "[algo]") {
    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto mesh = loadMesh(mesh_file);

    auto fields = makePressureVelocityFields(mesh);

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    MomentumEquations<div, laplacian, grad> meq(fields);

    algo::PISOControls controls;
    controls.outer_iterations = 1;
    controls.momentum_implicit_steps = 1;
    controls.pressure_correction_steps = 2;
    controls.non_ortho_correctors = 2;
    // Steady pipe (no time term): steady pressure-based solvers need under-relaxation to converge.
    controls.momentum_urf = 0.7;
    controls.pressure_urf = 0.3;
    controls.momentum_residual = 1e-7;
    controls.pressure_residual = 1e-7;

    auto momentum_eqns = meq.eqns();
    auto nSteps = 50;
    for (int i = 0; i < nSteps; ++i) {
        algo::PISO(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), fields.U, fields.mdot, fields.p);
    }

    auto inlet = flowRate(*fields.U, *mesh, "Inlet");
    auto outlet = flowRate(*fields.U, *mesh, "Outlet");
    f64 imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << ", imbalance = " << imbalance);
    REQUIRE(imbalance < 1e-6);
}
