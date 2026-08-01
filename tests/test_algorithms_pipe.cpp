#include <prism/algorithm/PISO.h>
#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <memory>

namespace {

auto flowRate(const prism::field::Velocity& U,
              const prism::mesh::PMesh& mesh,
              const std::string& patch_name) -> double {
    for (const auto& patch : mesh.boundaryPatches()) {
        if (patch.name() != patch_name) continue;
        double rate = 0.0;
        for (auto fid : patch.facesIds()) {
            const auto& face = mesh.face(fid);
            rate += U.valueAtFace(face).dot(face.areaVector());
        }
        return rate;
    }
    return 0.0;
}

} // anonymous namespace

TEST_CASE("SIMPLE algorithm mass conservation on coarse pipe", "[algo]") {
    using namespace prism;
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";

    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {0, 0, 0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 1e-3);
    auto mdot = makeShared<field::Velocity>(U->clone());

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    auto uEqn = eqn::Momentum(div(mdot, U->x()), laplacian(nu, U->x()), grad(p, VectorCoord::X));
    auto vEqn = eqn::Momentum(div(mdot, U->y()), laplacian(nu, U->y()), grad(p, VectorCoord::Y));

    algo::SIMPLEControls controls;
    uEqn.setUnderRelaxFactor(controls.momentum_urf);
    vEqn.setUnderRelaxFactor(controls.momentum_urf);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    auto nOuterIter = 200;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
    }

    auto inlet = flowRate(*U, *mesh, "Inlet");
    auto outlet = flowRate(*U, *mesh, "Outlet");
    double imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << ", imbalance = " << imbalance);
    REQUIRE(imbalance < 0.05);
}

TEST_CASE("PISO auto pRef on cavity_hex (closed domain)", "[algo]") {
    using namespace prism;
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    auto mesh_file = std::filesystem::path("tests/cases/cavity_hex/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";

    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {0, 0, 0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 0.01);
    auto mdot = makeShared<field::Velocity>(U->clone());

    using div = LinearUpwind;
    using lap = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    U->x()->setHistorySize(2);
    U->y()->setHistorySize(2);
    U->x()->updatePrevTimeSteps();
    U->y()->updatePrevTimeSteps();

    f64 dt = 0.005;
    auto uEqn = eqn::Momentum(temporal::BackwardEuler(U->x(), dt),
                              div(mdot, U->x()),
                              lap(nu, U->x()),
                              grad(p, VectorCoord::X));
    auto vEqn = eqn::Momentum(temporal::BackwardEuler(U->y(), dt),
                              div(mdot, U->y()),
                              lap(nu, U->y()),
                              grad(p, VectorCoord::Y));

    uEqn.setUnderRelaxFactor(0.7);
    vEqn.setUnderRelaxFactor(0.7);
    std::vector<eqn::Momentum*> meq {&uEqn, &vEqn};

    algo::PISOControls controls;
    controls.outer_iterations = 3;
    controls.pressure_correction_steps = 1;
    controls.non_ortho_correctors = 0;
    controls.momentum_urf = 0.7;
    controls.pressure_urf = 0.3;

    auto nSteps = 20;
    for (int ts = 0; ts < nSteps; ++ts) {
        algo::PISO(controls).step(std::span<eqn::Momentum*>(meq), U, mdot, p);
        U->x()->updatePrevTimeSteps();
        U->y()->updatePrevTimeSteps();
    }
    REQUIRE_FALSE(std::isnan(p->values().lpNorm<Eigen::Infinity>()));
    REQUIRE_FALSE(std::isinf(p->values().lpNorm<Eigen::Infinity>()));
}

TEST_CASE("SIMPLE ductSIMPLE (open domain, no pRef needed)", "[algo]") {
    using namespace prism;
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    auto mesh_file = std::filesystem::path("tests/cases/ductSIMPLE/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";

    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {0, 0, 0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 1e-3);
    auto mdot = makeShared<field::Velocity>(U->clone());

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    auto uEqn = eqn::Momentum(div(mdot, U->x()), laplacian(nu, U->x()), grad(p, VectorCoord::X));
    auto vEqn = eqn::Momentum(div(mdot, U->y()), laplacian(nu, U->y()), grad(p, VectorCoord::Y));

    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Transport>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Transport>>();

    algo::SIMPLEControls controls;
    uEqn.setUnderRelaxFactor(controls.momentum_urf);
    vEqn.setUnderRelaxFactor(controls.momentum_urf);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    auto nOuterIter = 500;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
    }

    auto inlet = flowRate(*U, *mesh, "inlet");
    auto outlet = flowRate(*U, *mesh, "outlet");
    double imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << "  imbalance = " << imbalance);
    REQUIRE(imbalance < 0.05);
}

TEST_CASE("PISO no-PRIME mass conservation on coarse pipe", "[algo]") {
    using namespace prism;
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";

    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {0, 0, 0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 1e-3);
    auto mdot = makeShared<field::Velocity>(U->clone());

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    auto uEqn = eqn::Momentum(div(mdot, U->x()), laplacian(nu, U->x()), grad(p, VectorCoord::X));
    auto vEqn = eqn::Momentum(div(mdot, U->y()), laplacian(nu, U->y()), grad(p, VectorCoord::Y));

    uEqn.setUnderRelaxFactor(0.7);
    vEqn.setUnderRelaxFactor(0.7);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    algo::PISOControls controls;
    controls.outer_iterations = 1;
    controls.momentum_implicit_steps = 1;
    controls.pressure_correction_steps = 1;
    controls.non_ortho_correctors = 2;
    controls.momentum_urf = 0.7;
    controls.pressure_urf = 0.3;
    controls.momentum_residual = 1e-7;
    controls.pressure_residual = 1e-7;

    auto nOuterIter = 200;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::PISO(controls).step(std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
    }

    auto inlet = flowRate(*U, *mesh, "Inlet");
    auto outlet = flowRate(*U, *mesh, "Outlet");
    double imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << ", imbalance = " << imbalance);
    REQUIRE(imbalance < 0.05);
}

TEST_CASE("PISO with PRIME mass conservation on coarse pipe", "[algo]") {
    using namespace prism;
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";

    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {0, 0, 0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 1e-3);
    auto mdot = makeShared<field::Velocity>(U->clone());

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    auto uEqn = eqn::Momentum(div(mdot, U->x()), laplacian(nu, U->x()), grad(p, VectorCoord::X));
    auto vEqn = eqn::Momentum(div(mdot, U->y()), laplacian(nu, U->y()), grad(p, VectorCoord::Y));

    uEqn.setUnderRelaxFactor(0.7);
    vEqn.setUnderRelaxFactor(0.7);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    algo::PISOControls controls;
    controls.outer_iterations = 1;
    controls.momentum_implicit_steps = 1;
    controls.pressure_correction_steps = 2;
    controls.non_ortho_correctors = 2;
    controls.momentum_urf = 0.7;
    controls.pressure_urf = 0.3;
    controls.momentum_residual = 1e-7;
    controls.pressure_residual = 1e-7;

    auto nSteps = 50;
    for (int i = 0; i < nSteps; ++i) {
        algo::PISO(controls).step(std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
    }

    auto inlet = flowRate(*U, *mesh, "Inlet");
    auto outlet = flowRate(*U, *mesh, "Outlet");
    double imbalance = std::abs(inlet + outlet) / std::abs(inlet);
    INFO("inlet flow = " << inlet << "  outlet flow = " << outlet << ", imbalance = " << imbalance);
    REQUIRE(imbalance < 0.05);
}
