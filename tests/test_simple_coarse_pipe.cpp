#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <fstream>
#include <memory>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

namespace {

auto loadReference(const std::filesystem::path& json_path)
    -> std::pair<prism::VectorXd, std::vector<prism::Vector3d>> {
    std::ifstream file(json_path);
    json data = json::parse(file);

    auto count = data["p"].size();
    prism::VectorXd p_ref(count);
    for (prism::size_t i = 0; i < count; ++i) {
        p_ref(i) = data["p"][i].get<double>();
    }

    std::vector<prism::Vector3d> U_ref;
    for (const auto& u : data["U"]) {
        U_ref.emplace_back(u[0].get<double>(), u[1].get<double>(), u[2].get<double>());
    }

    return {p_ref, U_ref};
}

auto l2NormRel(const prism::VectorXd& x, const prism::VectorXd& x_ref) -> double {
    return (x - x_ref).norm() / x_ref.norm();
}

} // anonymous namespace

TEST_CASE("SIMPLE algorithm converges on coarse pipe", "[SIMPLESolver]") {
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

    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

    algo::SIMPLEControls controls;
    uEqn.setUnderRelaxFactor(controls.momentum_urf);
    vEqn.setUnderRelaxFactor(controls.momentum_urf);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    auto nOuterIter = 200;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
    }

    auto [p_ref, U_ref] = loadReference(mesh_file.parent_path() / "foam_fields.json");

    // Compare pressure
    double p_error = l2NormRel(p->values(), p_ref);
    INFO("Pressure relative L2 error: " << p_error);
    REQUIRE(p_error < 0.01);

    // Compare Ux component
    VectorXd Ux_computed = U->x()->values();
    VectorXd Ux_ref(mesh->cellCount());
    for (size_t i = 0; i < mesh->cellCount(); ++i) Ux_ref(i) = U_ref[i].x();
    double Ux_error = l2NormRel(Ux_computed, Ux_ref);
    INFO("Ux relative L2 error: " << Ux_error);
    REQUIRE(Ux_error < 0.01);

    // Compare Uy component
    VectorXd Uy_computed = U->y()->values();
    VectorXd Uy_ref(mesh->cellCount());
    for (size_t i = 0; i < mesh->cellCount(); ++i) Uy_ref(i) = U_ref[i].y();
    double Uy_error = l2NormRel(Uy_computed, Uy_ref);
    INFO("Uy relative L2 error: " << Uy_error);
    REQUIRE(Uy_error < 0.1);
}
