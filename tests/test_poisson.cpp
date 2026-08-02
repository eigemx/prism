#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>
#include <memory>

#include "test_utils.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::test;

namespace {

auto poisson_analytical_solution(const auto& mesh) -> field::Scalar {
    return makeScalarField(mesh, "S", [](const auto& cell) {
        double x = cell.center().x();
        double y = cell.center().y();
        return std::sin(prism::PI * x) * std::cos(prism::PI * y);
    });
}

auto testPoissonWithMesh(const std::string& mesh_file_name) -> SharedPtr<field::Scalar> {
    auto mesh = loadMesh(mesh_file_name);

    auto P = makeShared<field::Scalar>("P", mesh, 0.0);

    VectorXd src_values =
        makeScalarField(mesh, "S", [](const auto& cell) {
            double x = cell.center().x();
            double y = cell.center().y();
            return 2 * std::pow(prism::PI, 2) * std::sin(prism::PI * x) * std::cos(prism::PI * y);
        }).values();

    auto c = makeShared<field::Tensor>("c", mesh, Matrix3d::Identity());

    using laplacian = diffusion::Corrected;
    auto source = makeShared<field::Scalar>("S", mesh, std::move(src_values));

    auto eqn = eqn::Transport(laplacian(c, P),                               // -∇.∇p
                              source::ConstantScalar<Sign::Positive>(source) // = S
    );

    solveWithBiCGSTAB(eqn);

    return P;
}

} // anonymous namespace

TEST_CASE("test poisson equation unstructured", "[poisson]") {
    auto P = testPoissonWithMesh("tests/cases/poisson/mesh.unv");
    double diff_norm = l2NormRel(P->values(), poisson_analytical_solution(P->mesh()).values());
    REQUIRE(diff_norm < 0.078); // for poisson/mesh.unv it should be = 0.0769
}

TEST_CASE("test poisson equation structured", "[poisson]") {
    auto P = testPoissonWithMesh("tests/cases/poisson/mesh_hex.unv");
    double diff_norm = l2NormRel(P->values(), poisson_analytical_solution(P->mesh()).values());
    REQUIRE(diff_norm < 0.0004); // for poisson/mesh_hex.unv it should be = 0.000323
}
