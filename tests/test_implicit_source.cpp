/*
 *
 * This test case is based on Numerical Analysis with Applications in Python book
 * by John S Butler (john.s.butler@tudublin.ie)
 * https://john-s-butler-dit.github.io/NumericalAnalysisBook/index.html
 *
 */
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

auto implicit_analytic_solution(const SharedPtr<mesh::PMesh>& mesh) -> field::Scalar {
    // y = sinh(2x + 1)
    return makeScalarField(mesh, "S", [](const auto& cell) {
        double x = cell.center().x();
        return std::sinh((2 * x) + 1);
    });
}

} // anonymous namespace

TEST_CASE("test implicit source", "[implicit-source]") {
    auto mesh = loadMesh("tests/cases/channel1d_coarse/mesh.unv");

    auto y = makeShared<field::Scalar>("y", mesh, 0.0);
    auto c = makeShared<field::Scalar>("c", mesh, 1.0);

    using laplacian = diffusion::NonCorrected;

    auto eqn = eqn::Transport(laplacian(c, y),                            // -∇.∇y
                              source::ImplicitField<Sign::Negative>(4, y) // = -4y
    );

    solveWithBiCGSTAB(eqn);

    auto norm = l2NormRel(y->values(), implicit_analytic_solution(mesh).values());
    REQUIRE(norm < 0.005);
}
