#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "test_utils.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::scheme::convection;
using namespace prism::test;

TEST_CASE("SIMPLE algorithm converges on coarse pipe", "[SIMPLESolver]") {
    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHex/mesh.unv");
    auto mesh = loadMesh(mesh_file);

    auto fields = makePressureVelocityFields(mesh);

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    MomentumEquations<div, laplacian, grad> meq(fields);

    algo::SIMPLEControls controls;

    auto momentum_eqns = meq.eqns();
    auto nOuterIter = 100;
    for (int iter = 0; iter < nOuterIter; ++iter) {
        algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), fields.U, fields.mdot, fields.p);
    }

    auto [p_ref, U_ref] = loadReference(mesh_file.parent_path() / "foam_fields.json");

    // Compare pressure
    double p_error = l2NormRel(fields.p->values(), p_ref);
    INFO("Pressure relative L2 error: " << p_error);
    REQUIRE(p_error < 0.005);

    // Compare Ux component
    double Ux_error = l2ComponentError(fields.U->x(), U_ref, VectorCoord::X);
    INFO("Ux relative L2 error: " << Ux_error);
    REQUIRE(Ux_error < 0.005);

    // Compare Uy component
    double Uy_error = l2ComponentError(fields.U->y(), U_ref, VectorCoord::Y);
    INFO("Uy relative L2 error: " << Uy_error);
    REQUIRE(Uy_error < 0.09);
}
