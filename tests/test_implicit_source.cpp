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

using namespace prism;
using namespace prism::scheme;

auto implicit_analytic_solution(const auto& mesh) -> prism::field::Scalar {
    // y = sinh(2x + 1)
    /// TODO: this can be done more efficiently using a field::updateCells() function
    prism::VectorXd sol;
    sol.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        double x = cell.center().x();
        double sol_cell = std::sinh(2 * x + 1);
        sol[cell.id()] = sol_cell;
    }

    return {"S", mesh, std::move(sol)};
}

auto l2NormRelative(const Vector3d& x, const Vector3d& x_ref) -> double {
    return (x - x_ref).norm() / x_ref.norm();
}

TEST_CASE("test implicit source", "[implicit-source]") {
    const auto* unv_file_name = "tests/cases/channel1d_coarse/mesh.unv";

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    auto y = field::Scalar("y", mesh, 0.0);
    auto c = field::UniformScalar("c", mesh, 1.0);

    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport<field::Scalar>(
        laplacian(c, y),                                           // -∇.∇y
        source::ImplicitField<Sign::Negative, field::Scalar>(4, y) // = -4y
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    auto nOuterIter = 5;

    for (int i = 0; i < nOuterIter; ++i) {
        solver.solve(eqn, 15, 1e-20);
    }
    auto norm = l2NormRelative(y.values(), implicit_analytic_solution(mesh).values());
    REQUIRE(norm < 0.004);
}
