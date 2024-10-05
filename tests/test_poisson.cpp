#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>

auto solution(const auto& mesh) -> prism::field::Scalar {
    prism::VectorXd sol;
    sol.resize(mesh.cellCount());

    for (const auto& cell : mesh.cells()) {
        double x = cell.center().x();
        double y = cell.center().y();
        double sol_cell = std::sin(prism::PI * x);
        sol_cell *= std::cos(prism::PI * y);
        sol[cell.id()] = sol_cell;
    }

    return prism::field::Scalar("S", mesh, std::move(sol));
}

TEST_CASE("test poisson equation", "[poisson]") {
    using namespace prism;
    using namespace prism::scheme;


    const auto* unv_file_name = "tests/cases/poisson/mesh.unv";

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    auto P = field::Scalar("P", mesh, 0.0);

    // create source term
    VectorXd src_values;
    src_values.resize(mesh.cellCount());

    for (const auto& cell : mesh.cells()) {
        double x = cell.center().x();
        double y = cell.center().y();
        double src = 2 * std::pow(prism::PI, 2);
        src *= std::sin(prism::PI * x);
        src *= std::cos(prism::PI * y);
        src_values[cell.id()] = src;
    }

    auto c = field::Tensor("c", mesh, Matrix3d::Identity());

    using laplacian =
        diffusion::Corrected<field::Tensor, nonortho::OverRelaxedCorrector, field::Scalar>;
    auto source = field::Scalar("S", mesh, std::move(src_values));

    auto eqn = eqn::Transport<field::Scalar>(
        laplacian(c, P),                                                            // -∇.∇p
        source::ConstantScalar<source::SourceSign::Positive, field::Scalar>(source) // = S
    );

    // solve
    auto solver =
        solver::BiCGSTAB<field::Scalar, solver::ImplicitUnderRelaxation<field::Scalar>>();
    solver.solve(eqn, 25, 1e-20, 1.0);

    VectorXd diff = P.values() - solution(mesh).values();
    double diff_norm = diff.norm();

    // TODO: this is a large l2-norm criteria, we need to check if the solution is correct
    REQUIRE(diff_norm < 0.1);
}
