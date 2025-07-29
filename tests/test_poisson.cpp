#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>

auto poisson_analytical_solution(const auto& mesh) -> prism::field::Scalar {
    prism::VectorXd sol;
    sol.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        double x = cell.center().x();
        double y = cell.center().y();
        double sol_cell = std::sin(prism::PI * x);
        sol_cell *= std::cos(prism::PI * y);
        sol[cell.id()] = sol_cell;
    }

    return prism::field::Scalar("S", mesh, std::move(sol));
}

auto l2NormRel(const prism::Vector3d& x, const prism::Vector3d& x_ref) -> double {
    return (x - x_ref).norm() / x_ref.norm();
}

auto testPoissonWithMesh(const std::string& mesh_file_name) -> prism::field::Scalar {
    using namespace prism;
    using namespace prism::scheme;

    // read mesh
    auto boundary_file = std::filesystem::path(mesh_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file_name, boundary_file).toPMesh();

    auto P = field::Scalar("P", mesh, 0.0);

    // create source term
    VectorXd src_values;
    src_values.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        double x = cell.center().x();
        double y = cell.center().y();
        double src = 2 * std::pow(prism::PI, 2);
        src *= std::sin(prism::PI * x);
        src *= std::cos(prism::PI * y);
        src_values[cell.id()] = src;
    }

    auto c = field::Tensor("c", mesh, Matrix3d::Identity());

    using laplacian = diffusion::Corrected<field::Tensor,
                                           scheme::diffusion::nonortho::OverRelaxedCorrector,
                                           field::Scalar>;
    auto source = field::Scalar("S", mesh, std::move(src_values));

    auto eqn = eqn::Transport<field::Scalar>(
        laplacian(c, P),                                              // -∇.∇p
        source::ConstantScalar<Sign::Positive, field::Scalar>(source) // = S
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    auto nOrthogonalCorrectors = 5;

    for (int i = 0; i < nOrthogonalCorrectors; ++i) {
        solver.solve(eqn, 15, 1e-20);
    }

    return P;
}

TEST_CASE("test poisson equation unstructured", "[poisson]") {
    auto P = testPoissonWithMesh("tests/cases/poisson/mesh.unv");
    double diff_norm = l2NormRel(P.values(), poisson_analytical_solution(P.mesh()).values());
    REQUIRE(diff_norm < 0.078); // for poisson/mesh.unv it should be = 0.0769
}

TEST_CASE("test poisson equation structured", "[poisson]") {
    auto P = testPoissonWithMesh("tests/cases/poisson/mesh_hex.unv");
    double diff_norm = l2NormRel(P.values(), poisson_analytical_solution(P.mesh()).values());
    REQUIRE(diff_norm < 0.0004); // for poisson/mesh_hex.unv it should be = 0.000323
}
