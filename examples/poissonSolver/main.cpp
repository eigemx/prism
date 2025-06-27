#include <fmt/core.h>
#include <prism/prism.h>

#include <cmath>
#include <string>
#include <vector>

#include "prism/scheme/nonortho.h"

using namespace prism;
using namespace prism::scheme;

auto solution(const auto& mesh) -> prism::field::Scalar {
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

auto l2NormRelative(const Vector3d& x, const Vector3d& x_ref) -> double {
    return (x - x_ref).norm() / x_ref.norm();
}

auto main(int argc, char* argv[]) -> int {
    log::setLevel(log::Level::Debug);

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        log::error("Usage: poisson [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    auto P = field::Scalar("P", mesh, 0.0);

    // create source term
    /// TODO: use a more generic way to create source terms using mapped scalar fields
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
    auto nNonOrthogonalCorrectors = 5;

    for (int i = 0; i < nNonOrthogonalCorrectors; ++i) {
        solver.solve(eqn, 15, 1e-20);
    }

    prism::export_field_vtu(P, "solution.vtu");
    prism::export_field_vtu(solution(mesh), "analytical.vtu");

    VectorXd diff = P.values() - solution(mesh).values();
    auto diff_field = field::Scalar("diff", mesh, diff);
    prism::export_field_vtu(diff_field, "diff.vtu");

    fmt::print("relative l2-norm: {}\n", l2NormRelative(P.values(), solution(mesh).values()));

    auto grad_p_anal_x = field::Scalar("grad_p_anal", mesh, 0.0);
    auto grad_p_anal_y = field::Scalar("grad_p_anal", mesh, 0.0);

    for (const auto& cell : mesh->cells()) {
        auto x = cell.center().x();
        auto y = cell.center().y();
        auto PI = prism::PI;
        grad_p_anal_x.values()[cell.id()] = PI * std::cos(x * PI) * std::cos(y * PI);
        grad_p_anal_y.values()[cell.id()] = -PI * std::sin(x * PI) * std::sin(y * PI);
    }

    prism::export_field_vtu(grad_p_anal_x, "grad_p_anal_x.vtu");
    prism::export_field_vtu(grad_p_anal_y, "grad_p_anal_y.vtu");
    prism::export_field_vtu(ops::grad(P, Coord::X), "grad_p_x.vtu");
    prism::export_field_vtu(ops::grad(P, Coord::Y), "grad_p_y.vtu");

    // print relative l2-norm of the analytical solution for each pressure gradient component
    fmt::print("relative l2-norm of the analytical solution for grad_p_x: {}\n",
               l2NormRelative(grad_p_anal_x.values(), ops::grad(P, Coord::X).values()));
    fmt::print("relative l2-norm of the analytical solution for grad_p_y: {}\n",
               l2NormRelative(grad_p_anal_y.values(), ops::grad(P, Coord::Y).values()));

    return 0;
}
