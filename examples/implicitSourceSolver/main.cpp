#include <prism/prism.h>

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "prism/field/scalar.h"
#include "prism/scheme/source.h"

using namespace prism;
using namespace prism::scheme;

auto solution(const auto& mesh) -> prism::field::Scalar {
    // y = sinh(2x + 1)
    prism::VectorXd sol;
    sol.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        double x = cell.center().x();
        double sol_cell = std::sinh(2 * x + 1);
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

    auto y = field::Scalar("y", mesh, 0.0);

    auto c = field::UniformScalar("c", mesh, -1.0);

    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport<field::Scalar>(
        laplacian(c, y) // ∇.∇y
        // source::ImplicitField<source::SourceSign::Positive, field::Scalar>(4, y) // = 4y
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    auto nOuterIter = 5;

    for (int i = 0; i < nOuterIter; ++i) {
        solver.solve(eqn, 15, 1e-20);
    }

    eqn.updateCoeffs();
    std::ofstream matrix("A.csv");
    matrix << eqn.matrix();
    matrix.close();

    std::ofstream rhs("b.csv");
    rhs << eqn.rhs();
    rhs.close();

    prism::export_field_vtu(y, "solution.vtu");
    prism::export_field_vtu(solution(mesh), "analytical.vtu");

    VectorXd diff = y.values() - solution(mesh).values();
    auto diff_field = field::Scalar("diff", mesh, diff);
    prism::export_field_vtu(diff_field, "diff.vtu");

    fmt::print("relative l2-norm: {}\n", l2NormRelative(y.values(), solution(mesh).values()));

    return 0;
}
