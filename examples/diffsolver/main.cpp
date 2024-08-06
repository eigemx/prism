#include <prism/prism.h>

#include <string>
#include <vector>

#include "fmt/core.h"
#include "prism/field/field.h"
#include "prism/gradient/gradient.h"
#include "prism/nonortho/nonortho.h"
#include "prism/schemes/diffusion.h"
#include "prism/schemes/source.h"
#include "spdlog/common.h"
#include "spdlog/spdlog.h"


auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    spdlog::set_level(spdlog::level::level_enum::debug);

    fmt::println("diffsolver - A steady state temperature diffusion solver");
    fmt::println("");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        fmt::println("Usage: diffsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "boundary.txt";
    fmt::print("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).to_pmesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = field::Scalar("temperature", mesh, 300.0);

    // define a source term
    prism::VectorXd source_field_data = VectorXd::Zero(mesh.nCells());
    for (const auto& cell : mesh.cells()) {
        const auto& center = cell.center();
        if (center.norm() <= 0.15) {
            source_field_data[cell.id()] = 100000.0;
        }
    }
    auto S = field::Scalar("S", mesh, source_field_data);

    // assemble the equation
    // solve for temperature diffision: -∇.(κ ∇T) = 0
    // where κ is the diffusion coefficient
    auto kappa = field::UniformScalar("kappa", mesh, 1e-5);
    auto eqn = TransportEquation(scheme::diffusion::CorrectedDiffusion(kappa, T),
                                 scheme::source::ConstantScalar(S));

    // solve
    auto solver =
        solver::BiCGSTAB<field::Scalar, solver::ImplicitUnderRelaxation<field::Scalar>>();
    solver.solve(eqn, 100, 1e-5, 1);

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    auto gradT_x = gradient::LeastSquares(T).gradient_field().x();
    prism::export_field_vtu(gradT_x, "gradT_x_ls.vtu");

    return 0;
}
