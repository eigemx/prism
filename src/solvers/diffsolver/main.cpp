#include <prism/core.h>

#include <iostream>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    fmt::println("diffsolver - A steady state temperature diffusion solver\n");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: diffsolver [mesh-file]");
        return 1;
    }

    // read mesh
    fmt::print("Loading mesh file {}...", args[1]);
    auto mesh = mesh::UnvToPMeshConverter(args[1]).to_pmesh();
    fmt::println("Okay.");

    fmt::print("Reordering mesh cells...");
    auto cm = mesh::CuthillMckee(mesh);
    cm.reorder();
    fmt::println("Okay.");
    fmt::println("");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);
    auto T_grad = gradient::create<gradient::LeastSquares>(T);

    // solve for temperature diffision: -∇.(κ ∇T) = 0
    // where κ is the diffusion coefficient
    auto diff = diffusion::Diffusion<diffusion::NonOrthoCorrection::OverRelaxed>(1, T, T_grad);

    // define a source term
    auto S = ScalarField("S", mesh).map([](const mesh::Cell& cell) {
        const auto& center = cell.center();
        if (center.norm() <= 0.15) {
            return 100000.0;
        }
        return 0.0;
    });
    auto source = source::ConstantScalar(S);

    // assemble the equation
    auto eqn =
        Equation(diffusion::Diffusion<diffusion::NonOrthoCorrection::OverRelaxed>(1, T, T_grad));
    eqn.add_scheme(source::ConstantScalar(S));

    // solve
    auto solver = solver::BiCGSTAB();
    solver.solve(eqn, 100, 1e-3);

    prism::export_field(eqn.scalar_field(), "solution.vtu");

    return 0;
}
