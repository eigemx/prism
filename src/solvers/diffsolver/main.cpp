#include <prism/core.h>

#include <iostream>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    fmt::print("diffsolver - A steady state temperature diffusion solver\n");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: diffsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    fmt::print("Loading mesh file {}...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name).to_pmesh();
    fmt::print("Okay.\n");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);


    // solve for temperature diffision: -∇.(κ ∇T) = 0
    // where κ is the diffusion coefficient
    auto diff = diffusion::Diffusion<diffusion::NonOrthoCorrection::None>(1, T);
    auto S = ScalarField("S", mesh).map([](const mesh::Cell& cell) {
        const auto& center = cell.center();
        if (center.norm() <= 0.15) {
            return 1000.0;
        }
        return 0.0;
    });
    auto source = source::ConstantScalar(S);
    auto sink = source::ImplicitPhi<source::SourceSign::Negative>(T);


    // assemble the equation
    auto eqn = Equation(T, {&diff, &source, &sink});

    // solve
    auto solver = solver::BiCGSTAB();

    solver.solve(eqn, 1000, 1e-3);

    prism::export_field(eqn.scalar_field(), "solution.vtu");

    return 0;
}
