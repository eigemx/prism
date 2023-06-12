#include <prism/core.h>

#include <iostream>
#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    print("diffsolver - A steady state temperature diffusion solver\n");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: diffsolver [unv-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    print("Loading mesh file {}...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name).to_pmesh();
    print("Okay.\n");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);

    // solve for temperature diffision: -∇.(κ ∇T) = 0
    // where κ is the diffusion coefficient
    auto diff = diffusion::Diffusion(4e-5, T);

    // assemble the equation
    auto eqn = Equation(T, {&diff});

    // solve
    auto solver = solver::BiCGSTAB();
    solver.solve_bicg(eqn, 1000, 1e-10);

    prism::export_field(eqn.scalar_field(), "solution.vtu");

    return 0;
}
