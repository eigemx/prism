#include <prism/core.h>

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
    auto mesh = mesh::UnvToPMesh(unv_file_name).to_pmesh();
    print("Okay.\n");

    // Set up the temperature field
    auto T = ScalarField("temperature", mesh, 300.0);

    // Solve for temperature diffision
    auto diff = diffusion::Linear(4e-5, T);

    // Construct the equation
    auto eqn = Equation(T, {&diff});

    // solve
    auto solver = solver::GaussSeidel();
    solver.solve(eqn, 5000, 1e-10);

    prism::export_field(eqn.scalar_field(), "sol.vtu");
}
