#include <prism/core.h>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    print("convsolver - A steady state temperature advection solver\n");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: convsolver [unv-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    print("Loading mesh file {}...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name).to_pmesh();
    print("Okay.\n");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);
    auto U = VectorField("velocity", mesh, Vector3d(0.00001, 0.0, 0.0));

    // solve for temperature convection: ∇.(ρuT) - ∇.(κ ∇T) = 0
    // where ρ is the density and u is the velocity
    auto diff = diffusion::Linear(0.958, T, gradient::GreenGauss(T));
    auto conv = convection::CentralDifference(1, U, T);

    // assemble the equation
    auto eqn = Equation(T, {&diff, &conv});

    // solve
    auto solver = solver::GaussSeidel();
    solver.solve(eqn, 1000, 1e-10);

    prism::export_field(eqn.scalar_field(), "solution.vtu");

    const auto& A = eqn.coeff_matrix();
    const auto& b = eqn.rhs_vector();

    // export the coefficient matrix and the RHS vector to 'matrix.csv' and 'rhs.csv' respectively
    std::ofstream matrix_file("matrix.csv");
    matrix_file << A;
    matrix_file.close();

    std::ofstream rhs_file("rhs.csv");
    rhs_file << b;
    rhs_file.close();


    return 0;
}