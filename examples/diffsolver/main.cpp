#include <prism/prism.h>

#include <string>
#include <vector>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    fmt::println("diffsolver - A steady state temperature diffusion solver");
    fmt::println("");

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

    //fmt::print("Reordering mesh cells...");
    //auto cm = mesh::CuthillMckee(mesh);
    //cm.reorder();
    //fmt::println("Okay.");
    //fmt::println("");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);
    auto T_grad = gradient::create<gradient::LeastSquares>(T);


    // define a source term
    auto S = ScalarField("S", mesh).map([](const mesh::Cell& cell) {
        const auto& center = cell.center();
        if (center.norm() <= 0.15) {
            return 100000.0;
        }
        return 0.0;
    });

    // assemble the equation
    // solve for temperature diffision: -∇.(κ ∇T) = 0
    // where κ is the diffusion coefficient
    auto eqn = TransportEquation(
        diffusion::Diffusion<diffusion::NonOrthoCorrection::None>(1e-5, T, T_grad),
        source::ConstantScalar(S));

    // solve
    auto solver = solver::BiCGSTAB<solver::ImplicitUnderRelaxation>();
    solver.solve(eqn, 100, 1e-6, 1);

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    return 0;
}
