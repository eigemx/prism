#include <prism/core.h>

#include <memory>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    print_header();
    fmt::print("convsolver - A steady state temperature advection solver\n");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        error("Usage: convsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    fmt::print("Loading mesh file {}...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name).to_pmesh();
    fmt::print("Okay.\n");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);
    auto T_grad = gradient::create<gradient::GreenGauss>(T);

    // density field
    auto rho = ScalarField("density", mesh, 1.18);

    // set up a unifform velocity field defined over the mesh
    // set the velocity of the field to be the same as the inlet value
    const auto& inlet_patch = std::find_if(
        mesh.boundary_patches().begin(), mesh.boundary_patches().end(), [](const auto& patch) {
            return patch.name() == "inlet";
        });

    if (inlet_patch == mesh.boundary_patches().end()) {
        error("No inlet patch found!");
        return 1;
    }

    Vector3d inlet_velocity = inlet_patch->get_vector_bc("velocity");
    auto U = VectorField("velocity", mesh, inlet_velocity);

    // solve for temperature convection: ∇.(ρuT) - ∇.(κ ∇T) = 0
    // where ρ is the density and u is the velocity
    auto eqn = Equation(
        diffusion::Diffusion<diffusion::NonOrthoCorrection::OverRelaxed>(1e-2, T, T_grad),
        convection::Convection<convection::Upwind>(rho, U, T, T_grad));

    // solve
    auto solver = solver::BiCGSTAB();
    solver.solve(eqn, 2000, 1e-3);

    prism::export_field(eqn.scalar_field(), "solution.vtu");

    return 0;
}