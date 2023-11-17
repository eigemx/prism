#include <fmt/core.h>
#include <prism/prism.h>

#include "prism/schemes/source.h"

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
    fmt::println("Okay.");

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = ScalarField("temperature", mesh, 300.0);

    // density field
    auto rho = ScalarField("rho", mesh, 1.18);

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

    // solve for temperature convection: ∇.(ρUT) - ∇.(κ ∇T) = 0
    // where ρ is the density and u is the velocity vector
    auto eqn = TransportEquation(
        diffusion::Diffusion(1e-2, T),
        convection::SecondOrderUpwind<>(rho, U, T),
        source::Divergence(U)); // This should not affect the solution, because ∇.U = 0

    // solve
    auto solver = solver::BiCGSTAB();
    solver.solve(eqn, 2000, 1e-3);

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    return 0;
}