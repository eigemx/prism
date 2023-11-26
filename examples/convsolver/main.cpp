#include <fmt/core.h>
#include <prism/prism.h>

#include <algorithm>

#include "prism/gradient/gradient.h"
#include "prism/nonortho/nonortho.h"
#include "prism/schemes/convection.h"


auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    fmt::println("convsolver - A steady state temperature advection solver");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        fmt::println("Usage: convsolver [mesh-file]");
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
        fmt::println(
            "Error: No boundary patch with name `inlet` was found, cannot set velocity field.");
        return 1;
    }

    // Set a uniform velocity field, with value equal to inlet velocity;
    Vector3d inlet_velocity = inlet_patch->get_vector_bc("velocity");
    auto U = VectorField("velocity", mesh, inlet_velocity);

    // A zero field, just to demonstrate how to add arbitray constant source terms
    auto useLessField = ScalarField("zero", mesh, 0.0);

    // solve for temperature advection: ∇.(ρUT) - ∇.(κ ∇T) = S
    // where ρ is the density and U is the velocity vector, and S is an arbitraty constant source
    auto eqn = TransportEquation(
        // Add discretization schemes
        convection::SecondOrderUpwind(rho, U, T), // ∇.(ρUT)
        diffusion::Diffusion(1e-2, T),            // - ∇.(κ ∇T)
        source::ConstantScalar(useLessField)      // S (sources are always added to the RHS)
    );

    // solve
    auto solver = solver::BiCGSTAB();
    solver.solve(eqn, 2000, 1e-3);

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    return 0;
}