#include <fmt/core.h>
#include <prism/prism.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <filesystem>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    spdlog::set_level(spdlog::level::level_enum::debug);


    fmt::println("convsolver - A steady state temperature advection solver");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        fmt::println("Usage: convsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "boundary.txt";
    spdlog::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).to_pmesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = field::Scalar("temperature", mesh, 300.0);

    // density field
    auto rho = field::Scalar("rho", mesh, 1.18);

    // set up a unifform velocity field defined over the mesh
    // set the velocity of the field to be the same as the inlet value
    const auto& inlet_patch = std::find_if(
        mesh.boundaryPatches().begin(), mesh.boundaryPatches().end(), [](const auto& patch) {
            return patch.name() == "inlet";
        });

    if (inlet_patch == mesh.boundaryPatches().end()) {
        fmt::println(
            "Error: No boundary patch with name `inlet` was found, cannot set velocity field.");
        return 1;
    }

    // Set a uniform velocity field, with value equal to inlet velocity;
    Vector3d inlet_velocity = inlet_patch->getVectorBoundaryCondition("velocity");
    auto U = field::Vector("velocity", mesh, inlet_velocity);

    // A zero field, just to demonstrate how to add arbitray constant source terms
    auto useLessField = field::Scalar("zero", mesh, 0.0);

    auto kappa = field::UniformScalar("kappa", mesh, 1e-2);

    // solve for temperature advection: ∇.(ρUT) - ∇.(κ ∇T) = S
    // where ρ is the density and U is the velocity vector, and S is an arbitraty constant source
    auto eqn = TransportEquation(
        scheme::convection::Upwind(rho, U, T),              // ∇.(ρUT)
        scheme::diffusion::NonCorrectedDiffusion(kappa, T), // - ∇.(κ ∇T)
        scheme::source::ConstantScalar(useLessField) // S (sources are always added to the RHS)
    );

    // solve
    auto solver =
        solver::BiCGSTAB<field::Scalar, solver::ImplicitUnderRelaxation<field::Scalar>>();
    solver.solve(eqn, 100, 1e-5, 1);

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    auto div = ops::div(U);
    prism::export_field_vtu(div, "div.vtu");

    return 0;
}