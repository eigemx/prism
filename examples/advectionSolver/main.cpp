#include <prism/prism.h>

#include <algorithm>
#include <filesystem>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    log::setLevel(log::Level::Info);

    fmt::println("convsolver - A steady state temperature advection solver");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        fmt::println("Usage: convsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = field::Scalar("temperature", mesh, 300.0);

    // density field
    auto rho = field::UniformScalar("rho", mesh, 1.18);

    // set up a uniform velocity field defined over the mesh
    // set the velocity of the field to be the same as the inlet value
    const auto& inlet_patch = std::find_if(
        mesh->boundaryPatches().begin(), mesh->boundaryPatches().end(), [](const auto& patch) {
            return patch.name() == "inlet";
        });

    if (inlet_patch == mesh->boundaryPatches().end()) {
        fmt::println(
            "Error: No boundary patch with name `inlet` was found, cannot set velocity field.");
        return 1;
    }

    // Set a uniform velocity field, with value equal to inlet velocity;
    Vector3d inlet_velocity = inlet_patch->getVectorBoundaryCondition("U");

    log::info("Setting velocity field to {} [m/s]", inlet_velocity.norm());
    auto U = field::Velocity("U", mesh, inlet_velocity);
    field::Velocity rhoU = rho * U;
    auto kappa = field::UniformScalar("kappa", mesh, 1e-2);

    // solve for temperature advection: ∇.(ρUT) - ∇.(κ ∇T) = 0
    // where ρ is the density and U is the velocity vector, and S is an arbitraty constant source
    using div = scheme::convection::LinearUpwind<field::Velocity, field::Scalar>;
    using laplacian =
        scheme::diffusion::Corrected<field::UniformScalar,
                                     scheme::diffusion::nonortho::OverRelaxedCorrector,
                                     field::Scalar>;

    auto eqn = eqn::Transport(div(rhoU, T),       // ∇.(ρUT)
                              laplacian(kappa, T) // - ∇.(κ ∇T)
    );

    eqn.setUnderRelaxFactor(0.95);

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    auto nOuterIterations = 10;

    for (int i = 0; i < nOuterIterations; ++i) {
        solver.solve(eqn, 10, 1e-20);
    }

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    auto div_U = ops::div(U);
    prism::export_field_vtu(div_U, "div.vtu");

    export_field_vtu(U.x().clone(), "U_x.vtu");
    export_field_vtu(U.clone().y(), "U_y.vtu");
    return 0;
}
