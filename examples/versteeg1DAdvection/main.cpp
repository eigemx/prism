#include <prism/prism.h>

#include <algorithm>
#include <filesystem>

// This example is based on section 5.6 from "An Introduction to Computational Fluid Dynamics -The
// Finite Volume Method" - Second Edition by H K Versteeg and W Malalasekera
// it's also implemented as a test case.

auto analytical_solution(double u, const prism::SharedPtr<prism::mesh::PMesh>& mesh)
    -> prism::field::Scalar {
    prism::VectorXd sol;
    sol.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        double x = cell.center().x();
        double y = cell.center().y();
        double sol_cell = -(std::exp(u * x / 0.1) - 1) / (std::exp(u / 0.1) - 1);
        sol[cell.id()] = sol_cell + 1;
    }

    return {"analytical_solution", mesh, std::move(sol)};
}

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    log::setLevel(log::Level::Info);

    fmt::println("versteeg1DAdvection - A steady state temperature advection solver");

    // silence clang-tidy pointer arithmetic warnings
    std::vector<std::string> args(argv, argv + argc);

    if (argc < 2) {
        fmt::println("Usage: versteeg1DAdvection [mesh-file]");
        return 1;
    }

    auto unv_file_name = args[1];

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    auto T = field::Scalar("T", mesh, 0.0);
    auto rho = field::UniformScalar("rho", mesh, 1.0);

    // set up a uniform velocity field defined over the mesh
    // set the velocity of the field to be the same as the inlet value
    const auto& inlet_patch = std::find_if(
        mesh->boundaryPatches().begin(), mesh->boundaryPatches().end(), [](const auto& patch) {
            return patch.name() == "Inlet";
        });

    if (inlet_patch == mesh->boundaryPatches().end()) {
        fmt::println(
            "Error: No boundary patch with name `Inlet` was found, cannot set velocity field.");
        return 1;
    }

    // Set a uniform velocity field, with value equal to inlet velocity;
    Vector3d inlet_velocity = inlet_patch->getVectorBoundaryCondition("U");

    log::info("Setting velocity field to {} [m/s]", inlet_velocity.norm());
    auto U = field::Velocity("U", mesh, inlet_velocity);
    field::Velocity rhoU = rho * U;
    auto kappa = field::UniformScalar("kappa", mesh, 0.1);

    log::info("Peclet number = {}", inlet_velocity.x() * 0.2 / 0.1);

    // solve for temperature advection: ∇.(ρUT) - ∇.(κ ∇T) = 0
    using div = scheme::convection::Upwind<field::Velocity, field::Scalar>;
    using laplacian = scheme::diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport(div(rhoU, T),       // ∇.(ρUT)
                              laplacian(kappa, T) // - ∇.(κ ∇T)
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();

    solver.solve(eqn, 5, 1e-20);
    VectorXd diff = eqn.field().values().array() -
                    analytical_solution(inlet_velocity.x(), mesh).values().array();
    auto diff_norm = diff.norm();
    log::info("diff norm = {}", diff_norm);

    prism::exportToVTU(eqn.field(), "solution.vtu");

    return 0;
}
