#include <prism/prism.h>

#include <filesystem>

#include "prism/scheme/temporal/backward_euler.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::field;
namespace fs = std::filesystem;

auto main(int argc, char* argv[]) -> int {
    log::setLevel(log::Level::Info);

    log::info("transLaplacianSolver - A transient diffusion equation solver");

    if (argc < 2) {
        log::error("Usage: diffsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = argv[1]; // NOLINT

    // read mesh
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto boundary_file = fs::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = Scalar("T", mesh, 0.0);

    // diffusion coefficient
    auto kappa = UniformScalar("kappa", mesh, 1e-3);

    auto dt = 1e-2;

    using diffusion::nonortho::OverRelaxedCorrector;
    auto eqn = eqn::Transport(temporal::BackwardEuler<Scalar>(T, dt), // dT/dt
                              diffusion::NonCorrected<UniformScalar, Scalar>(kappa, T));

    // solve
    auto solver = solver::BiCGSTAB<Scalar>();
    auto nNonOrthoIter = 1;
    auto nTimesteps = 9000;

    for (auto timestep = 0; timestep < nTimesteps; timestep++) {
        log::info("Solving timestep {}/{} at time = {}", timestep + 1, nTimesteps, dt * timestep);
        solver.solve(eqn, 10, 1e-20);
    }

    prism::exportToVTU(eqn.field(), "solution.vtu");

    return 0;
}
