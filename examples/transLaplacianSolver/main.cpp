#include <prism/prism.h>

#include <filesystem>

#include "prism/scheme/temporal/adam_moulton.h"

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
    T.setHistorySize(1);     // enable history with a single time step in the past
    T.updatePrevTimeSteps(); // sets the initial value of the field at t = 0

    // diffusion coefficient
    auto kappa = UniformScalar("kappa", mesh, 1e-3);

    auto dt = 2;

    // solve
    auto solver = solver::BiCGSTAB<Scalar>();
    auto nNonOrthoIter = 2;
    auto nTimesteps = 200;

    using diffusion::nonortho::OverRelaxedCorrector;
    auto eqn = eqn::Transport(
        temporal::BackwardEuler<Scalar>(T, dt), // dT/dt
        diffusion::Corrected<UniformScalar, OverRelaxedCorrector, Scalar>(kappa, T));

    for (auto timestep = 0; timestep < nTimesteps; timestep++) {
        log::info("Solving timestep {}/{} at time = {}", timestep + 1, nTimesteps, dt * timestep);

        for (auto i = 0; i < nNonOrthoIter; i++) {
            solver.solve(eqn, 10, 1e-20);
        }

        // update the field time history
        T.updatePrevTimeSteps();

        log::info("T = {}", T.values().mean());
        log::info("T_prev = {}", T.prevValues().value().mean());
        // log::info("T_prev_prev = {}", T.prevPrevValues().value().mean());
    }
    prism::exportToVTU(eqn.field(), "solution.vtu");

    return 0;
}
