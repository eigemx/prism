#include <prism/prism.h>

#include <filesystem>

using namespace prism;
using namespace prism::scheme;
using namespace prism::field;
namespace fs = std::filesystem;

auto main(int argc, char* argv[]) -> int {
    log::setLevel(log::Level::Info);

    log::info("laplacianSolver - A steady state heat equation solver");

    if (argc < 2) {
        log::error("Usage: laplacianSolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = argv[1]; // NOLINT

    // read mesh
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto boundary_file = fs::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = Scalar("T", mesh, 300.0);

    // diffusion coefficient
    // Note: this does not have to be a tensor, but it is just a demonstration of how to use
    // different diffusion coefficient fields.
    auto kappa = Tensor("kappa", mesh, Matrix3d::Identity() * 1e-5);

    auto eqn = eqn::Transport(
        diffusion::Corrected<Tensor, diffusion::nonortho::OverRelaxedCorrector, Scalar>(kappa,
                                                                                        T));

    // solve
    auto solver = solver::BiCGSTAB<Scalar>();
    auto nIter = 40;
    auto nNonOrthoIter = 2;

    for (int iter = 0; iter < nIter; iter++) {
        for (int j = 0; j < nNonOrthoIter; j++) {
            solver.solve(eqn, 10, 1e-20);
        }
    }

    prism::exportToVTU(eqn.field(), "solution.vtu");
    auto gradT_x = ops::grad(T, Coord::X);
    prism::exportToVTU(gradT_x, "gradT_x_ls.vtu");

    return 0;
}
