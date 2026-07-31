#include <prism/prism.h>

#include <filesystem>

#include "prism/algorithm/SIMPLE.h"

using namespace prism;
namespace fs = std::filesystem;

auto main(int argc, char* argv[]) -> int {
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    log::setLevel(log::Level::Info);

    if (argc < 2) {
        log::error("usage: {} [mesh-file]", argv[1]); // NOLINT
        return 1;
    }

    const auto* unv_file_name = argv[1]; // NOLINT

    // read mesh
    log::info("Reading `fields.json` file...");
    auto boundary_file = fs::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set mesh fields
    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {.0, .0, .0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 1e-3);

    using div = LinearUpwind;
    using laplacian = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;

    auto nOuterIter = 50;
    auto mdot = makeShared<field::Velocity>(U->clone());

    algo::SIMPLEControls controls;

    auto uEqn = eqn::Momentum(div(mdot, U->x()),      // ∇.(Uu)
                              laplacian(nu, U->x()),  // -∇.(ν∇u)
                              grad(p, VectorCoord::X) // = -∂p/∂x
    );

    auto vEqn = eqn::Momentum(div(mdot, U->y()),      // ∇.(Uv)
                              laplacian(nu, U->y()),  // -∇.(ν∇v)
                              grad(p, VectorCoord::Y) // = -∂p/∂y
    );

    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

    uEqn.setUnderRelaxFactor(controls.momentum_urf);
    vEqn.setUnderRelaxFactor(controls.momentum_urf);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    for (auto outer_iteration = 0; outer_iteration < nOuterIter; ++outer_iteration) {
        log::info("Outer iteration: {}", outer_iteration);
        auto reports = algo::IncompressibleSIMPLE(controls).step(
            std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
        for (const auto& entry : reports) {
            log::info("Residuals: Initial = {:.4e} | Final: {:.4e} (nIterations = {}) | field: {}",
                      entry.initial_residual,
                      entry.final_residual,
                      entry.n_iterations,
                      entry.field_name);
        }
    }

    auto scalars = std::vector<const field::IScalar*> {p.get(), nu.get()};
    auto vectors = std::vector<const field::IVector*> {U.get()};
    output::toVtu(scalars, vectors, {}, "solution.vtu");
}
