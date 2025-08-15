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
    auto U = field::Velocity("U", mesh, {.0, .0, .0});
    auto p = field::Pressure("P", mesh, 0.0);
    auto nu = field::UniformScalar("nu", mesh, 1e-3);

    using div = LinearUpwind<field::Velocity, field::VelocityComponent>;
    using laplacian = diffusion::Corrected<field::UniformScalar,
                                           diffusion::nonortho::OverRelaxedCorrector,
                                           field::VelocityComponent>;
    using grad = source::Gradient<Sign::Negative, field::Pressure>;

    auto nOuterIter = 50;
    auto mdot = U.clone();

    algo::SIMPLEParameters params;

    auto uEqn = eqn::Momentum(div(mdot, U.x()),     // ∇.(Uu)
                              laplacian(nu, U.x()), // -∇.(ν∇u)
                              grad(p, Coord::X)     // = -∂p/∂x
    );

    auto vEqn = eqn::Momentum(div(mdot, U.y()),     // ∇.(Uv)
                              laplacian(nu, U.y()), // -∇.(ν∇v)
                              grad(p, Coord::Y)     // = -∂p/∂y
    );
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

    uEqn.setUnderRelaxFactor(params.momentum_urf);
    vEqn.setUnderRelaxFactor(params.momentum_urf);

    std::vector<eqn::Momentum*> momentum_eqns {&uEqn, &vEqn};

    for (auto outer_iteration = 0; outer_iteration < nOuterIter; ++outer_iteration) {
        log::info("Outer iteration: {}", outer_iteration);
        auto SIMPLE = algo::IncompressibleSIMPLE(params);
        SIMPLE.step(std::span<eqn::Momentum*>(momentum_eqns), U, mdot, p);
    }

    exportToVTU(U.x(), "solution_x.vtu");
    exportToVTU(U.y(), "solution_y.vtu");
    exportToVTU(p, "pressure.vtu");
}
