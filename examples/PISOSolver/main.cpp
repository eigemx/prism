#include <prism/algorithm/PISO.h>
#include <prism/prism.h>
#include <prism/scheme/temporal/adam_moulton.h>

#include <filesystem>

using namespace prism;

auto main(int argc, char* argv[]) -> int {
    using namespace prism::scheme;
    using namespace prism::scheme::convection;
    log::setLevel(log::Level::Info);

    if (argc < 2) {
        log::error("usage: {} [mesh-file]", argv[0]);
        return 1;
    }
    const auto* unv_file_name = argv[1];

    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    auto U = makeShared<field::Velocity>("U", mesh, Vector3d {.0, .0, .0});
    auto p = makeShared<field::Pressure>("P", mesh, 0.0);
    auto nu = makeShared<field::Scalar>("nu", mesh, 0.01);
    auto mdot = makeShared<field::Velocity>(U->clone());

    U->x()->setHistorySize(2);
    U->y()->setHistorySize(2);
    U->x()->updatePrevTimeSteps();
    U->y()->updatePrevTimeSteps();

    f64 dt = 0.005;
    f64 end_time = 1.0;
    auto n_time_steps = size_t(end_time / dt);

    using div = LinearUpwind;
    using lap = diffusion::Corrected;
    using grad = source::Gradient<Sign::Negative>;
    using ddt = temporal::AdamMoulton;

    auto uEqn =
        eqn::Momentum(ddt(U->x(), dt), div(mdot, U->x()), lap(nu, U->x()), grad(p, VectorCoord::X));
    auto vEqn =
        eqn::Momentum(ddt(U->y(), dt), div(mdot, U->y()), lap(nu, U->y()), grad(p, VectorCoord::Y));

    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Transport>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Transport>>();

    std::vector<eqn::Momentum*> meq {&uEqn, &vEqn};

    algo::PISOControls controls;
    controls.outer_iterations = 1;
    controls.pressure_correction_steps = 3; // SIMPLE predictor + 2 PRIME corrector
    controls.non_ortho_correctors = 1;

    f64 time = 0.0;

    output::TimeWriter writer("result", /*interval=*/2, output::WriteMode::Ascii);
    writer.add(p);
    writer.add(U);

    for (size_t ts = 0; ts < n_time_steps; ++ts) {
        time = (ts + 1) * dt;
        log::info("t={:.5f}  |U|max={:.4e}  |p|max={:.4e}",
                  time,
                  U->x()->values().lpNorm<Eigen::Infinity>(),
                  p->values().lpNorm<Eigen::Infinity>());
        auto reps = algo::PISO(controls).step(std::span<eqn::Momentum*>(meq), U, mdot, p);
        for (auto& r : reps) {
            if (std::isnan(r.final_residual)) {
                log::warn("NaN residual in field {}", r.field_name);
            }
        }
        U->x()->updatePrevTimeSteps();
        U->y()->updatePrevTimeSteps();
        writer.writeAndAdvance(dt);
    }

    return 0;
}
