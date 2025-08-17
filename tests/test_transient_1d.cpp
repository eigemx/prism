#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>

#include "prism/field/scalar.h"
#include "prism/scheme/temporal/adam_moulton.h"
#include "prism/types.h"

auto transientAnalyticalSolution(prism::f64 x, prism::f64 t, prism::f64 kappa) -> prism::f64 {
    auto n_modes = 100;
    prism::f64 sum = 0.0;
    for (auto n = 0; n < n_modes; n++) {
        auto lambda_n = ((2 * n) - 1) * prism::PI / 2;
        prism::f64 mode_val = std::pow(-1, n + 1) * std::exp(-kappa * t * std::pow(lambda_n, 2));
        mode_val *= std::cos(lambda_n * x);
        mode_val /= ((2 * n) - 1);
        sum += mode_val;
    }
    sum *= 4.0 / prism::PI;
    return sum * 100;
}


TEST_CASE("solve transient diffusion equation 1D", "[transiet]") {
    using namespace prism;
    log::setLevel(log::Level::Error);

    // read mesh
    const auto* mesh_file = "tests/cases/versteeg_trans_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto T = field::Scalar("T", mesh, 200.0);
    T.setHistorySize(2); // enable history with a single time step in the past
    T.update(T.values());

    // diffusion coefficient
    auto kappa = prism::field::UniformScalar("kappa", mesh, 1e-3);

    auto dt = 2;

    // solve
    auto solver = solver::BiCGSTAB<prism::field::Scalar>();
    auto nNonOrthoIter = 2;
    auto nTimesteps = 200;

    using prism::scheme::diffusion::nonortho::OverRelaxedCorrector;
    auto eqn =
        eqn::Transport(prism::scheme::temporal::AdamMoulton<prism::field::Scalar>(T, dt), // dT/dt
                       prism::scheme::diffusion::Corrected<prism::field::UniformScalar,
                                                           OverRelaxedCorrector,
                                                           prism::field::Scalar>(kappa, T));

    
    std::vector<double> diff_norm_vec;

    for (auto timestep = 0; timestep < nTimesteps; timestep++) {
        log::info("Solving timestep {}/{} at time = {}", timestep + 1, nTimesteps, dt * timestep);

        // update the field time history
        T.update(T.values());

        for (auto i = 0; i < nNonOrthoIter; i++) {
            solver.solve(eqn, 10, 1e-20);
        }

        VectorXd analytical_sol(mesh->cellCount());
        analytical_sol.setZero();

        std::for_each(mesh->cells().begin(), mesh->cells().end(), [&](const auto& cell) {
            auto x = cell.center().x();
            auto t = dt * timestep;
            analytical_sol[cell.id()] =
                transientAnalyticalSolution(x, t, kappa.valueAtCell(cell.id()));
        });

        auto t = dt * timestep;

        if (t == 20 || t == 40 || t == 60) {
            auto diff_norm = (T.values() - analytical_sol).norm() / analytical_sol.norm();
            diff_norm_vec.push_back(diff_norm);
        }
    }
    for (auto diff_norm : diff_norm_vec) {
        REQUIRE(diff_norm < 0.15);
    }
}
