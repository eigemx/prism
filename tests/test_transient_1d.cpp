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


TEST_CASE("solve transient diffusion equation 1D", "[transient]") {
    using namespace prism;
    log::setLevel(log::Level::Error);

    // read mesh
    const auto* mesh_file = "tests/cases/versteeg_trans_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto T = std::make_shared<field::Scalar>("T", mesh, 200.0);
    T->setHistorySize(2);
    T->update(T->values());

    auto kappa = std::make_shared<field::Scalar>("kappa", mesh, 1e-3);

    auto dt = 2;

    auto solver = solver::BiCGSTAB<field::Scalar>();
    auto nNonOrthoIter = 2;
    auto nTimesteps = 200;

    using laplacian = scheme::diffusion::Corrected<field::Scalar>;
    auto eqn = eqn::Transport(scheme::temporal::AdamMoulton(T, dt), laplacian(kappa, T));


    std::vector<double> diff_norm_vec;
    double final_diff_norm = 0.0;

    for (auto timestep = 0; timestep < nTimesteps; timestep++) {
        log::info("Solving timestep {}/{} at time = {}", timestep + 1, nTimesteps, dt * timestep);

        T->update(T->values());

        for (auto i = 0; i < nNonOrthoIter; i++) {
            solver.solve(eqn, 10, 1e-20);
        }

        VectorXd analytical_sol(mesh->cellCount());
        analytical_sol.setZero();

        std::for_each(mesh->cells().begin(), mesh->cells().end(), [&](const auto& cell) {
            auto x = cell.center().x();
            auto t = dt * timestep;
            analytical_sol[cell.id()] =
                transientAnalyticalSolution(x, t, kappa->valueAtCell(cell.id()));
        });

        auto t = dt * timestep;
        auto diff_norm = (T->values() - analytical_sol).norm() / analytical_sol.norm();

        if (t == 40 || t == 60) {
            diff_norm_vec.push_back(diff_norm);
        }

        if (timestep == nTimesteps - 1) {
            final_diff_norm = diff_norm;
        }
    }

    for (std::size_t i = 1; i < diff_norm_vec.size(); i++) {
        REQUIRE(diff_norm_vec[i] < diff_norm_vec[i - 1]);
    }

    REQUIRE(diff_norm_vec[0] < 0.10);
    REQUIRE(diff_norm_vec[1] < 0.10);

    REQUIRE(final_diff_norm < 0.02);
}
