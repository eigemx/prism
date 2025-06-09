#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>

auto peclet_number(double rho, double kappa, double u, double L) {
    return (rho * u * L) / kappa;
}

auto advection_1d(double x, double u) -> double {
    double phi_west = 1;
    double phi_east = 0;
    double pecklet = peclet_number(1.18, 1e-2, u, 1.0);

    auto result = (std::exp(pecklet * (x)) - 1) / (std::exp(pecklet) - 1);
    result *= phi_east - phi_west;
    result += phi_west;
    return result;
}

TEST_CASE("solve advection equation at u = 0.05 m/s, Pe ~= 5", "[advection]") {
    using namespace prism;

    log::setLevel(log::Level::Error);

    const auto* mesh_file = "tests/cases/channel1d_coarse/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto T = field::Scalar("temperature", mesh, 0.5);
    auto rho = field::Scalar("rho", mesh, 1.18);

    const auto& inlet_patch = std::find_if(
        mesh.boundaryPatches().begin(), mesh.boundaryPatches().end(), [](const auto& patch) {
            return patch.name() == "inlet";
        });


    Vector3d inlet_velocity = inlet_patch->getVectorBoundaryCondition("U");
    auto U = field::Velocity("U", mesh, inlet_velocity);

    auto kappa = field::UniformScalar("kappa", mesh, 1e-2);
    auto eqn = eqn::Transport(scheme::convection::SecondOrderUpwind(rho, U, T), // ∇.(ρUT)
                              scheme::diffusion::NonCorrected(kappa, T) // - ∇.(κ ∇T) = 0
    );

    auto solver = solver::BiCGSTAB<field::Scalar>();
    solver.solve(eqn, 100, 1e-20, 1);

    VectorXd analytical_solution;
    analytical_solution.resize(T.mesh().cellCount());
    for (const auto& cell : T.mesh().cells()) {
        analytical_solution[cell.id()] = advection_1d(cell.center().x(), inlet_velocity.x());
    }

    VectorXd diff_vec = analytical_solution.array() - T.values().array();


    // TODO: delete this
    auto analytic_sol_field = field::Scalar("analytical_solution", mesh, analytical_solution);
    export_field_vtu(T, "T_1d.vtu");
    export_field_vtu(analytic_sol_field, "analytical_solution.vtu");

    REQUIRE(diff_vec.norm() < 0.1);

    std::vector<double> T_vec;
    for (const auto& cell : T.mesh().cells()) {
        T_vec.push_back(T.valueAtCell(cell.id()));
    }

    REQUIRE(std::is_sorted(T_vec.rbegin(), T_vec.rend()));
}
