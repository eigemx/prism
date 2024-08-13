#include <fmt/core.h>
#include <prism/prism.h>
#include <spdlog/spdlog.h>

#include <filesystem>


auto main(int argc, char* argv[]) -> int {
    using namespace prism;

    spdlog::set_level(spdlog::level::level_enum::debug);
    if (argc < 2) {
        fmt::println("usage: {} [mesh-file]", argv[1]); // NOLINT
        return 1;
    }

    const auto* unv_file_name = argv[1]; // NOLINT

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "boundary.txt";
    fmt::println("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).to_pmesh();

    auto mu = field::UniformScalar("viscosity", mesh, 1e-6);
    auto rho = field::Scalar("density", mesh, 1.18);
    auto U = field::Velocity("velocity", mesh, {0.05, 0.05, 0.0});
    auto P = field::Pressure("pressure", mesh, 1.0);

    auto uEqn = eqn::Momentum(
        scheme::convection::Upwind<field::VelocityComponent>(rho, U, U.x()), // ∇.(ρUu)
        scheme::diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>(
            mu, U.x()) //, // - ∇.(μ∇u)
        // scheme::source::Gradient<scheme::source::SourceSign::Negative, field::Pressure>(
        // P, Coord::X) // ∂p/∂x
    );

    auto vEqn = eqn::Momentum(
        scheme::convection::Upwind<field::VelocityComponent>(rho, U, U.y()),
        scheme::diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>(mu,
                                                                                        U.y()) //,
        // scheme::source::Gradient<scheme::source::SourceSign::Negative, field::Pressure>(
        // P, Coord::Y)
    );

    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();


    auto solver = solver::BiCGSTAB<field::VelocityComponent,
                                   solver::ImplicitUnderRelaxation<field::VelocityComponent>>();

    auto p_solver =
        solver::BiCGSTAB<field::Pressure, solver::ImplicitUnderRelaxation<field::Pressure>>();

    for (auto i = 0; i < 2; ++i) {
        uEqn.updateCoeffs();
        vEqn.updateCoeffs();

        // calculate coefficients for the pressure equation
        const auto vol_vec = mesh::cells_volume_vec(mesh);
        const auto& uEqn_diag = uEqn.matrix().diagonal();
        const auto& vEqn_diag = vEqn.matrix().diagonal();

        auto D_data = std::vector<prism::Matrix3d>();
        D_data.reserve(mesh.nCells());

        auto Du = vol_vec.array() / (uEqn_diag.array() + prism::EPSILON);
        auto Dv = vol_vec.array() / (vEqn_diag.array() + prism::EPSILON);

        for (std::size_t i = 0; i < mesh.nCells(); ++i) {
            // clang-format off
            Matrix3d Di;
            Di << Du[i], 0,     0,
                  0,     Dv[i], 0,
                  0,     0,     0;
            // clang-format on
            D_data.emplace_back(std::move(Di));
        }

        auto D = prism::field::Tensor("D", mesh, D_data);

        spdlog::info("Solving y-eqn::Momentum equation");
        solver.solve(vEqn, 2, 1e-3, 0.9);

        spdlog::info("Solving x-eqn::Momentum equation");
        solver.solve(uEqn, 2, 1e-3, 0.9);

        // Rhie-Chow interpolation for velocity face values
        ops::correctRhieChow(U, D, P);

        // pressure equation
        // Few things missing to implement:
        // 1) it's actually ρD not just D
        // 2) ∇.(ρU) not ∇.U
        auto P_prime = field::Pressure("pressure", mesh, 0.0);
        auto pEqn = eqn::Transport<field::Pressure>(
            scheme::diffusion::NonCorrected<field::Tensor, field::Pressure>(D,
                                                                            P_prime) //,
            // scheme::source::Divergence<scheme::source::SourceSign::Negative,
            // field::Velocity>(U)
        );

        pEqn.updateCoeffs();
        spdlog::info("Solving pressure correction equation");
        p_solver.solve(pEqn, 10, 1e-5, 1);

        // update velocity fields
        auto p_grad = gradient::LeastSquares(P_prime);
        for (const auto& cell : mesh.cells()) {
            auto correction = -D.valueAtCell(cell) * p_grad.gradAtCell(cell);
            U.x()[cell.id()] += correction[0];
            U.y()[cell.id()] += correction[1];
        }

        P.values() = P.values().array() + (0.85 * P_prime.values().array());
    }

    prism::export_field_vtu(uEqn.field(), "solution.vtu");
}