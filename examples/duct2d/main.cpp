#include "fmt/core.h"
#include "prism/constants.h"
#include "prism/equation.h"
#include "prism/field.h"
#include "prism/field/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/unv.h"
#include "prism/mesh/utilities.h"
#include "prism/nonortho/nonortho.h"
#include "prism/numerics/relax.h"
#include "prism/numerics/solver.h"
#include "prism/operations/rhie_chow.h"
#include "prism/schemes/convection.h"
#include "prism/schemes/diffusion.h"
#include "prism/schemes/source.h"
#include "prism/types.h"
#include "spdlog/common.h"
#include "spdlog/spdlog.h"


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
    fmt::print("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).to_pmesh();

    auto rho = field::Scalar("density", mesh, 1.18);
    auto U = field::Vector("velocity", mesh, {0.05, 0.05, 0.0});
    auto P = field::Scalar("pressure", mesh, 1.0);

    auto uEqn = TransportEquation(
        scheme::convection::Upwind(rho, U, U.x()),             // ∇.(ρUu)
        scheme::diffusion::NonCorrectedDiffusion(1e-6, U.x()), // - ∇.(μ∇u)
        scheme::source::Gradient<scheme::source::SourceSign::Negative>(P, Coord::X), // ∂p/∂x
        scheme::source::Laplacian(1e-6, U.x()) // - ∇.(μ∇u^T)
    );

    auto vEqn = TransportEquation(
        scheme::convection::Upwind(rho, U, U.y()),
        scheme::diffusion::NonCorrectedDiffusion(1e-6, U.y()),
        scheme::source::Gradient<scheme::source::SourceSign::Negative>(P, Coord::Y),
        scheme::source::Laplacian(1e-6, U.y()));

    auto solver = solver::BiCGSTAB();

    for (auto i = 0; i < 2; ++i) {
        uEqn.update_coeffs();
        vEqn.update_coeffs();

        // calculate coefficients for the pressure equation
        const auto vol_vec = mesh::cells_volume_vec(mesh);
        const auto& uEqn_diag = uEqn.matrix().diagonal();
        const auto& vEqn_diag = vEqn.matrix().diagonal();

        auto D_data = std::vector<prism::Matrix3d>();
        D_data.reserve(mesh.n_cells());

        auto Du = vol_vec.array() / (uEqn_diag.array() + prism::EPSILON);
        auto Dv = vol_vec.array() / (vEqn_diag.array() + prism::EPSILON);

        for (std::size_t i = 0; i < mesh.n_cells(); ++i) {
            // clang-format off
            Matrix3d Di;
            Di << Du[i], 0,     0,
                  0,     Dv[i], 0,
                  0,     0,     0;
            // clang-format on
            D_data.emplace_back(std::move(Di));
        }

        auto D = prism::field::Tensor("D", mesh, D_data);

        spdlog::info("Solving y-momentum equation");
        solver.solve(vEqn, 2, 1e-3, 0.9);

        spdlog::info("Solving x-momentum equation");
        solver.solve(uEqn, 2, 1e-3, 0.9);

        // Rhie-Chow interpolation for velocity face values
        ops::rhie_chow_correct(U, D, P);

        // pressure equation
        // Few things missing to implement:
        // 1) it's actually ρD not just D
        // 2) ∇.(ρU) not ∇.U
        auto P_prime = field::Scalar("pressure", mesh, 0.0);
        auto pEqn = TransportEquation(
            diffusion::Diffusion<field::Tensor, nonortho::OverRelaxedCorrector<>>(D, P_prime),
            source::Divergence<source::SourceSign::Negative>(U));

        pEqn.update_coeffs();
        spdlog::info("Solving pressure correction equation");
        solver.solve(pEqn, 10, 1e-5, 1);

        // update velocity fields
        auto p_grad = gradient::LeastSquares(P_prime);
        for (const auto& cell : mesh.cells()) {
            auto correction = -D.value_at_cell(cell) * p_grad.gradient_at_cell(cell);
            U.x()[cell.id()] += correction[0];
            U.y()[cell.id()] += correction[1];
        }

        P.data() = P.data().array() + (0.85 * P_prime.data().array());
    }
}