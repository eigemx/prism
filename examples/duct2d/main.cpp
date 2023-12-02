#include "fmt/core.h"
#include "prism/constants.h"
#include "prism/equation.h"
#include "prism/field.h"
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

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    if (argc < 2) {
        fmt::println("usage: {} [mesh-file]", argv[1]); // NOLINT
        return 1;
    }

    auto mesh = mesh::UnvToPMeshConverter(argv[1]).to_pmesh(); // NOLINT

    auto rho = ScalarField("density", mesh, 1.18);
    auto U = VectorField("velocity", mesh, Vector3d {0.05, 0.05, 0.0});
    auto P = ScalarField("pressure", mesh, 1.0);
    auto P_prime = ScalarField("pressure", mesh, 1.0);

    auto uEqn = TransportEquation(
        convection::Upwind(rho, U, U.x()),                                           // ∇.(ρUu)
        diffusion::Diffusion<double, nonortho::OverRelaxedCorrector<>>(1e-6, U.x()), // - ∇.(μ∇u)
        source::Gradient<source::SourceSign::Negative>(P, Coord::X)                  // ∂p/∂x
        //source::Laplacian(1e-6, U.x()) // - ∇.(μ∇u^T)
    );

    auto vEqn = TransportEquation(
        convection::Upwind(rho, U, U.y()),
        diffusion::Diffusion<double, nonortho::OverRelaxedCorrector<>>(1e-6, U.y()),
        source::Gradient<source::SourceSign::Negative>(P, Coord::Y)
        //source::Laplacian(1e-6, U.y())
    );

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

        auto D = prism::TensorField("D", mesh, D_data);

        fmt::println("Solving y-momentum equation");
        solver.solve(vEqn, 2, 1e-3, 0.9);

        fmt::println("Solving x-momentum equation");
        solver.solve(uEqn, 2, 1e-3, 0.9);

        // Rhie-Chow interpolation for velocity face values
        ops::rhie_chow_correct(U, D, P_prime);

        // pressure equation
        // Few things missing to implement:
        // 1) it's actually ρD not just D
        // 2) ∇.(ρU) not ∇.U
        auto pEqn = TransportEquation(
            diffusion::Diffusion<TensorField, nonortho::OverRelaxedCorrector<>>(D, P_prime),
            source::Divergence<source::SourceSign::Negative>(U));

        pEqn.update_coeffs();
        fmt::println("Solving pressure correction equation");
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