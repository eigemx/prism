#include <fmt/core.h>
#include <foamJSONToPMesh.h>
#include <prism/prism.h>

#include <filesystem>

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    using namespace prism::scheme;


    log::setLevel(log::Level::Info);
    if (argc < 2) {
        fmt::println("usage: {} [mesh-file]", argv[1]); // NOLINT
        return 1;
    }

    const auto* unv_file_name = argv[1]; // NOLINT

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    fmt::println("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    auto mu = field::UniformScalar("mu", mesh, 1e-3);
    auto U = field::Velocity("U", mesh, {0.01, 0.01, 0});
    auto P = field::Pressure("P", mesh, 0.0);
    auto rho = field::UniformSalar("rho", mesh, 1000);

    using div = convection::Upwind<field::UniformSalar, field::VelocityComponent>;
    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>;
    using grad = source::Gradient<scheme::source::SourceSign::Negative, field::Pressure>;

    auto uEqn = eqn::Momentum(div(rho, U, U.x()),   // ∇.(ρUu)
                              laplacian(mu, U.x()), // - ∇.(μ∇u)
                              grad(P, Coord::X)     // = -∂p/∂x
    );

    auto vEqn = eqn::Momentum(div(rho, U, U.y()),   // ∇.(ρUv)
                              laplacian(mu, U.y()), // - ∇.(μ∇v)
                              grad(P, Coord::Y)     // = -∂p/∂y
    );

    uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
    vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();


    auto U_solver = solver::BiCGSTAB<field::VelocityComponent,
                                     solver::ImplicitUnderRelaxation<field::VelocityComponent>>();

    auto p_solver =
        solver::BiCGSTAB<field::Pressure, solver::ImplicitUnderRelaxation<field::Pressure>>();

    for (auto i = 0; i < 2; ++i) {
        log::info("Solving x-momentum equation");
        U_solver.solve(uEqn, 10, 1e-4, 0.8);

        log::info("Solving y-momentum equation");
        U_solver.solve(vEqn, 10, 1e-4, 0.8);

        uEqn.zeroOutCoeffs();
        vEqn.zeroOutCoeffs();

        uEqn.updateCoeffs();
        vEqn.updateCoeffs();

        // calculate coefficients for the pressure equation
        const auto& vol_vec = mesh.cellsVolumeVector();
        const auto& uEqn_diag = uEqn.matrix().diagonal();
        const auto& vEqn_diag = vEqn.matrix().diagonal();

        auto D_data = std::vector<Matrix3d>();
        D_data.reserve(mesh.cellCount());

        auto Du = vol_vec.array() / (uEqn_diag.array() + prism::EPSILON);
        auto Dv = vol_vec.array() / (vEqn_diag.array() + prism::EPSILON);

        for (const auto& cell : mesh.cells()) {
            auto i = cell.id();
            // clang-format off
            Matrix3d Di;
            Di << Du[i], 0,     0,
                  0,     Dv[i], 0,
                  0,     0,     0;
            // clang-format on
            D_data.emplace_back(Di);
        }

        auto D = prism::field::Tensor("D", mesh, D_data);

        // Rhie-Chow interpolation for velocity face values
        ops::correctRhieChow(U, D, P);

        // pressure correction field created with same name as pressure field to get same boundary
        // conditions without having to define P_prime in fields.json file.
        auto P_prime = field::Pressure("P", mesh, 0.0);

        using laplacian_p = diffusion::NonCorrected<field::Tensor, field::Pressure>;
        using div_p = source::Divergence<source::SourceSign::Negative, field::Velocity>;

        // pressure correction equation (density is dropped from both sides of the equation due to
        // incompressibility assumption)
        auto pEqn = eqn::Transport<field::Pressure>(laplacian_p(D, P_prime), // - ∇.(D ∇P_prime)
                                                    div_p(U)                 // == - (∇.U)
        );

        log::info("Solving pressure correction equation");
        p_solver.solve(pEqn, 5, 1e-5, 1.0);

        // update velocity fields
        for (const auto& cell : mesh.cells()) {
            auto correction = D.valueAtCell(cell) * P_prime.gradAtCell(cell);
            U.x()[cell.id()] += -correction[0];
            U.y()[cell.id()] += -correction[1];
        }

        P.values() = P.values().array() + (0.85 * P_prime.values().array());

        uEqn.zeroOutCoeffs();
        vEqn.zeroOutCoeffs();
    }
}
