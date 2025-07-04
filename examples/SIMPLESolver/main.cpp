#include <prism/prism.h>

#include <filesystem>

/// TODO: The corrector should reset to zero the correction ﬁeld at every iteration and
/// should also apply a zero value at all boundaries for which a Dirichlet (fixed) boundary
/// condition is used for the pressure.
using namespace prism;

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
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set mesh fields
    auto U = field::Velocity("U", mesh, {.0, .0, .0});
    auto p = field::Pressure("P", mesh, 0.0);
    auto rho = field::UniformScalar("rho", mesh, 1.0);
    auto nu = field::UniformScalar("nu", mesh, 1e-3);

    using div = LinearUpwind<field::Velocity, field::VelocityComponent>;
    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>;
    using grad = source::Gradient<Sign::Negative, field::Pressure>;

    auto momentum_solver = solver::BiCGSTAB<field::VelocityComponent>();
    auto p_solver = solver::BiCGSTAB<field::Pressure>();

    auto nNonOrthCorrectiors = 2;
    auto nOuterIter = 50;
    auto momentumURF = 0.7;
    auto pressureURF = 0.3;
    auto mDot = rho * U;

    for (auto outer_iteration = 0; outer_iteration < nOuterIter; ++outer_iteration) {
        log::info("Outer iteration: {}", outer_iteration);
        auto uEqn = eqn::Momentum(div(mDot, U.x()),     // ∇.(Uu)
                                  laplacian(nu, U.x()), // -∇.(ν∇u)
                                  grad(p, Coord::X)     // = -∂p/∂x
        );

        auto vEqn = eqn::Momentum(div(mDot, U.y()),     // ∇.(Uv)
                                  laplacian(nu, U.y()), // -∇.(ν∇v)
                                  grad(p, Coord::Y)     // = -∂p/∂y
        );

        uEqn.setUnderRelaxFactor(momentumURF);
        vEqn.setUnderRelaxFactor(momentumURF);

        log::info("Solving x-momentum equations");
        momentum_solver.solve(uEqn, 15, 1e-20);

        log::info("Solving y-momentum equations");
        momentum_solver.solve(vEqn, 15, 1e-20);

        uEqn.updateCoeffs();
        uEqn.relax();
        vEqn.updateCoeffs();
        vEqn.relax();

        // calculate coefficients for the pressure equation
        const VectorXd& vol_vec = mesh->cellsVolumeVector();
        const VectorXd& uEqn_diag = uEqn.matrix().diagonal();
        const VectorXd& vEqn_diag = vEqn.matrix().diagonal();

        auto D_data = std::vector<Matrix3d>(mesh->cellCount(), Matrix3d::Zero());

        auto Du = vol_vec.array() / (uEqn_diag.array() + EPSILON);
        auto Dv = vol_vec.array() / (vEqn_diag.array() + EPSILON);

        for (const auto& cell : mesh->cells()) {
            auto i = cell.id();
            // clang-format off
            Matrix3d Di;
            Di << Du[i], 0,     0,
                  0,     Dv[i], 0,
                  0,     0,     1;
            // clang-format on
            D_data[i] = Di;
        }

        auto D = field::Tensor("D", mesh, D_data);

        // Rhie-Chow interpolation for velocity face values
        log::debug("Correcting faces velocities using Rhie-Chow interpolation");
        mDot = rho * U; // update mDot with the new U values

        mDot.updateInteriorFaces(
            [&](const mesh::Face& face) { return ops::rhieChowCorrectFace(face, mDot, D, p); });

        // pressure correction field created with same name as pressure field to get same boundary
        // conditions without having to define P_prime in fields.json file.
        auto P_prime = field::Pressure("P", mesh, 0.0);

        using laplacian_p = diffusion::
            Corrected<field::Tensor, diffusion::nonortho::OverRelaxedCorrector, field::Pressure>;
        using div_U = source::Divergence<Sign::Negative, field::Velocity>;

        auto pEqn = eqn::Transport<field::Pressure>(laplacian_p(D, P_prime), // - ∇.(D ∇P_prime)
                                                    div_U(mDot)              // == - (∇.U)
        );

        log::info("Solving pressure correction equation");
        p_solver.solve(pEqn, 3, 1e-16);
        for (auto i = 0; i < nNonOrthCorrectiors; ++i) {
            p_solver.solve(pEqn, 3, 1e-16);
        }

        // For some reason, the following code produces a different result (and wrong) if we used
        // U.valueAtCell directly after the return statement. It seems that we must first evaluate
        // the value of U at the cell then use it.
        /// TODO: fix this.

        // update velocity field
        U.updateCells([&](const mesh::Cell& cell) {
            auto U_cell = U.valueAtCell(cell);
            return U_cell - (D.valueAtCell(cell) * P_prime.gradAtCell(cell));
        });

        // update mass flow rate at faces
        mDot.updateInteriorFaces([&](const mesh::Face& face) {
            return mDot.valueAtFace(face) - (D.valueAtFace(face) * P_prime.gradAtFace(face));
        });

        // update pressure
        p.values() = p.values().array() + (pressureURF * P_prime.values().array());

        uEqn.zeroOutCoeffs();
        vEqn.zeroOutCoeffs();
    }

    export_field_vtu(U.x(), "solution_x.vtu");
    export_field_vtu(U.y(), "solution_y.vtu");
    export_field_vtu(p, "pressure.vtu");
}
