#include <fmt/base.h>
#include <prism/prism.h>

#include <filesystem>

#include "helpers.h"
#include "prism/constants.h"
#include "prism/export.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/log.h"
#include "prism/operations/rhie_chow.h"

using namespace prism;

auto main(int argc, char* argv[]) -> int {
    using namespace prism::scheme;
    using namespace prism::scheme::convection;

    log::setLevel(log::Level::Info);
    if (argc < 2) {
        fmt::println("usage: {} [mesh-file]", argv[1]); // NOLINT
        return 1;
    }

    const auto* unv_file_name = argv[1]; // NOLINT

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // read pressure field
    auto fields = readFields(fs::path(unv_file_name).parent_path() / "foam_fields.json");
    auto pressure_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(fields.pressure.data(),
                                                                      fields.pressure.size());
    // read velocity field
    auto raw_velocity_array = readVelocityComponents(fields.velocity);
    auto Ux = field::VelocityComponent("U_x", mesh, raw_velocity_array[0], Coord::X);
    auto Uy = field::VelocityComponent("U_y", mesh, raw_velocity_array[1], Coord::Y);
    auto Uz = field::VelocityComponent("U_z", mesh, raw_velocity_array[2], Coord::Z);
    auto components = std::array {Ux, Uy, Uz};

    // set mesh fields
    auto nu = field::UniformScalar("mu", mesh, 1e-5);
    auto U = field::Velocity("U", mesh, {0.0, 0.0, 0.0});
    auto p = field::Pressure("P", mesh, pressure_vec);
    auto rho = field::UniformScalar("rho", mesh, 1.0);
    // rhoU = rho * U;

    using div = Upwind<field::Velocity, field::VelocityComponent>;
    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>;
    using grad = source::Gradient<Sign::Negative, field::Pressure>;

    auto momentumSolver = solver::BiCGSTAB<field::VelocityComponent>();

    auto nNonOrthCorrectiors = 3;

    for (auto nOuterIter = 0; nOuterIter < 1; ++nOuterIter) {
        auto mDot = rho * U; // remove this in true run

        auto uEqn = eqn::Momentum(div(mDot, U.x()),     // ∇.(Uu)
                                  laplacian(nu, U.x()), // -∇.(ν∇u)
                                  grad(p, Coord::X)     // = -∂p/∂x
        );

        auto vEqn = eqn::Momentum(div(mDot, U.y()),     // ∇.(Uv)
                                  laplacian(nu, U.y()), // -∇.(ν∇v)
                                  grad(p, Coord::Y)     // = -∂p/∂y
        );

        // uEqn.setUnderRelaxFactor(0.9);
        // vEqn.setUnderRelaxFactor(0.9);
        // uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
        // vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();

        log::info("Outer iteration {}", nOuterIter);
        log::info("Solving y-momentum equations");
        momentumSolver.solve(vEqn, 5, 1e-16);
        log::info("Solving x-momentum equations");
        momentumSolver.solve(uEqn, 5, 1e-16);

        uEqn.updateCoeffs();
        vEqn.updateCoeffs();

        writeToCSVfile("A.csv", Eigen::MatrixXd(uEqn.matrix()));
        writeToCSVfile("B.csv", Eigen::MatrixXd(vEqn.matrix()));
        writeToCSVfile("diffusion_matrix.csv", Eigen::MatrixXd(uEqn.diffusionScheme()->matrix()));
        writeToCSVfile("convection_matrix.csv",
                       Eigen::MatrixXd(vEqn.convectionScheme()->matrix()));

        /*
        export_field_vtu(U.x(), "solution_x_b.vtu");
        export_field_vtu(U.y(), "solution_y_b.vtu");

        uEqn.updateCoeffs();
        uEqn.relax();
        vEqn.updateCoeffs();
        vEqn.relax();

        // calculate coefficients for the pressure equation
        const auto& vol_vec = mesh->cellsVolumeVector();
        const auto& uEqn_diag = uEqn.matrix().diagonal();
        const auto& vEqn_diag = vEqn.matrix().diagonal();

        auto D_data = std::vector<Matrix3d>();
        D_data.resize(mesh->cellCount());

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
        log::info("Correcting faces velocities using Rhie-Chow interpolation");
        auto U_corrected = U.clone();
        U_corrected.updateInteriorFaces(
            [&](const mesh::Face& face) { return ops::rhieChowCorrectFace(face, U, D, P); });

        // pressure correction field created with same name as pressure field to get same boundary
        // conditions without having to define P_prime in fields.json file.
        auto P_prime = field::Pressure("P", mesh, 0.0);

        // NOTE: The corrector should reset to zero the correction ﬁeld at every iteration and
        // should also apply a zero value at all boundaries for which a Dirichlet (fixed) boundary
        // condition is used for the pressure.

        // pressure correction equation (density is dropped from both sides of the equation due to
        // incompressibility assumption)
        using laplacian_p =
            diffusion::Corrected<field::Tensor,
                                 scheme::diffusion::nonortho::OverRelaxedCorrector,
                                 field::Pressure>;
        using div_U = source::Divergence<Sign::Negative, field::Velocity>;

        auto pEqn = eqn::Transport<field::Pressure>(laplacian_p(D, P_prime), // - ∇.(D ∇P_prime)
                                                    div_U(U_corrected)       // == - (∇.U)
        );
        auto p_solver = solver::BiCGSTAB<field::Pressure>();

        log::info("Solving pressure correction equation");
        for (auto i = 0; i < nNonOrthCorrectiors; ++i) {
            p_solver.solve(pEqn, 3, 1e-16);
        }
        export_field_vtu(pEqn.field(), "pressure_correction.vtu");


        // update velocity fields
        U.updateCells([&](const mesh::Cell& cell) {
            return U.valueAtCell(cell) - (D.valueAtCell(cell) * P_prime.gradAtCell(cell));
        });

        // update mass flow rate at faces
        U_corrected.updateInteriorFaces([&](const mesh::Face& face) {
            return U_corrected.valueAtFace(face) -
                   (D.valueAtFace(face) * P_prime.gradAtFace(face));
        });

        // update pressure
        P.values() = P.values().array() + (0.75 * P_prime.values().array());

        rhoU = rho * U_corrected;

        uEqn.zeroOutCoeffs();
        vEqn.zeroOutCoeffs();
        */
    }

    export_field_vtu(U.x(), "solution_x.vtu");
    export_field_vtu(U.y(), "solution_y.vtu");
    export_field_vtu(p, "pressure.vtu");
    auto diff = field::Scalar("diff", mesh, components[0].values() - U.x().values());
    export_field_vtu(diff, "diff.vtu");
}
