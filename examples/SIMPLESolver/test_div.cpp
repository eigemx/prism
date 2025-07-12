#include <prism/prism.h>

#include <filesystem>

#include "helpers.h"
#include "prism/equation/transport.h"
#include "prism/export.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/log.h"
#include "prism/mesh/utilities.h"
#include "prism/numerics/solver.h"
#include "prism/operations/operations.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::scheme::convection;

template <typename Field>
void inspectCell(const mesh::Cell& cell, eqn::Transport<Field>& eqn, const field::Velocity& U) {
    log::info("Inspecting cell: {}", cell.id());
    log::info("Inspecting convection coefficients...");
    auto ac = eqn.convectionScheme()->matrix().coeff(cell.id(), cell.id());
    auto b = eqn.convectionScheme()->rhs().coeff(cell.id());
    log::info("ac: {}, b: {}", ac, b);

    double ac_mine = 0.0;
    double b_mine = 0.0;
    for (const auto face_id : cell.facesIds()) {
        const auto& face = U.mesh()->face(face_id);

        if (face.isBoundary()) {
            auto patch = U.mesh()->boundaryPatch(face);
            if (patch.isEmpty()) {
                continue;
            }

            if (patch.name() == "UpperWall") {
                continue;
            }
            auto mdot_f = U.fluxAtFace(face);
            auto phi_wall = patch.getScalarBoundaryCondition("T");
            ac_mine += std::max(mdot_f, 0.0);
            b_mine += std::max(-mdot_f, 0.0) * phi_wall;

        } else {
            auto sf = mesh::outwardAreaVector(face, cell);
            auto uf = U.valueAtFace(face);
            auto mdot = uf.dot(sf);
            log::info("***************************************");
            log::info("face: {} is an interior face, owner to= {}, neighbor to {}",
                      face.id(),
                      face.owner(),
                      face.neighbor().value());
            log::info("uf = [{}, {}, {}]", uf.x(), uf.y(), uf.z());
            log::info("mdot: {}", mdot);

            if (mdot > 0) {
                ac_mine += mdot;
            }
        }
    }
    log::info("ac_mine: {}, b_mine: {}", ac_mine, b_mine);
}

auto main(int argc, char* argv[]) -> int {
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

    auto foam_fields =
        readFields(std::filesystem::path(unv_file_name).parent_path() / "foam_fields.json");

    auto U_raw_components = readVelocityComponents(foam_fields.velocity);
    auto U_x = field::VelocityComponent("U_x", mesh, U_raw_components[0]);
    auto U_y = field::VelocityComponent("U_y", mesh, U_raw_components[1]);
    auto U_z = field::VelocityComponent("U_z", mesh, U_raw_components[2]);
    auto components = std::array {U_x, U_y, U_z};

    // convert pressure from std::vector<double> to Eigen::VectorXd
    auto p_vec = VectorXd::Map(foam_fields.pressure.data(), foam_fields.pressure.size());
    auto p = field::Pressure("P", mesh, p_vec);
    auto U = field::Velocity("U", mesh, components);
    auto T = field::Scalar("T", mesh, 0.0);
    auto rho = field::UniformScalar("rho", mesh, 1.0);
    auto nu = field::UniformScalar("nu", mesh, 1e-3);

    using div = LinearUpwind<field::Velocity, field::VelocityComponent>;
    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>;
    using grad = source::Gradient<Sign::Negative, field::Pressure>;

    auto momentum_solver = solver::BiCGSTAB<field::VelocityComponent>();
    auto p_solver = solver::BiCGSTAB<field::Pressure>();

    auto nNonOrthCorrectiors = 2;
    auto nOuterIter = 1;
    auto momentumURF = 1.0;
    auto pressureURF = 0.3;
    auto mDot = rho * U;

    auto kappa = field::UniformScalar("kappa", mesh, 4e-5);
    using div2 = scheme::convection::Upwind<field::Velocity, field::Scalar>;
    using laplacian2 = scheme::diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport(div2(mDot, T),       // ∇.(ρUT)
                              laplacian2(kappa, T) // - ∇.(κ ∇T)
    );


    eqn.updateCoeffs();
    writeToCSVfile("conv_matrix.csv", Eigen::MatrixXd(eqn.convectionScheme()->matrix()));
    writeToCSVfile("conv_b.csv", eqn.convectionScheme()->rhs());
    auto b_field = field::Scalar("b", mesh, eqn.convectionScheme()->rhs());
    export_field_vtu(b_field, "b.vtu");
    inspectCell(mesh->cell(110), eqn, U);
    eqn.zeroOutCoeffs();


    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();

    for (auto outer_iteration = 0; outer_iteration < nOuterIter; ++outer_iteration) {
        /*
        mDot = rho * U; // update mDot with the new U values

        log::info("Outer iteration: {}", outer_iteration);
        auto uEqn = eqn::Momentum(div(mDot, U.x()),     // ∇.(Uu)
                                  laplacian(nu, U.x()), // -∇.(ν∇u)
                                  grad(p, Coord::X)     // = -∂p/∂x
        );

        auto vEqn = eqn::Momentum(div(mDot, U.y()),     // ∇.(Uv)
                                  laplacian(nu, U.y()), // -∇.(ν∇v)
                                  grad(p, Coord::Y)     // = -∂p/∂y
        );

        uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
        uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
        vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
        vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();

        // uEqn.setUnderRelaxFactor(momentumURF);
        // vEqn.setUnderRelaxFactor(momentumURF);

        log::info("Solving x-momentum equations");
        momentum_solver.solve(uEqn, 15, 1e-20);

        log::info("Solving y-momentum equations");
        momentum_solver.solve(vEqn, 15, 1e-20);
        */


        solver.solve(eqn, 10, 1e-20);
    }
    export_field_vtu(U.x(), "solution_x.vtu");
    export_field_vtu(U.y(), "solution_y.vtu");
    prism::export_field_vtu(T, "T.vtu");

    export_field_vtu(ops::grad(p, Coord::X), "gradP_x.vtu");
    export_field_vtu(ops::grad(p, Coord::Y), "gradP_y.vtu");
}
