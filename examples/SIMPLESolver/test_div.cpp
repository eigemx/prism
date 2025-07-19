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
#include "prism/scheme/convection.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::scheme::convection;

template <typename Field>
void inspectCell(const mesh::Cell& cell, eqn::Transport<Field>& eqn, const field::Velocity& U) {
    log::info("Inspecting cell: {}", cell.id());
    log::info("Inspecting convection coefficients...");
    auto ac_conv = eqn.convectionScheme()->matrix().coeff(cell.id(), cell.id());
    auto b_conv = eqn.convectionScheme()->rhs().coeff(cell.id());
    auto ac_diff = eqn.diffusionScheme()->matrix().coeff(cell.id(), cell.id());
    auto b_diff = eqn.diffusionScheme()->rhs().coeff(cell.id());
    auto ac = eqn.matrix().coeff(cell.id(), cell.id());
    auto b = eqn.rhs().coeff(cell.id());
    log::info("ac_conv: {}, b_conv: {}", ac_conv, b_conv);
    log::info("ac_diff: {}, b_diff: {}", ac_diff, b_diff);
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
            auto neig = U.mesh()->otherSharingCell(cell, face);
            auto dx = (cell.center() - neig.center()).norm();
            auto uf = U.valueAtFace(face);
            auto mdot = uf.dot(sf);
            log::info("***************************************");
            log::info("face: {} is an interior face, owner to= {}, neighbor to {}",
                      face.id(),
                      face.owner(),
                      face.neighbor().value());
            log::info("uf = [{}, {}, {}]", uf.x(), uf.y(), uf.z());
            log::info("mdot: {}", mdot);
            log::info("sf = [{}, {}, {}]", sf.x(), sf.y(), sf.z());
            log::info("dx: {}", dx);

            if (mdot > 0) {
                log::info("ac = {}", mdot);
                ac_mine += mdot;
            } else {
                log::info("an = {}", -mdot);
            }

            log::info("an_diffusion = {}",
                      eqn.diffusionScheme()->matrix().coeff(cell.id(), neig.id()));
            log::info("an_convection = {}",
                      eqn.convectionScheme()->matrix().coeff(cell.id(), neig.id()));
            log::info("an_eqn = {}", eqn.matrix().coeff(cell.id(), neig.id()));
        }
    }
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

    auto U = field::Velocity("U", mesh, 0.0);
    auto T = field::Scalar("T", mesh, 0);
    auto rho = field::UniformScalar("rho", mesh, 1.0);
    auto nu = field::UniformScalar("nu", mesh, 1e-3);

    using div = Upwind<field::Velocity, field::VelocityComponent>;
    using laplacian = diffusion::Corrected<field::UniformScalar,
                                           diffusion::nonortho::OverRelaxedCorrector,
                                           field::VelocityComponent>;
    using grad = source::Gradient<Sign::Negative, field::Pressure>;

    auto momentum_solver = solver::BiCGSTAB<field::VelocityComponent>();
    auto p_solver = solver::BiCGSTAB<field::Pressure>();

    auto nNonOrthCorrectiors = 2;
    auto nOuterIter = 15;
    auto momentumURF = 1.0;
    auto pressureURF = 0.3;
    auto mDot = rho * U;

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();

    for (auto iter = 0; iter < nOuterIter; ++iter) {
        mDot = rho * U; // update mDot with the new U values

        log::info("Outer iteration: {}", iter);
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

        // auto kappa = field::UniformScalar("kappa", mesh, 1e-2);
        // using div2 = scheme::convection::Upwind<field::Velocity, field::Scalar>;
        // using laplacian2 = scheme::diffusion::NonCorrected<field::UniformScalar,
        // field::Scalar>;
        //  auto eqn = eqn::Transport(div2(U, T),          // ∇.(ρUT)
        //                            laplacian2(kappa, T) // - ∇.(κ ∇T)
        //);
        //  solver.solve(eqn, 10, 1e-20);
    }
    export_field_vtu(U.x(), "solution_x.vtu");
    export_field_vtu(U.y(), "solution_y.vtu");
    prism::export_field_vtu(T, "T.vtu");

    // export_field_vtu(ops::grad(p, Coord::X), "gradP_x.vtu");
    // export_field_vtu(ops::grad(p, Coord::Y), "gradP_y.vtu");
}
