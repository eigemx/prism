#include <fmt/base.h>
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
#include "prism/scheme/source.h"

using namespace prism;
using namespace prism::scheme;
using namespace prism::scheme::convection;


void inspectGradient(const mesh::Cell& cell, const field::Pressure& P) {
    Vector3d grad {.0, .0, .0};
    fmt::println("Value of P at this cell = {}", P.valueAtCell(cell));
    for (const auto face_id : cell.facesIds()) {
        const auto& face = P.mesh()->face(face_id);

        if (face.isBoundary()) {
            const auto& patch = P.mesh()->boundaryPatch(face);
            if (patch.isEmpty()) {
                continue;
            }

            fmt::println("Face {} is a boundary face on patch {}", face_id, patch.name());
            fmt::println("Value of face at this boundary = {}", P.valueAtFace(face));
            grad += P.valueAtFace(face) * face.areaVector();
        } else {
            grad += P.valueAtFace(face) * mesh::outwardAreaVector(face, cell);
            fmt::println("Value of face at this interior = {}", P.valueAtFace(face));
        }
    }
    grad /= cell.volume();
    fmt::println("grad at cell = [{}, {}, {}]", grad.x(), grad.y(), grad.z());
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

    using div = LinearUpwind<field::Velocity, field::VelocityComponent>;
    using laplacian = diffusion::Corrected<field::UniformScalar,
                                           diffusion::nonortho::OverRelaxedCorrector,
                                           field::VelocityComponent>;
    using grad = source::Gradient<Sign::Negative, field::Pressure>;
    using laplacian_s =
        source::Laplacian<Sign::Positive, field::UniformScalar, field::VelocityComponent>;

    auto momentum_solver = solver::BiCGSTAB<field::VelocityComponent>();
    auto p_solver = solver::BiCGSTAB<field::Pressure>();

    auto nNonOrthCorrectors = 3;
    auto nOuterIter = 30;
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
                                                        // laplacian_s(nu, U.x()) // + ∇.(U∇u)
        );

        auto vEqn = eqn::Momentum(div(mDot, U.y()),     // ∇.(Uv)
                                  laplacian(nu, U.y()), // -∇.(ν∇v)
                                  grad(p, Coord::Y)     // = -∂p/∂y
                                                        // laplacian_s(nu, U.y()) // + ∇.(U∇v)
        );

        uEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
        uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
        uEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

        vEqn.boundaryHandlersManager().addHandler<eqn::boundary::NoSlip<eqn::Momentum>>();
        vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Symmetry<eqn::Momentum>>();
        vEqn.boundaryHandlersManager().addHandler<eqn::boundary::Outlet<eqn::Momentum>>();

        // uEqn.setUnderRelaxFactor(momentumURF);
        // vEqn.setUnderRelaxFactor(momentumURF);

        for (auto n = 0; n < nNonOrthCorrectors; ++n) {
            log::info("Solving x-momentum equations");
            momentum_solver.solve(uEqn, 15, 1e-20);

            log::info("Solving y-momentum equations");
            momentum_solver.solve(vEqn, 15, 1e-20);
        }

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
    prism::export_field_vtu(p, "p.vtu");

    export_field_vtu(ops::grad(p, Coord::X), "gradP_x.vtu");
    export_field_vtu(ops::grad(p, Coord::Y), "gradP_y.vtu");

    inspectGradient(mesh->cell(900), p);
}
