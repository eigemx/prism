#include <fmt/core.h>
#include <prism/prism.h>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>
#include <stdexcept>

#include "prism/constants.h"
#include "prism/field/scalar.h"

using json = nlohmann::json;
using namespace prism;

namespace fs = std::filesystem;

struct FoamFields {
    std::vector<prism::Vector3d> velocity;
    std::vector<double> pressure;
    std::vector<double> temperature;
};

auto fileToJson(const fs::path& path) -> json {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File " + path.string() + " does not exist!");
    }
    auto file = std::ifstream(path);
    return json::parse(file);
}

auto readFields(const std::filesystem::path& fields_file) -> FoamFields {
    auto doc = fileToJson(fields_file);

    auto velocity = doc["U"];
    auto pressure = doc["p"].get<std::vector<double>>();

    std::vector<prism::Vector3d> velocity_vec;
    for (const auto& v : velocity) {
        velocity_vec.emplace_back(v[0], v[1], v[2]);
    }

    return {velocity_vec, pressure, {}};
}

auto readVelocityComponents(const std::vector<prism::Vector3d>& velocity)
    -> std::array<VectorXd, 3> {
    VectorXd u, v, w; // NOLINT
    u.resize(velocity.size());
    v.resize(velocity.size());
    w.resize(velocity.size());

    for (size_t i = 0; i < velocity.size(); i++) {
        u[i] = velocity[i].x();
        v[i] = velocity[i].y();
        w[i] = velocity[i].z();
    }

    return {u, v, w};
}

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
    auto mu = field::UniformScalar("mu", mesh, 1e-5);
    // auto U = field::Velocity("U", mesh, {0.0, 0.0, 0.0});
    auto U = field::Velocity("U", mesh, components);
    auto P = field::Pressure("P", mesh, 0.0);
    auto rho = field::UniformScalar("rho", mesh, 1.0);

    using div = Upwind<field::UniformScalar, field::VelocityComponent>;
    using laplacian = diffusion::NonCorrected<field::UniformScalar, field::VelocityComponent>;
    using grad = source::Gradient<source::SourceSign::Negative, field::Pressure>;

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

    auto U_solver = solver::BiCGSTAB<field::VelocityComponent>();

    auto nNonOrthCorrectiors = 4;
    for (auto nOuterIter = 0; nOuterIter < 10; ++nOuterIter) {
        /*
        log::info("Outer iteration {}", nOuterIter);

        for (auto i = 0; i < nNonOrthCorrectiors; ++i) {
            log::info("Solving y-momentum equations");
            U_solver.solve(vEqn, 5, 1e-9);

            log::info("Solving x-momentum equations");
            U_solver.solve(uEqn, 5, 1e-9);
        }
        */

        uEqn.updateCoeffs();
        vEqn.updateCoeffs();
        // uEqn.relax();
        // vEqn.relax();

        // calculate coefficients for the pressure equation
        const auto& vol_vec = mesh.cellsVolumeVector();
        const auto& uEqn_diag = uEqn.matrix().diagonal();
        const auto& vEqn_diag = vEqn.matrix().diagonal();

        auto D_data = std::vector<Matrix3d>();
        D_data.resize(mesh.cellCount());

        auto Du = vol_vec.array() / (uEqn_diag.array() + EPSILON);
        auto Dv = vol_vec.array() / (vEqn_diag.array() + EPSILON);

        for (const auto& cell : mesh.cells()) {
            auto i = cell.id();
            // clang-format off
            Matrix3d Di;
            Di << Du[i], 0,     0,
                  0,     Dv[i], 0,
                  0,     0,     0;
            // clang-format on
            D_data[i] = Di;
        }

        auto D = field::Tensor("D", mesh, D_data);

        // Rhie-Chow interpolation for velocity face values
        log::info("Correcting faces velocities using Rhie-Chow interpolation");
        auto U_rh = ops::rhieChowCorrect(U, D, P);

        auto divU_rh = ops::div(U_rh);
        auto divU = ops::div(U);
        export_field_vtu(divU, "divU.vtu");
        export_field_vtu(divU_rh, "divU_rh.vtu");

        // pressure correction field created with same name as pressure field to get same boundary
        // conditions without having to define P_prime in fields.json file.
        auto P_prime = field::Pressure("P", mesh, 0.0);

        // NOTE: The corrector should reset to zero the correction ﬁeld at every iteration and
        // should also apply a zero value at all boundaries for which a Dirichlet (fixed) boundary
        // condition is used for the pressure.


        // pressure correction equation (density is dropped from both sides of the equation due to
        // incompressibility assumption)
        using laplacian_p = diffusion::NonCorrected<field::Tensor, field::Pressure>;
        using div_U = source::Divergence<source::SourceSign::Negative, field::Velocity>;

        auto pEqn = eqn::Transport<field::Pressure>(laplacian_p(D, P_prime), // - ∇.(D ∇P_prime)
                                                    div_U(U_rh)              // == - (∇.U)
        );
        auto p_solver = solver::BiCGSTAB<field::Pressure>();

        log::info("Solving pressure correction equation");
        for (auto i = 0; i < nNonOrthCorrectiors; ++i) {
            p_solver.solve(pEqn, 3, 1e-16);
        }
        export_field_vtu(pEqn.field(), "pressure_correction.vtu");

        // update velocity fields
        /*
        for (const auto& cell : mesh.cells()) {
            prism::Vector3d correction = -D.valueAtCell(cell) * P_prime.gradAtCell(cell);
            U.x()[cell.id()] += correction.x();
            U.y()[cell.id()] += correction.y();
        }
        */
        P.values() = P.values().array() + (0.75 * P_prime.values().array());

        uEqn.zeroOutCoeffs();
        vEqn.zeroOutCoeffs();
    }

    export_field_vtu(U.x(), "solution_x.vtu");
    export_field_vtu(U.y(), "solution_y.vtu");
    export_field_vtu(P, "pressure.vtu");
    export_field_vtu(components[0], "solution_x_true.vtu");

    auto diff = field::Scalar("diff", mesh, components[0].values() - U.x().values());
    export_field_vtu(diff, "diff.vtu");

    auto vol_field = field::Scalar("volume", mesh, mesh.cellsVolumeVector());
    export_field_vtu(vol_field, "volume.vtu");
}
