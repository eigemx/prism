#include <fmt/core.h>
#include <prism/prism.h>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

#include "prism/scheme/nonortho.h"

using json = nlohmann::json;
using prism::export_field_vtu;
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

auto main(int argc, char* argv[]) -> int {
    using namespace prism::scheme;

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

    // read pressure field from json file
    auto fields = readFields(fs::path(unv_file_name).parent_path() / "foam_fields.json");
    auto pressure_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(fields.pressure.data(),
                                                                      fields.pressure.size());

    // set mesh fields
    auto mu = field::UniformScalar("mu", mesh, 1e-5);
    auto U = field::Velocity("U", mesh, {0.0, 0.0, 0.0});
    auto P = field::Pressure("P", mesh, pressure_vec);
    auto rho = field::UniformScalar("rho", mesh, 1.0);

    using div = convection::Upwind<field::UniformScalar, field::VelocityComponent>;
    using laplacian = diffusion::
        Corrected<field::UniformScalar, nonortho::OverRelaxedCorrector, field::VelocityComponent>;
    using grad = source::Gradient<scheme::source::SourceSign::Negative, field::Pressure>;

    auto uEqn = eqn::Momentum(div(rho, U, U.x()),   // ∇.(ρUu)
                              laplacian(mu, U.x()), // - ∇.(μ∇u)
                              grad(P, Coord::X)     // = -∂p/∂x
    );

    auto vEqn = eqn::Momentum(div(rho, U, U.y()),   // ∇.(ρUv)
                              laplacian(mu, U.y()), // - ∇.(μ∇v)
                              grad(P, Coord::Y)     // = -∂p/∂y
    );


    auto U_solver = solver::BiCGSTAB<field::VelocityComponent,
                                     solver::ImplicitUnderRelaxation<field::VelocityComponent>>();
    for (auto i = 0; i < 10; ++i) {
        log::info("Solving y-momentum equation");
        U_solver.solve(vEqn, 10, 1e-4, 1.0);

        log::info("Solving x-momentum equation");
        U_solver.solve(uEqn, 10, 1e-4, 1.0);
    }

    export_field_vtu(uEqn.field(), "solution_x.vtu");
    export_field_vtu(vEqn.field(), "solution_y.vtu");
    export_field_vtu(ops::grad(P, Coord::X), "gradP_x.vtu");
    /*
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
    }*/
}
