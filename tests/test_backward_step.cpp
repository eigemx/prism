#include <prism/field/velocity.h>
#include <prism/prism.h>

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using std::filesystem::path;
using namespace prism;

struct FoamFields {
    std::vector<prism::Vector3d> velocity;
    std::vector<double> pressure;
    std::vector<double> temperature;
};

auto fileToJson(const path& path) -> json {
    // check if file exists
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File " + path.string() + " does not exist!");
    }
    auto file = std::ifstream(path);
    return json::parse(file);
}

auto readFields(const std::filesystem::path& fields_file) -> FoamFields {
    /**
     * Reads the internal fields from the json file "inetrnal_fields.json" in the directory
     * "tests/cases/pitzDailyFoam/" and returns them as a PitzDailyFields struct.
     */
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
    VectorXd u, w, v; // NOLINT
    u.resize(velocity.size());
    w.resize(velocity.size());
    v.resize(velocity.size());

    for (size_t i = 0; i < velocity.size(); i++) {
        u[i] = velocity[i].x();
        w[i] = velocity[i].y();
        v[i] = velocity[i].z();
    }

    return {u, w, v};
}

auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    log::setLevel(log::Level::Debug);

    if (argc < 2) {
        fmt::println("Usage: convsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = argv[1]; // NOLINT

    // read mesh
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = field::Scalar("T", mesh, 300.0);

    // density field
    auto rho = field::UniformScalar("rho", mesh, 1);

    // velocity field
    auto foam_fields = std::filesystem::path(unv_file_name).parent_path() / "foam_fields.json";
    auto fields = readFields(foam_fields);
    auto raw_velocity_array = readVelocityComponents(fields.velocity);
    auto Ux = field::VelocityComponent("U_x", mesh, raw_velocity_array[0]);
    auto Uy = field::VelocityComponent("U_y", mesh, raw_velocity_array[1]);
    auto Uz = field::VelocityComponent("U_z", mesh, raw_velocity_array[2]);
    auto components = std::array {Ux, Uy, Uz};
    auto U = field::Velocity("U", mesh, components);

    // diffusion coefficient
    auto kappa = field::UniformScalar("kappa", mesh, 1e-2);

    // solve for temperature advection: ∇.(ρUT) - ∇.(κ ∇T) = 0
    // where ρ is the density and U is the velocity vector, and S is an arbitraty constant source
    using div = scheme::convection::QUICK<field::UniformScalar, field::Scalar>;
    using laplacian = scheme::diffusion::NonCorrected<field::UniformScalar, field::Scalar>;

    auto eqn = eqn::Transport(div(rho, U, T),     // ∇.(ρUT)
                              laplacian(kappa, T) // - ∇.(κ ∇T)
    );

    // solve
    auto solver = solver::BiCGSTAB<field::Scalar>();
    solver.solve(eqn, 5, 1e-5, 1);

    prism::export_field_vtu(eqn.field(), "solution.vtu");
    auto div_U = ops::div(U);
    prism::export_field_vtu(div_U, "div.vtu");
}
