#include <fmt/core.h>
#include <prism/prism.h>

#include <fstream>
#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;
namespace fs = std::filesystem;

auto fileToJson(const fs::path& path) -> json {
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File " + path.string() + " does not exist!");
    }
    auto file = std::ifstream(path);
    return json::parse(file);
}

auto readFields(const std::filesystem::path& fields_file) -> std::vector<double> {
    auto doc = fileToJson(fields_file);
    auto temp = doc["T"];
    return temp;
}

auto jsonToPrismField(const fs::path& fields_file, const prism::mesh::PMesh& mesh)
    -> prism::field::Scalar {
    auto doc = fileToJson(fields_file);
    auto temp = doc["T"].get<std::vector<double>>();
    auto temp_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(temp.data(), temp.size());

    return {"T", mesh, temp_vec};
}


auto main(int argc, char* argv[]) -> int {
    using namespace prism;
    using namespace prism::scheme;
    using namespace prism::field;

    log::setLevel(log::Level::Info);

    log::info("diffsolver - A steady state diffusion equation solver");

    if (argc < 2) {
        log::error("Usage: diffsolver [mesh-file]");
        return 1;
    }

    auto unv_file_name = argv[1]; // NOLINT

    // read mesh
    log::info("Loading mesh file `{}`...", unv_file_name);
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    // set up the temperature field defined over the mesh, with an initial value of 300.0 [K]
    auto T = Scalar("T", mesh, 300.0);

    // define a source term
    VectorXd source_field_data = VectorXd::Zero(mesh.cellCount());
    for (const auto& cell : mesh.cells()) {
        const auto& center = cell.center();
        if (center.norm() <= 0.15) {
            source_field_data[cell.id()] = 100000.0;
        }
    }

    // diffusion coefficient
    auto kappa = Tensor("kappa", mesh, Matrix3d::Identity() * 1e-5);

    auto eqn = eqn::Transport(
        diffusion::Corrected<Tensor, nonortho::OverRelaxedCorrector, Scalar>(kappa, T));

    // solve
    auto solver = solver::BiCGSTAB<Scalar>();
    solver.solve(eqn, 20, 1e-20, 1);

    prism::export_field_vtu(eqn.field(), "solution.vtu");

    auto gradT_x = ops::grad(T, Coord::X);
    prism::export_field_vtu(gradT_x, "gradT_x_ls.vtu");

    // this is just for diagnostics, using torus_2d test case
    auto foam_temp = jsonToPrismField(fs::path(unv_file_name).parent_path() / "foam.json", mesh);
    VectorXd diff_abs = (foam_temp.values().array() - T.values().array()).abs();
    log::info("l2-norm difference between our solution and the one from the FOAM solver: {}",
              diff_abs.norm());

    auto diff_abs_field = field::Scalar("diff", mesh, diff_abs);
    prism::export_field_vtu(diff_abs_field, "diff_abs.vtu");

    return 0;
}
