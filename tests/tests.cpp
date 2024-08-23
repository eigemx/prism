#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

TEST_CASE("field::UniformScalar as uniform cell and face values", "[UniformScalar]") {
    using namespace prism;

    const auto* mesh_file = "tests/cases/duct/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).to_pmesh();

    auto T = field::UniformScalar("T", mesh, 1.0);

    REQUIRE(T.valueAtCell(0) == 1.0);
    REQUIRE(T.gradAtFace(mesh.face(0)) == Vector3d {0.0, 0.0, 0.0});
    REQUIRE(T.gradAtCell(mesh.cell(0)) == Vector3d {0.0, 0.0, 0.0});
    REQUIRE(T.gradAtCellStored(mesh.cell(0)) == Vector3d {0.0, 0.0, 0.0});
}