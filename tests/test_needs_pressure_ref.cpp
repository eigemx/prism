#include <prism/algorithm/SIMPLE.h>
#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

TEST_CASE("needsPressureReference cavity_hex (all zero-gradient) returns true", "[algo]") {
    using namespace prism;

    auto mesh_file = std::filesystem::path("tests/cases/cavity_hex/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto p = std::make_shared<field::Pressure>("P", mesh, 0.0);
    REQUIRE(algo::needsPressureReference(p));
}

TEST_CASE("needsPressureReference coarsePipeHex (fixed outlet) returns false", "[algo]") {
    using namespace prism;

    auto mesh_file = std::filesystem::path("tests/cases/coarsePipeHexSymmetrical/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto p = std::make_shared<field::Pressure>("P", mesh, 0.0);
    REQUIRE_FALSE(algo::needsPressureReference(p));
}

TEST_CASE("needsPressureReference ductSIMPLE (fixed outlet) returns false", "[algo]") {
    using namespace prism;

    auto mesh_file = std::filesystem::path("tests/cases/ductSIMPLE/mesh.unv");
    auto boundary_file = mesh_file.parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    auto p = std::make_shared<field::Pressure>("P", mesh, 0.0);
    REQUIRE_FALSE(algo::needsPressureReference(p));
}
