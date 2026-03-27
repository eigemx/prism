#include <prism/prism.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <map>
#include <stdexcept>
#include <string>

using namespace prism;
using namespace prism::mesh;

TEST_CASE("FieldInfo basic construction", "[boundary-files][FieldInfo]") {
    FieldInfo info("T", "scalar");
    REQUIRE(info.name() == "T");
    REQUIRE(info.type() == "scalar");
    REQUIRE(info.gradScheme() == std::nullopt);
    REQUIRE(info.units() == std::nullopt);
}

TEST_CASE("FieldInfo with gradScheme", "[boundary-files][FieldInfo]") {
    FieldInfo info("U", "velocity", "greenGauss", std::nullopt);
    REQUIRE(info.name() == "U");
    REQUIRE(info.type() == "velocity");
    REQUIRE(info.gradScheme().has_value());
    REQUIRE(*info.gradScheme() == "greenGauss");
    REQUIRE(info.units() == std::nullopt);
}

TEST_CASE("FieldInfo with units", "[boundary-files][FieldInfo]") {
    std::vector<double> units_vec = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    FieldInfo info("P", "scalar", std::nullopt, units_vec);
    REQUIRE(info.name() == "P");
    REQUIRE(info.type() == "scalar");
    REQUIRE(info.gradScheme() == std::nullopt);
    REQUIRE(info.units().has_value());
    REQUIRE(info.units()->size() == 9);
}

TEST_CASE("FieldInfo with both optionals", "[boundary-files][FieldInfo]") {
    std::vector<double> units_vec = {1.0, 2.0, 3.0};
    FieldInfo info("T", "temperature", "leastSquares", units_vec);
    REQUIRE(info.name() == "T");
    REQUIRE(info.type() == "temperature");
    REQUIRE(info.gradScheme().has_value());
    REQUIRE(*info.gradScheme() == "leastSquares");
    REQUIRE(info.units().has_value());
    REQUIRE(info.units()->size() == 3);
}

TEST_CASE("BoundaryCondition Nil kind", "[boundary-files][BoundaryCondition]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Nil, 0.0, "empty");
    REQUIRE(bc.valueKind() == BoundaryConditionValueKind::Nil);
    REQUIRE(bc.kindString() == "empty");
}

TEST_CASE("BoundaryCondition Scalar kind", "[boundary-files][BoundaryCondition]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Scalar, 1.5, "fixedValue");
    REQUIRE(bc.valueKind() == BoundaryConditionValueKind::Scalar);
    REQUIRE(bc.kindString() == "fixedValue");
    REQUIRE(std::get<double>(bc.value()) == 1.5);
}

TEST_CASE("BoundaryCondition Vector kind", "[boundary-files][BoundaryCondition]") {
    Vector3d vec {1.0, 2.0, 3.0};
    BoundaryCondition bc(BoundaryConditionValueKind::Vector, vec, "fixedValue");
    REQUIRE(bc.valueKind() == BoundaryConditionValueKind::Vector);
    REQUIRE(bc.kindString() == "fixedValue");
    Vector3d extracted = std::get<Vector3d>(bc.value());
    REQUIRE(extracted[0] == 1.0);
    REQUIRE(extracted[1] == 2.0);
    REQUIRE(extracted[2] == 3.0);
}

TEST_CASE("bad_Patch empty map throws", "[boundary-files][BoundaryPatch]") {
    std::map<std::string, BoundaryCondition> empty_map;
    REQUIRE_THROWS_AS(BoundaryPatch("Wall", empty_map), std::runtime_error);
}

TEST_CASE("Patch with empty bc type is empty", "[boundary-files][BoundaryPatch]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Nil, 0.0, "empty");
    std::map<std::string, BoundaryCondition> m = {{"T", bc}};
    BoundaryPatch patch("EmptyPatch", m);
    REQUIRE(patch.isEmpty() == true);
    REQUIRE(patch.name() == "EmptyPatch");
}

TEST_CASE("Patch with non-empty bc type is not empty", "[boundary-files][BoundaryPatch]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Scalar, 300.0, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"T", bc}};
    BoundaryPatch patch("Wall", m);
    REQUIRE(patch.isEmpty() == false);
}

TEST_CASE("getBC direct lookup", "[boundary-files][BoundaryPatch]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Scalar, 300.0, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"T", bc}};
    BoundaryPatch patch("Wall", m);
    const auto& retrieved = patch.getBoundaryCondition("T");
    REQUIRE(retrieved.kindString() == "fixed");
}

TEST_CASE("getBC component field fallback", "[boundary-files][BoundaryPatch]") {
    Vector3d vec {1.0, 2.0, 3.0};
    BoundaryCondition bc(BoundaryConditionValueKind::Vector, vec, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"U", bc}};
    BoundaryPatch patch("Inlet", m);
    const auto& retrieved = patch.getBoundaryCondition("U_x");
    REQUIRE(retrieved.kindString() == "fixed");
}

TEST_CASE("bad_getBC unknown field throws", "[boundary-files][BoundaryPatch]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Scalar, 300.0, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"T", bc}};
    BoundaryPatch patch("Wall", m);
    REQUIRE_THROWS_AS(patch.getBoundaryCondition("U"), std::runtime_error);
}

TEST_CASE("getScalarBC scalar field", "[boundary-files][BoundaryPatch]") {
    BoundaryCondition bc(BoundaryConditionValueKind::Scalar, 350.0, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"T", bc}};
    BoundaryPatch patch("Inlet", m);
    REQUIRE(patch.getScalarBoundaryCondition("T") == 350.0);
}

TEST_CASE("getScalarBC component from vector", "[boundary-files][BoundaryPatch]") {
    Vector3d vec {1.0, 2.0, 3.0};
    BoundaryCondition bc(BoundaryConditionValueKind::Vector, vec, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"U", bc}};
    BoundaryPatch patch("Inlet", m);
    REQUIRE(patch.getScalarBoundaryCondition("U_x") == 1.0);
    REQUIRE(patch.getScalarBoundaryCondition("U_y") == 2.0);
    REQUIRE(patch.getScalarBoundaryCondition("U_z") == 3.0);
}

TEST_CASE("bad_getScalarBC not scalar throws", "[boundary-files][BoundaryPatch]") {
    Vector3d vec {1.0, 2.0, 3.0};
    BoundaryCondition bc(BoundaryConditionValueKind::Vector, vec, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"U", bc}};
    BoundaryPatch patch("Inlet", m);
    REQUIRE_THROWS_AS(patch.getScalarBoundaryCondition("U"), std::runtime_error);
}

TEST_CASE("bad_getScalarBCSubfield invalid component throws", "[boundary-files][BoundaryPatch]") {
    Vector3d vec {1.0, 2.0, 3.0};
    BoundaryCondition bc(BoundaryConditionValueKind::Vector, vec, "fixed");
    std::map<std::string, BoundaryCondition> m = {{"U", bc}};
    BoundaryPatch patch("Inlet", m);
    REQUIRE_THROWS_AS(patch.getScalarBoundaryCondition("U_w"), std::runtime_error);
}

TEST_CASE("good_one_field loads correctly", "[boundary-files][MeshBoundary]") {
    auto fields_file = std::filesystem::path("tests/cases/boundary/good_one_field/fields.json");
    MeshBoundary mb(fields_file);
    const auto& fields = mb.fields();
    REQUIRE(fields.size() == 1);
    REQUIRE(fields[0].name() == "T");
    REQUIRE(fields[0].type() == "scalar");
    REQUIRE(mb.patches().size() == 3);
}

TEST_CASE("good_two_fields parses scalar and vector", "[boundary-files][MeshBoundary]") {
    auto fields_file = std::filesystem::path("tests/cases/boundary/good_two_fields/fields.json");
    MeshBoundary mb(fields_file);
    const auto& fields = mb.fields();
    REQUIRE(fields.size() == 2);
    REQUIRE(fields[0].name() == "T");
    REQUIRE(fields[0].gradScheme().has_value());
    REQUIRE(*fields[0].gradScheme() == "greenGauss");
    REQUIRE(fields[0].units().has_value());
    REQUIRE(fields[1].name() == "U");
    REQUIRE(*fields[1].gradScheme() == "leastSquares");
    REQUIRE(mb.patches().size() == 3);
    const auto& patch_inlet = *std::find_if(mb.patches().begin(),
                                            mb.patches().end(),
                                            [](const auto& p) { return p.name() == "Inlet"; });
    REQUIRE(patch_inlet.getScalarBoundaryCondition("T") == 350.0);
    Vector3d u_inlet = patch_inlet.getVectorBoundaryCondition("U");
    REQUIRE(u_inlet[0] == 1.0);
    REQUIRE(u_inlet[1] == 0.0);
    REQUIRE(u_inlet[2] == 0.0);
}

TEST_CASE("good_empty_fields returns empty", "[boundary-files][MeshBoundary]") {
    auto fields_file =
        std::filesystem::path("tests/cases/boundary/good_empty_fields/fields.json");
    MeshBoundary mb(fields_file);
    REQUIRE(mb.fields().empty());
    REQUIRE(mb.patches().empty());
}

TEST_CASE("bad_missing_name throws", "[boundary-files][MeshBoundary]") {
    auto fields_file = std::filesystem::path("tests/cases/boundary/bad_missing_name/fields.json");
    REQUIRE_THROWS_AS(MeshBoundary(fields_file), std::runtime_error);
}

TEST_CASE("bad_missing_type throws", "[boundary-files][MeshBoundary]") {
    auto fields_file = std::filesystem::path("tests/cases/boundary/bad_missing_type/fields.json");
    REQUIRE_THROWS_AS(MeshBoundary(fields_file), std::runtime_error);
}

TEST_CASE("bad_missing_field_file throws", "[boundary-files][MeshBoundary]") {
    auto fields_file =
        std::filesystem::path("tests/cases/boundary/bad_missing_field_file/fields.json");
    REQUIRE_THROWS_AS(MeshBoundary(fields_file), std::runtime_error);
}

TEST_CASE("bad_inconsistent_patches throws", "[boundary-files][MeshBoundary]") {
    auto fields_file =
        std::filesystem::path("tests/cases/boundary/bad_inconsistent_patches/fields.json");
    REQUIRE_THROWS_AS(MeshBoundary(fields_file), std::runtime_error);
}