#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "prism/field/scalar.h"
#include "prism/mesh/unv.h"

TEST_CASE("Scalar updateInteriorFaces", "[field][update]") {
    const auto* mesh_file = "tests/cases/versteeg_advection_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = prism::mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    const std::size_t num_cells = mesh->cellCount();
    const std::size_t num_faces = mesh->faceCount();

    prism::VectorXd cell_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::field::Scalar T("T", mesh, cell_values);

    SECTION("initializes face values if not set") {
        REQUIRE_FALSE(T.hasFaceValues());

        std::size_t interior_face_count = 0;
        for (const auto& face : mesh->interiorFaces()) {
            interior_face_count++;
        }

        T.updateInteriorFaces([]([[maybe_unused]] const prism::mesh::Face& face) { return 2.0; });

        REQUIRE(T.hasFaceValues());
        for (const auto& face : mesh->interiorFaces()) {
            REQUIRE(T.valueAtFace(face) == 2.0);
        }
    }

    SECTION("updates only interior faces") {
        T.updateInteriorFaces([]([[maybe_unused]] const prism::mesh::Face& face) { return 5.0; });

        for (const auto& patch : mesh->boundaryPatches()) {
            if (patch.isEmpty()) {
                continue;
            }
            for (const auto& face_id : patch.facesIds()) {
                INFO("Boundary face " << face_id << " should still be 1.0");
                REQUIRE(T.valueAtFace(face_id) == 1.0);
            }
        }
    }

    SECTION("callback receives correct face IDs") {
        T.updateInteriorFaces(
            [&mesh](const prism::mesh::Face& face) { return static_cast<prism::f64>(face.id()); });

        for (const auto& face : mesh->interiorFaces()) {
            REQUIRE(T.valueAtFace(face) == static_cast<prism::f64>(face.id()));
        }
    }
}


TEST_CASE("Scalar updateFaces", "[field][update]") {
    const auto* mesh_file = "tests/cases/versteeg_advection_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = prism::mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    const std::size_t num_cells = mesh->cellCount();

    prism::VectorXd cell_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::field::Scalar T("T", mesh, cell_values);

    SECTION("updates all faces including boundary") {
        T.updateFaces([]([[maybe_unused]] const prism::mesh::Face& face) { return 10.0; });

        for (std::size_t i = 0; i < mesh->faceCount(); ++i) {
            REQUIRE(T.valueAtFace(i) == 10.0);
        }
    }

    SECTION("callback receives correct face IDs for all faces") {
        T.updateFaces(
            [](const prism::mesh::Face& face) { return static_cast<prism::f64>(face.id() * 2); });

        for (std::size_t i = 0; i < mesh->faceCount(); ++i) {
            REQUIRE(T.valueAtFace(i) == static_cast<prism::f64>(i * 2));
        }
    }
}


TEST_CASE("Scalar updateCells", "[field][update]") {
    const auto* mesh_file = "tests/cases/versteeg_advection_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = prism::mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    const std::size_t num_cells = mesh->cellCount();

    prism::VectorXd cell_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::field::Scalar T("T", mesh, cell_values);

    SECTION("updates all cell values") {
        T.updateCells([]([[maybe_unused]] const prism::mesh::Cell& cell) { return 7.0; });

        for (const auto& cell : mesh->cells()) {
            REQUIRE(T.valueAtCell(cell) == 7.0);
        }
    }

    SECTION("callback receives correct cell IDs") {
        T.updateCells(
            [](const prism::mesh::Cell& cell) { return static_cast<prism::f64>(cell.id() + 1); });

        for (const auto& cell : mesh->cells()) {
            REQUIRE(T.valueAtCell(cell) == static_cast<prism::f64>(cell.id() + 1));
        }
    }

    SECTION("can use lambda with captured state") {
        prism::f64 multiplier = 3.0;
        T.updateCells(
            [multiplier](const prism::mesh::Cell& cell) { return multiplier * cell.id(); });

        for (const auto& cell : mesh->cells()) {
            REQUIRE(T.valueAtCell(cell) == multiplier * cell.id());
        }
    }
}