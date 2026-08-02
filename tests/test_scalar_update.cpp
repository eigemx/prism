#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "prism/field/scalar.h"
#include "prism/mesh/unv.h"
#include "test_utils.h"

TEST_CASE("Scalar mapInteriorFaceValues", "[field][update]") {
    auto mesh = prism::test::loadMesh("tests/cases/versteeg_advection_1d/mesh.unv");

    const size_t num_cells = mesh->cellCount();
    const size_t num_faces = mesh->faceCount();

    prism::VectorXd cell_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::field::Scalar T("T", mesh, cell_values);

    SECTION("initializes face values if not set") {
        REQUIRE_FALSE(T.hasFaceValues());

        size_t interior_face_count = 0;
        for (const auto& face : mesh->interiorFaces()) {
            interior_face_count++;
        }

        T.mapInteriorFaceValues([]([[maybe_unused]] const prism::mesh::Face& face) { return 2.0; });

        REQUIRE(T.hasFaceValues());
        for (const auto& face : mesh->interiorFaces()) {
            REQUIRE(T.valueAtFace(face) == 2.0);
        }
    }

    SECTION("updates only interior faces") {
        T.mapInteriorFaceValues([]([[maybe_unused]] const prism::mesh::Face& face) { return 5.0; });

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
        T.mapInteriorFaceValues(
            [&mesh](const prism::mesh::Face& face) { return static_cast<prism::f64>(face.id()); });

        for (const auto& face : mesh->interiorFaces()) {
            REQUIRE(T.valueAtFace(face) == static_cast<prism::f64>(face.id()));
        }
    }
}


TEST_CASE("Scalar mapFaceValues", "[field][update]") {
    auto mesh = prism::test::loadMesh("tests/cases/versteeg_advection_1d/mesh.unv");

    const size_t num_cells = mesh->cellCount();

    prism::VectorXd cell_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::field::Scalar T("T", mesh, cell_values);

    SECTION("updates all faces including boundary") {
        T.mapFaceValues([]([[maybe_unused]] const prism::mesh::Face& face) { return 10.0; });

        for (size_t i = 0; i < mesh->faceCount(); ++i) {
            REQUIRE(T.valueAtFace(i) == 10.0);
        }
    }

    SECTION("callback receives correct face IDs for all faces") {
        T.mapFaceValues(
            [](const prism::mesh::Face& face) { return static_cast<prism::f64>(face.id() * 2); });

        for (size_t i = 0; i < mesh->faceCount(); ++i) {
            REQUIRE(T.valueAtFace(i) == static_cast<prism::f64>(i * 2));
        }
    }
}


TEST_CASE("Scalar mapCellValues", "[field][update]") {
    auto mesh = prism::test::loadMesh("tests/cases/versteeg_advection_1d/mesh.unv");

    const size_t num_cells = mesh->cellCount();

    prism::VectorXd cell_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::field::Scalar T("T", mesh, cell_values);

    SECTION("updates all cell values") {
        T.mapCellValues([]([[maybe_unused]] const prism::mesh::Cell& cell) { return 7.0; });

        for (const auto& cell : mesh->cells()) {
            REQUIRE(T.valueAtCell(cell) == 7.0);
        }
    }

    SECTION("callback receives correct cell IDs") {
        T.mapCellValues(
            [](const prism::mesh::Cell& cell) { return static_cast<prism::f64>(cell.id() + 1); });

        for (const auto& cell : mesh->cells()) {
            REQUIRE(T.valueAtCell(cell) == static_cast<prism::f64>(cell.id() + 1));
        }
    }

    SECTION("can use lambda with captured state") {
        prism::f64 multiplier = 3.0;
        T.mapCellValues(
            [multiplier](const prism::mesh::Cell& cell) { return multiplier * cell.id(); });

        for (const auto& cell : mesh->cells()) {
            REQUIRE(T.valueAtCell(cell) == multiplier * cell.id());
        }
    }
}