#include <catch2/catch_test_macros.hpp>

#include "prism/field/scalar.h"
#include "prism/mesh/unv.h"

#include <filesystem>

TEST_CASE("GeneralScalar History with Real Mesh", "[field][history]") {
    // Load mesh from test_advection_1d.cpp
    const auto* mesh_file = "tests/cases/versteeg_advection_1d/mesh.unv";
    auto boundary_file = std::filesystem::path(mesh_file).parent_path() / "fields.json";
    auto mesh = prism::mesh::UnvToPMeshConverter(mesh_file, boundary_file).toPMesh();

    const std::size_t num_cells = mesh->cellCount();

    // Define distinct field values for timesteps
    prism::VectorXd t0_values = prism::VectorXd::Constant(num_cells, 1.0);
    prism::VectorXd t1_values = prism::VectorXd::Constant(num_cells, 2.0);
    prism::VectorXd t2_values = prism::VectorXd::Constant(num_cells, 3.0);
    prism::VectorXd t3_values = prism::VectorXd::Constant(num_cells, 4.0);
    prism::VectorXd t4_values = prism::VectorXd::Constant(num_cells, 5.0);
    prism::VectorXd t5_values = prism::VectorXd::Constant(num_cells, 6.0);

    prism::field::Scalar T("T", mesh, t0_values);

    SECTION("Basic History Tracking and Updates") {
        T.setHistorySize(3); // Keep track of t-1, t-2, t-3

        // Initial state: T = initial_values
        REQUIRE(T.values().isApprox(t0_values));
        REQUIRE_FALSE(T.prevValues().has_value());

        // Update 1: T = t1_values
        T.update(t1_values);
        REQUIRE(T.values().isApprox(t1_values));
        REQUIRE(T.prevValues().has_value());
        REQUIRE(T.prevValues()->isApprox(t0_values));
        REQUIRE_FALSE(T.prevPrevValues().has_value());

        // Update 2: T = t2_values
        T.update(t2_values);
        REQUIRE(T.values().isApprox(t2_values));
        REQUIRE(T.prevValues().has_value());
        REQUIRE(T.prevValues()->isApprox(t1_values));
        REQUIRE(T.prevPrevValues().has_value());
        REQUIRE(T.prevPrevValues()->isApprox(t0_values));
        REQUIRE_FALSE(T.getHistory(3).has_value()); // Index 3 is out of bounds

        // Update 3: T = t3_values (history is now full)
        T.update(t3_values);
        REQUIRE(T.values().isApprox(t3_values));
        REQUIRE(T.prevValues()->isApprox(t2_values));
        REQUIRE(T.prevPrevValues()->isApprox(t1_values));
        REQUIRE(T.getHistory(2)->isApprox(t0_values));
        REQUIRE_FALSE(T.getHistory(3).has_value());

        // Update 4: T = t4_values (oldest history should be dropped)
        T.update(t4_values);
        REQUIRE(T.values().isApprox(t4_values));
        REQUIRE(T.prevValues()->isApprox(t3_values));
        REQUIRE(T.prevPrevValues()->isApprox(t2_values));
        REQUIRE(T.getHistory(2)->isApprox(t1_values));
        REQUIRE_FALSE(T.getHistory(3).has_value());
    }

    SECTION("Resizing History - Truncation") {
        T.setHistorySize(4);
        T.update(t1_values);
        T.update(t2_values);
        T.update(t3_values);
        T.update(t4_values);

        // Current: t4_values. History: [t3_values, t2_values, t1_values, initial_values]
        REQUIRE(T.getHistory(0)->isApprox(t3_values));
        REQUIRE(T.getHistory(1)->isApprox(t2_values));
        REQUIRE(T.getHistory(2)->isApprox(t1_values));
        REQUIRE(T.getHistory(3)->isApprox(t0_values));

        // Resize to smaller (2), should truncate oldest history
        T.setHistorySize(2);
        REQUIRE(T.getHistory(0)->isApprox(t3_values));
        REQUIRE(T.getHistory(1)->isApprox(t2_values));
        REQUIRE_FALSE(T.getHistory(2).has_value());

        // Update to see if it respects new size
        T.update(t5_values);
        REQUIRE(T.values().isApprox(t5_values));
        REQUIRE(T.prevValues()->isApprox(t4_values));
        REQUIRE(T.prevPrevValues()->isApprox(t3_values));
        REQUIRE_FALSE(T.getHistory(2).has_value());
    }

    SECTION("Resizing History - Expansion") {
        T.setHistorySize(1);
        T.update(t1_values);
        T.update(t2_values);

        // Current: t2_values. History: [t1_values]
        REQUIRE(T.prevValues()->isApprox(t1_values));
        REQUIRE_FALSE(T.prevPrevValues().has_value());

        // Resize to larger (3), should preserve existing and allow more history
        T.setHistorySize(3);
        REQUIRE(T.prevValues()->isApprox(t1_values)); // Still t1_values
        REQUIRE_FALSE(T.prevPrevValues().has_value()); // Still no t-2

        // Update to fill new capacity
        T.update(t3_values);
        REQUIRE(T.values().isApprox(t3_values));
        REQUIRE(T.prevValues()->isApprox(t2_values));
        REQUIRE(T.prevPrevValues()->isApprox(t1_values));
        REQUIRE_FALSE(T.getHistory(3).has_value());

        T.update(t4_values);
        REQUIRE(T.values().isApprox(t4_values));
        REQUIRE(T.prevValues()->isApprox(t3_values));
        REQUIRE(T.prevPrevValues()->isApprox(t2_values));
        REQUIRE(T.getHistory(2)->isApprox(t1_values));
    }

    SECTION("Clearing History") {
        T.setHistorySize(2);
        T.update(t1_values);
        T.update(t2_values);
        REQUIRE(T.prevValues().has_value());

        // Set size to 0 to clear
        T.setHistorySize(0);
        REQUIRE_FALSE(T.prevValues().has_value());
        REQUIRE_FALSE(T.prevPrevValues().has_value());

        // Updating should not populate history
        T.update(t3_values);
        REQUIRE_FALSE(T.prevValues().has_value());
    }
}
