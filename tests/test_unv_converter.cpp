#include <prism/prism.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <filesystem>

TEST_CASE("test UNV converter owner-neighbor shared face normal direction",
          "[unv-face-normals]") {
    using namespace prism;

    const auto* unv_file_name = "tests/cases/cylinder/mesh.unv";
    auto boundary_file = std::filesystem::path(unv_file_name).parent_path() / "fields.json";
    auto mesh = mesh::UnvToPMeshConverter(unv_file_name, boundary_file).toPMesh();

    std::size_t n_wrong_normals = 0;

    for (const auto& face : mesh.interiorFaces()) {
        const auto& owner = mesh.cell(face.owner());
        const auto& neighbor = mesh.cell(face.neighbor().value());

        // Check if the face normal is pointing from owner to neighbor
        const auto& face_normal = face.normal();

        // vector joining the centers of owner and neighbor cells
        const Vector3d d_CF = neighbor.center() - owner.center();
        const Vector3d n = d_CF.normalized();

        // The dot product should be positive if the face normal is correct
        // (i.e., it points from owner to neighbor)
        // If the dot product is negative, it means the face normal is pointing in the wrong
        // direction (i.e., it points from neighbor to owner)
        if (face_normal.dot(n) < 0) {
            n_wrong_normals++;
        }
    }

    REQUIRE(n_wrong_normals == 0);

    double volume = 0.0;
    for (const auto& cell : mesh.cells()) {
        volume += cell.volume();
    }

    // "cylinder" mesh should have a volume of 31.375
    REQUIRE(std::abs(volume - 31.375) < 1e-4);

    double surface_area = 0.0;
    for (const auto& face : mesh.boundaryFaces()) {
        surface_area += face.area();
    }
    // "cylinder" mesh should have a surface area of 69.0842
    REQUIRE(std::abs(surface_area - 69.0842) < 1e-4);
}
