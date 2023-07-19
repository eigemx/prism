#pragma once

#include <memory>
#include <type_traits>

#include "../field.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::gradient {

// Base class for gradient schemes for explicity calculating the cell gradient of a scalar field.
// All gradient schemes should inherit from this class and define gradient() function.
class GradientSchemeBase {
  public:
    virtual auto gradient_at_cell(const mesh::Cell& c) -> Vector3d = 0;
    virtual auto gradient_at_face(const mesh::Face& f) -> Vector3d = 0;
    virtual auto gradient_field() -> VectorField = 0;

    static auto gradient_at_boundary_face(const mesh::Face& f, const ScalarField& field)
        -> Vector3d;
};

class GreenGauss : public GradientSchemeBase {
  public:
    GreenGauss(const ScalarField& field)
        : _field(field), _cell_gradients(MatrixX3d::Zero(field.mesh().n_cells(), 3)) {}

    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;
    auto gradient_at_face(const mesh::Face& face) -> Vector3d override;
    auto gradient_field() -> VectorField override;

  private:
    const ScalarField& _field;
    MatrixX3d _cell_gradients;
};

class LeastSquares : public GradientSchemeBase {
  public:
    LeastSquares(const ScalarField& field);
    auto gradient_at_cell(const mesh::Cell& cell) -> Vector3d override;
    auto gradient_at_face(const mesh::Face& face) -> Vector3d override;
    auto gradient_field() -> VectorField override;

  private:
    void set_lsq_matrices();
    auto boundary_face_phi(const mesh::Face& face) -> double;

    const ScalarField& _field;
    MatrixX3d _cell_gradients;

    // least-squares optimization matrices
    std::vector<MatrixX3d> _distance_matrices;
    std::vector<MatrixX3d> _pinv_matrices; // pseudo-inverse matrices
};

template <typename G>
auto create(const ScalarField& field)
    -> std::enable_if_t<std::is_base_of_v<GradientSchemeBase, G>, std::shared_ptr<G>> {
    return std::make_shared<G>(field);
}

auto inline GradientSchemeBase::gradient_at_boundary_face(const mesh::Face& face,
                                                          const ScalarField& field) -> Vector3d {
    const auto& mesh = field.mesh();

    auto face_boundary_patch_id = face.boundary_patch_id().value();
    const auto& boundary_patch = mesh.boundary_patches()[face_boundary_patch_id];
    const auto& boundary_condition = boundary_patch.get_bc(field.name());
    auto patch_type = boundary_condition.patch_type();

    switch (patch_type) {
        case mesh::BoundaryPatchType::Empty:
        case mesh::BoundaryPatchType::Outlet:
        case mesh::BoundaryPatchType::Symmetry: {
            return Vector3d {0., 0., 0.};
        }

        case mesh::BoundaryPatchType::Fixed:
        case mesh::BoundaryPatchType::Inlet: {
            const auto& phi_name = field.name();
            auto phi = boundary_patch.get_scalar_bc(phi_name);
            return phi * face.area_vector();
        }

        case mesh::BoundaryPatchType::FixedGradient: {
            const auto& phi_name = field.name();
            auto flux = boundary_patch.get_scalar_bc(phi_name + "-flux");
            return flux * face.area_vector();
        }

        default:
            throw std::runtime_error(
                fmt::format("gradient/gradient.cpp gradient_at_boundary_face(): "
                            "Non-implemented boundary type for boundary patch: '{}'",
                            boundary_patch.name()));
    }
}

} // namespace prism::gradient