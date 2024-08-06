#pragma once

#include "boundary.h"
#include "prism/operations/operations.h"
#include "spdlog/spdlog.h"

namespace prism::scheme::convection {
// forward declarations
template <typename GradScheme>
class IConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <typename G>
class Fixed<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename G>
class Symmetry<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename G>
class Outlet<convection::IConvection<G>>
    : public FVSchemeBoundaryHandler<convection::IConvection<G>> {
  public:
    void apply(convection::IConvection<G>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};

template <typename G>
void Fixed<convection::IConvection<G>>::apply(convection::IConvection<G>& scheme,
                                              const mesh::BoundaryPatch& patch) {
    assert(scheme.field().has_value());

    const auto phi = scheme.field().value();
    const auto& mesh = phi.mesh();

    for (const auto face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const double phi_wall = patch.getScalarBoundaryCondition(phi.name());

        const Vector3d& S_f = face.area_vector();
        const Vector3d U_f = scheme.U().valueAtFace(face);

        // TODO: check if this is correct
        // TODO: warn user if mass flow rate is not entering the patch (negative flow rate)
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::face_mdot(rho_f, U_f, S_f);

        // TODO: this assumes an upwind based scheme, this is wrong for central schemes
        // and should be generalized to work for all schemes.

        // in case owner cell is an upstream cell
        // TODO: face value should be the same regardless of the upstream cell mass flow rate, so,
        // applying the right hand side should not depend on a std::max() call and should be
        // definite.
        const std::size_t cell_id = owner.id();
        scheme.insert(cell_id, cell_id, std::max(m_dot_f, 0.0));
        scheme.rhs(cell_id) += std::max(-m_dot_f * phi_wall, 0.0);
    }
}

template <typename G>
void Outlet<convection::IConvection<G>>::apply(convection::IConvection<G>& scheme,
                                               const mesh::BoundaryPatch& patch) {
    assert(scheme.field().has_value());
    _n_reverse_flow_faces = 0;

    const auto phi = scheme.field().value();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const std::size_t cell_id = owner.id();

        // face area vector
        const Vector3d& S_f = face.area_vector();

        // use owner cell velocity as the velocity at the outlet face centroid
        const Vector3d U_f = scheme.U().valueAtFace(face);
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::face_mdot(rho_f, U_f, S_f);

        if (m_dot_f < 0.0) {
            _n_reverse_flow_faces++;
        }
        // TODO: this assumes an upwind based scheme, this is wrong for central schemes
        // and should be generalized to work for all schemes.

        // Because this is an outlet, owner `cell` is an upstream cell
        scheme.insert(cell_id, cell_id, m_dot_f);
    }

    if (_n_reverse_flow_faces > 0) {
        spdlog::warn("Reverse flow detected in {} faces in outlet flow patch '{}'",
                     _n_reverse_flow_faces,
                     patch.name());
    }
}
} // namespace prism::scheme::boundary
