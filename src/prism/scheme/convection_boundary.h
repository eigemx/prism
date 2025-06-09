#pragma once

#include "boundary.h"
#include "prism/operations/operations.h"

namespace prism::scheme::convection {
// forward declarations
template <typename RhoType, typename Field>
class IConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <typename Rho, typename F>
class Fixed<convection::IConvection<Rho, F>>
    : public ISchemeBoundaryHandler<convection::IConvection<Rho, F>> {
  public:
    void apply(convection::IConvection<Rho, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Rho, typename F>
class VelocityInlet<convection::IConvection<Rho, F>>
    : public ISchemeBoundaryHandler<convection::IConvection<Rho, F>> {
  public:
    void apply(convection::IConvection<Rho, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename Rho, typename F>
class NoSlip<convection::IConvection<Rho, F>>
    : public ISchemeBoundaryHandler<convection::IConvection<Rho, F>> {
  public:
    void apply(convection::IConvection<Rho, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <typename Rho, typename F>
class Symmetry<convection::IConvection<Rho, F>>
    : public ISchemeBoundaryHandler<convection::IConvection<Rho, F>> {
  public:
    void apply(convection::IConvection<Rho, F>& scheme,
               const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename Rho, typename F>
class ZeroGradient<convection::IConvection<Rho, F>>
    : public ISchemeBoundaryHandler<convection::IConvection<Rho, F>> {
  public:
    // we treat zero-gradient boundary condition as outlet condition.
    void apply(convection::IConvection<Rho, F>& scheme,
               const mesh::BoundaryPatch& patch) override; 
    auto inline name() const -> std::string override { return "zero-gradient"; }
};

template <typename Rho, typename F>
class Outlet<convection::IConvection<Rho, F>>
    : public ISchemeBoundaryHandler<convection::IConvection<Rho, F>> {
  public:
    void apply(convection::IConvection<Rho, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};

template <typename Rho, typename F>
void Fixed<convection::IConvection<Rho, F>>::apply(convection::IConvection<Rho, F>& scheme,
                                                   const mesh::BoundaryPatch& patch) {
    const auto& phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const double phi_wall = patch.getScalarBoundaryCondition(phi.name());

        const Vector3d& S_f = face.areaVector();
        const Vector3d U_f = scheme.U().valueAtFace(face);
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);
        scheme.rhs(owner.id()) += -m_dot_f * phi_wall;
    }
}

template <typename Rho, typename F>
void NoSlip<convection::IConvection<Rho, F>>::apply(convection::IConvection<Rho, F>& scheme,
                                                    const mesh::BoundaryPatch& patch) {
    Fixed<convection::IConvection<Rho, F>> fixed;
    return fixed.apply(scheme, patch);
}

template <typename Rho, typename F>
void ZeroGradient<convection::IConvection<Rho, F>>::apply(convection::IConvection<Rho, F>& scheme,
                                                    const mesh::BoundaryPatch& patch) {
    Outlet<convection::IConvection<Rho, F>> outlet;
    return outlet.apply(scheme, patch);
}

template <typename Rho, typename F>
void Outlet<convection::IConvection<Rho, F>>::apply(convection::IConvection<Rho, F>& scheme,
                                                    const mesh::BoundaryPatch& patch) {
    _n_reverse_flow_faces = 0;

    const auto& phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const std::size_t owner_id = owner.id();

        // face area vector
        const Vector3d& S_f = face.areaVector();

        // use owner cell velocity as the velocity at the outlet face centroid
        const Vector3d U_f = scheme.U().valueAtFace(face);
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

        if (m_dot_f < 0.0) {
            _n_reverse_flow_faces++;
        }

        scheme.insert(owner_id, owner_id, m_dot_f);
    }

    if (_n_reverse_flow_faces > 0) {
        log::warn("Reverse flow detected in {} faces in outlet flow patch '{}'",
                  _n_reverse_flow_faces,
                  patch.name());
    }
}
} // namespace prism::scheme::boundary
