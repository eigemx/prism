#pragma once

#include "boundary.h"
#include "prism/operations/operations.h"

namespace prism::scheme::convection {
// forward declarations
template <typename Field>
class IConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <typename F>
class Fixed<convection::IConvection<F>>
    : public ISchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename F>
class VelocityInlet<convection::IConvection<F>>
    : public ISchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <typename F>
class NoSlip<convection::IConvection<F>>
    : public ISchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <typename F>
class Symmetry<convection::IConvection<F>>
    : public ISchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <typename F>
class ZeroGradient<convection::IConvection<F>>
    : public ISchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "zero-gradient"; }
};

template <typename F>
class Outlet<convection::IConvection<F>>
    : public ISchemeBoundaryHandler<convection::IConvection<F>> {
  public:
    void apply(convection::IConvection<F>& scheme, const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};

template <typename F>
void Fixed<convection::IConvection<F>>::apply(convection::IConvection<F>& scheme,
                                              const mesh::BoundaryPatch& patch) {
    const auto phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const double phi_wall = patch.getScalarBoundaryCondition(phi.name());

        const Vector3d& S_f = face.areaVector();
        const Vector3d U_f = scheme.U().valueAtFace(face);

        // TODO: check if this is correct
        // TODO: warn user if mass flow rate is not entering the patch (negative flow rate)
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

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

template <typename F>
void NoSlip<convection::IConvection<F>>::apply(convection::IConvection<F>& scheme,
                                               const mesh::BoundaryPatch& patch) {
    Fixed<convection::IConvection<F>> fixed;
    return fixed.apply(scheme, patch);
}

template <typename F>
void Outlet<convection::IConvection<F>>::apply(convection::IConvection<F>& scheme,
                                               const mesh::BoundaryPatch& patch) {
    _n_reverse_flow_faces = 0;

    const auto phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh.face(face_id);
        const mesh::Cell& owner = mesh.cell(face.owner());
        const std::size_t cell_id = owner.id();

        // face area vector
        const Vector3d& S_f = face.areaVector();

        // use owner cell velocity as the velocity at the outlet face centroid
        const Vector3d U_f = scheme.U().valueAtFace(face);
        const double rho_f = scheme.rho().valueAtCell(owner);
        const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

        if (m_dot_f < 0.0) {
            _n_reverse_flow_faces++;
        }
        // TODO: this assumes an upwind based scheme, this is wrong for central schemes
        // and should be generalized to work for all schemes.

        // Because this is an outlet, owner `cell` is an upstream cell
        scheme.insert(cell_id, cell_id, m_dot_f);
    }

    if (_n_reverse_flow_faces > 0) {
        log::warn("Reverse flow detected in {} faces in outlet flow patch '{}'",
                  _n_reverse_flow_faces,
                  patch.name());
    }
}
} // namespace prism::scheme::boundary
