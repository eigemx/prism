#pragma once

#include "boundary.h"
#include "prism/field/ifield.h"

namespace prism::scheme::convection {
// forward declarations
template <field::IVectorBased ConvectiveField, typename Field>
class IAppliedConvection;

} // namespace prism::scheme::convection

namespace prism::scheme::boundary {
template <field::IVectorBased ConvectiveField, typename F>
class Fixed<convection::IAppliedConvection<ConvectiveField, F>>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection<ConvectiveField, F>> {
  public:
    void apply(convection::IAppliedConvection<ConvectiveField, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "fixed"; }
};

template <field::IVectorBased ConvectiveField, typename F>
class NoSlip<convection::IAppliedConvection<ConvectiveField, F>>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection<ConvectiveField, F>> {
  public:
    void apply(convection::IAppliedConvection<ConvectiveField, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "no-slip"; }
};

template <field::IVectorBased ConvectiveField, typename F>
class Symmetry<convection::IAppliedConvection<ConvectiveField, F>>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection<ConvectiveField, F>> {
  public:
    void apply(convection::IAppliedConvection<ConvectiveField, F>& scheme,
               const mesh::BoundaryPatch& patch) override {}
    auto inline name() const -> std::string override { return "symmetry"; }
};

template <field::IVectorBased ConvectiveField, typename F>
class ZeroGradient<convection::IAppliedConvection<ConvectiveField, F>>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection<ConvectiveField, F>> {
  public:
    // we treat zero-gradient boundary condition as outlet condition.
    void apply(convection::IAppliedConvection<ConvectiveField, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "zero-gradient"; }
};

template <field::IVectorBased ConvectiveField, typename F>
class Outlet<convection::IAppliedConvection<ConvectiveField, F>>
    : public ISchemeBoundaryHandler<convection::IAppliedConvection<ConvectiveField, F>> {
  public:
    void apply(convection::IAppliedConvection<ConvectiveField, F>& scheme,
               const mesh::BoundaryPatch& patch) override;
    auto inline name() const -> std::string override { return "outlet"; }

  private:
    std::size_t _n_reverse_flow_faces {0};
};

template <field::IVectorBased ConvectiveField, typename F>
void Fixed<convection::IAppliedConvection<ConvectiveField, F>>::apply(
    convection::IAppliedConvection<ConvectiveField, F>& scheme,
    const mesh::BoundaryPatch& patch) {
    const auto& phi = scheme.field();
    const auto& mesh = phi.mesh();

    for (const auto face_id : patch.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);
        const mesh::Cell& owner = mesh->cell(face.owner());
        const double phi_wall = patch.getScalarBoundaryCondition(phi.name());
        const auto mdot_f = scheme.U().fluxAtFace(face);

        scheme.insert(owner.id(), owner.id(), std::max(mdot_f, 0.0));
        scheme.rhs(owner.id()) += std::max(-mdot_f, 0.0) * phi_wall;
    }
}

template <field::IVectorBased ConvectiveField, typename F>
void NoSlip<convection::IAppliedConvection<ConvectiveField, F>>::apply(
    convection::IAppliedConvection<ConvectiveField, F>& scheme,
    const mesh::BoundaryPatch& patch) {
    // Check section 12.9.3, in Moukallad et. al, normal velocity at walls is zero, and there is
    // no convection flux.
}

template <field::IVectorBased ConvectiveField, typename F>
void ZeroGradient<convection::IAppliedConvection<ConvectiveField, F>>::apply(
    convection::IAppliedConvection<ConvectiveField, F>& scheme,
    const mesh::BoundaryPatch& patch) {
    Outlet<convection::IAppliedConvection<ConvectiveField, F>> outlet;
    return outlet.apply(scheme, patch);
}

template <field::IVectorBased ConvectiveField, typename F>
void Outlet<convection::IAppliedConvection<ConvectiveField, F>>::apply(
    convection::IAppliedConvection<ConvectiveField, F>& scheme,
    const mesh::BoundaryPatch& patch) {
    _n_reverse_flow_faces = 0;

    const auto& mesh = scheme.field().mesh();

    for (const auto& face_id : patch.facesIds()) {
        const mesh::Face& face = mesh->face(face_id);
        const mesh::Cell& owner = mesh->cell(face.owner());
        const Vector3d& S_f = face.areaVector();

        // use owner cell velocity as the velocity at the outlet face centroid
        const Vector3d U_f = scheme.U().valueAtFace(face);
        const double m_dot_f = U_f.dot(S_f);

        if (m_dot_f < 0.0) {
            _n_reverse_flow_faces++;
        }

        /// TODO: is this correct?
        scheme.insert(owner.id(), owner.id(), m_dot_f);
    }

    if (_n_reverse_flow_faces > 0) {
        log::warn("Reverse flow detected in {} faces in outlet flow patch '{}'",
                  _n_reverse_flow_faces,
                  patch.name());
    }
}
} // namespace prism::scheme::boundary
