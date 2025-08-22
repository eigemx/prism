#include "convection.h"

namespace prism::scheme::boundary {
template <>
void Fixed<convection::IAppliedConvection>::apply(convection::IAppliedConvection& scheme,
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

template <>
void NoSlip<convection::IAppliedConvection>::apply(convection::IAppliedConvection& scheme,
                                                   const mesh::BoundaryPatch& patch) {
    // Check section 12.9.3, in Moukalled et. al, normal velocity at walls is zero, and there is
    // no convection flux.
}

template <>
void ZeroGradient<convection::IAppliedConvection>::apply(convection::IAppliedConvection& scheme,
                                                         const mesh::BoundaryPatch& patch) {
    /// TODO: implement ZeroGradient separate from Outlet
    Outlet<convection::IAppliedConvection> outlet;
    return outlet.apply(scheme, patch);
}

template <>
void Outlet<convection::IAppliedConvection>::apply(convection::IAppliedConvection& scheme,
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
