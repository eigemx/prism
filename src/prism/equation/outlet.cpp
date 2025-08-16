#include "boundary.h"
#include "prism/equation/transport.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/log.h"
#include "prism/types.h"

namespace prism::eqn::boundary {
void Outlet<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {
    field::VelocityComponent field = eqn.field();
    const auto& mesh = field.mesh();

    // Momentum's conserved field is always a VelocityComponent
    using F = field::VelocityComponent;
    using Convection = scheme::convection::IAppliedConvection<field::Velocity, F>;
    auto conv_scheme = castScheme<Convection>(eqn.convectionScheme());
    const auto& U = conv_scheme->U();
    auto& phi = eqn.field();

    for (std::size_t face_id : patch.facesIds()) {
        const auto& face = mesh->face(face_id);
        const auto& owner = mesh->cell(face.owner());
        const Vector3d grad_phi_c = phi.gradAtCell(owner);
        const Vector3d e = (face.center() - owner.center()).normalized();
        const Vector3d grad_phi_b = grad_phi_c - grad_phi_c.dot(e) * e;
        Vector3d d_Cb = face.center() - owner.center();

        eqn.rhs(owner.id()) += -U.fluxAtFace(face) * (grad_phi_b.dot(d_Cb));
    }
}
} // namespace prism::eqn::boundary
