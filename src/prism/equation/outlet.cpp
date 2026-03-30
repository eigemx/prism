#include "boundary.h"
#include "prism/equation/transport.h"
#include "prism/types.h"

namespace prism::eqn::boundary {
void Outlet<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {
    auto field = eqn.field();
    const auto& mesh = field->mesh();

    const auto& U = eqn.convectionScheme()->U();
    auto phi = eqn.field();

    for (std::size_t face_id : patch.facesIds()) {
        const auto& face = mesh->face(face_id);
        const auto& owner = mesh->cell(face.owner());
        const Vector3d grad_phi_c = phi->gradAtCell(owner);
        const Vector3d e = (face.center() - owner.center()).normalized();
        const Vector3d grad_phi_b = grad_phi_c - grad_phi_c.dot(e) * e;
        Vector3d d_Cb = face.center() - owner.center();

        eqn.rhs(owner.id()) += -U->fluxAtFace(face) * (grad_phi_b.dot(d_Cb));
    }
}
} // namespace prism::eqn::boundary
