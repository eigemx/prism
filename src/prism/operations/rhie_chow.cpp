#include "rhie_chow.h"

namespace prism::ops::detail {

auto correctGrad(const mesh::PMesh& mesh,
                 const mesh::Face& face,
                 const field::Pressure& P,
                 const Vector3d& grad_p_f) -> Vector3d {
    // This function is based on Eqn (15.62)

    // TODO: owner and neighbor cells finders, and cell centroid distance calculators are all
    // over the codebase this should be a utility function
    const auto& owner = mesh.cell(face.owner());
    const auto& neighbor = mesh.cell(face.neighbor().value());

    auto P_owner = P.valueAtCell(owner);
    auto P_neigh = P.valueAtCell(neighbor);

    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();
    Vector3d e_CF = d_CF / d_CF_norm;

    auto P_correction = (P_neigh - P_owner) / d_CF_norm;
    P_correction -= grad_p_f.dot(e_CF);

    return P_correction * e_CF;
}


} // namespace prism::ops::detail