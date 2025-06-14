#include "rhie_chow.h"

namespace prism::ops::detail {

auto pressureGradCalculated(const mesh::PMesh& mesh,
                            const mesh::Face& face,
                            const field::Pressure& P,
                            const Vector3d& gradp_avg) -> Vector3d {
    // This function is based on Eqn (15.62)
    const auto& owner = mesh.cell(face.owner());
    const auto& neighbor = mesh.cell(face.neighbor().value());
    auto p_owner = P.valueAtCell(owner);
    auto p_neigh = P.valueAtCell(neighbor);

    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();
    Vector3d e_CF = d_CF / d_CF_norm;

    auto gradp_ortho = (p_neigh - p_owner) / d_CF_norm;

    // return gradp_avg + ((gradp_ortho - gradp_avg.dot(e_CF)) * e_CF); // book
    return gradp_ortho * e_CF; // aidan
    // return gradp_avg - (gradp_avg.dot(e_CF) * e_CF) + (gradp_ortho * e_CF); // foam
}

} // namespace prism::ops::detail
