#include "courant.h"

#include <cmath>

namespace prism::monitor {

auto courantNumber(const field::IVector& mdot, f64 dt) -> field::Scalar {
    const auto& mesh = mdot.mesh();
    VectorXd co = VectorXd::Zero(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        f64 flux_sum = 0.0;
        for (auto face_id : cell.facesIds()) {
            const auto& face = mesh->face(face_id);
            if (face.isBoundary() && mesh->boundaryPatch(face).isEmpty()) {
                continue;
            }
            flux_sum += std::abs(mdot.fluxAtFace(face_id));
        }
        co[cell.id()] = 0.5 * dt * flux_sum / cell.volume();
    }

    return {"Co", mesh, std::move(co)};
}

} // namespace prism::monitor
