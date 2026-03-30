#include <utility>

#include "boundary.h"
#include "prism/equation/transport.h"
#include "prism/types.h"

namespace prism::eqn::boundary {
/// TODO: many logic here is copied from the NoSlip boundary condition, we should refactor it to
/// avoid code duplication.

namespace {
auto contribution(VectorCoord coord, const Vector3d& Uc, const Vector3d& n) -> std::pair<f64, f64> {
    // n
    f64 nx = n.x();
    f64 ny = n.y();
    f64 nz = n.z();

    // Uc
    f64 uc = Uc.x();
    f64 vc = Uc.y();
    f64 wc = Uc.z();

    switch (coord) {
        case VectorCoord::X: {
            // Eqn (15.154)
            f64 ac = nx * nx;
            f64 b = ((vc * ny) + (wc * nz)) * nx;
            return {ac, b};
        }
        case VectorCoord::Y: {
            // Eqn (15.1545)
            f64 ac = ny * ny;
            f64 b = ((uc * nx) + (wc * nz)) * ny;
            return {ac, b};
        }
        case VectorCoord::Z: {
            // Eqn (15.156)
            f64 ac = nz * nz;
            f64 b = ((uc * nx) + (vc * ny)) * nz;
            return {ac, b};
        }
        default: {
            return {0.0, 0.0};
        }
    }
}
} // namespace

void Symmetry<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {
    auto field = eqn.field();
    const auto& mesh = eqn.field()->mesh();

    auto conv_scheme = eqn.convectionScheme();
    const auto& U = conv_scheme->U();

    const auto& diff_scheme = eqn.diffusionScheme();
    const auto& mu = diff_scheme->kappa();

    LinearSystem sys(mesh->cellCount());

    for (size_t face_id : patch.facesIds()) {
        const auto& face = mesh->face(face_id);
        const auto& owner = mesh->cell(face.owner());
        const auto owner_id = owner.id();

        Vector3d n = face.normal();
        Vector3d d_Cb = face.center() - owner.center();
        f64 d_normal = d_Cb.dot(n);

        Vector3d Uc = U->valueAtCell(owner);

        Vector3d Sf = face.areaVector();
        f64 g = 2 * mu->multiply(Sf, face).norm() / d_normal;
        auto [ac, b] = contribution(field->coord().value(), Uc, n);

        sys.insert(owner_id, owner_id, g * ac);
        sys.rhs(owner_id) += -g * b;
    }

    sys.collect();
    eqn.matrix() += sys.matrix();
    eqn.rhs() += sys.rhs();
}

} // namespace prism::eqn::boundary
