#include <utility>

#include "boundary.h"
#include "prism/equation/transport.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/types.h"

namespace prism::eqn::boundary {

auto contribution(Coord coord,
                  const Vector3d& Uc,
                  const Vector3d& Ub,
                  const Vector3d& n) -> std::pair<f64, f64> {
    // n
    f64 nx = n.x();
    f64 ny = n.y();
    f64 nz = n.z();

    // Uc
    f64 uc = Uc.x();
    f64 vc = Uc.y();
    f64 wc = Uc.z();

    // Ub
    f64 ub = Ub.x();
    f64 vb = Ub.y();
    f64 wb = Ub.z();

    switch (coord) {
        case Coord::X: {
            // Eqn (15.124)
            f64 ac = 1 - (nx * nx);
            f64 b = ub * (1 - (nx * nx));
            b += (vc - vb) * ny * nx;
            b += -(wc - wb) * nz * nx;
            return {ac, b};
        }

        case Coord::Y: {
            // Eqn (15.125)
            f64 ac = 1 - (ny * ny);
            f64 b = (uc - ub) * nx * ny;
            b += vb * (1 - (ny * ny));
            b += (wc - wb) * nz * ny;
            return {ac, b};
        }
        case Coord::Z: {
            // Eqn (15.126)
            f64 ac = 1 - (nz * nz);
            f64 b = (uc - ub) * nx * nz;
            b += (vc - vb) * ny * nz;
            b += wb * (1 - (nz * nz));
            return {ac, b};
        }
        default: {
            // Return zeros or throw if invalid Coord
            return {0.0, 0.0};
        }
    }
}

void NoSlip<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {
    auto field = eqn.field();
    const auto& mesh = eqn.field()->mesh();

    // Momentum's conserved field is always a VelocityComponent
    using F = field::VelocityComponent;
    using Convection = scheme::convection::IAppliedConvection;
    auto conv_scheme = castScheme<Convection>(eqn.convectionScheme());
    const auto& U = conv_scheme->U();

    /// TODO: this is a bit of a hack, what if kappa is not scalar?
    using Diffusion = scheme::diffusion::IAppliedDiffusion<field::Scalar>;
    auto diff_scheme = castScheme<Diffusion>(eqn.diffusionScheme());
    const auto& mu = diff_scheme->kappa();

    LinearSystem sys(mesh->cellCount());

    for (size_t face_id : patch.facesIds()) {
        const auto& face = mesh->face(face_id);
        const auto& owner = mesh->cell(face.owner());
        const auto owner_id = owner.id();

        // we need to calculate the normal distance from the centroid of the boundary
        // element to the face (wall patch)
        Vector3d n = face.normal();
        Vector3d d_CB = face.center() - owner.center();
        f64 d_normal = d_CB.dot(n);

        Vector3d Uc = U->valueAtCell(owner);
        Vector3d Ub = U->valueAtFace(face);

        f64 g = mu->valueAtFace(face) * face.area() / d_normal;
        /// TODO: at a test round, ac & b seems to be zero and this apply() function does not do
        /// anything. Check this again.
        auto [ac, b] = contribution(field->coord().value(), Uc, Ub, n);

        sys.insert(owner_id, owner_id, g * ac);
        sys.rhs(owner_id) += g * b;
    }

    sys.collect();

    eqn.matrix() += sys.matrix();
    eqn.rhs() += sys.rhs();
}

} // namespace prism::eqn::boundary
