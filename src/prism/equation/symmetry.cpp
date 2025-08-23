#include <utility>

#include "boundary.h"
#include "prism/equation/transport.h"
#include "prism/field/scalar.h"
#include "prism/field/velocity.h"
#include "prism/types.h"

namespace prism::eqn::boundary {
/// TODO: many logic here is copied from the NoSlip boundary condition, we should refactor it to
/// avoid code duplication.

auto contribution(Coord coord,
                  const Vector3d& Uc,
                  const Vector3d& n) -> std::pair<double, double> {
    // n
    double nx = n.x();
    double ny = n.y();
    double nz = n.z();

    // Uc
    double uc = Uc.x();
    double vc = Uc.y();
    double wc = Uc.z();

    switch (coord) {
        case Coord::X: {
            // Eqn (15.154)
            double ac = nx * nx;
            double b = ((vc * ny) + (wc * nz)) * nx;
            return {ac, b};
        }
        case Coord::Y: {
            // Eqn (15.1545)
            double ac = ny * ny;
            double b = ((uc * nx) + (wc * nz)) * ny;
            return {ac, b};
        }
        case Coord::Z: {
            // Eqn (15.156)
            double ac = nz * nz;
            double b = ((uc * nx) + (vc * ny)) * nz;
            return {ac, b};
        }
        default: {
            return {0.0, 0.0};
        }
    }
}

void Symmetry<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {
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

        Vector3d n = face.normal();
        Vector3d d_Cb = face.center() - owner.center();
        double d_normal = d_Cb.dot(n);

        Vector3d Uc = U->valueAtCell(owner);

        double g = 2 * mu->valueAtFace(face) * face.area() / d_normal;
        auto [ac, b] = contribution(field->coord().value(), Uc, n);

        sys.insert(owner_id, owner_id, g * ac);
        sys.rhs(owner_id) += -g * b;
    }

    sys.collect();
    eqn.matrix() += sys.matrix();
    eqn.rhs() += sys.rhs();
}

} // namespace prism::eqn::boundary
