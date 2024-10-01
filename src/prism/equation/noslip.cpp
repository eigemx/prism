#include <utility>

#include "boundary.h"
#include "prism/field/velocity.h"
#include "prism/types.h"
#include "transport.h"

namespace prism::eqn::boundary {

template <typename Field, typename To>
auto castScheme(const SharedPtr<scheme::IScheme>& ptr) -> SharedPtr<To> {
    return std::dynamic_pointer_cast<To>(ptr);
}

auto contribution(Coord coord,
                  const Vector3d& Uc,
                  const Vector3d& Ub,
                  const Vector3d& n) -> std::pair<double, double> {
    // n
    double nx = n.x();
    double ny = n.y();
    double nz = n.z();

    // Uc
    double uc = Uc.x();
    double vc = Uc.y();
    double wc = Uc.z();

    // Ub
    double ub = Ub.x();
    double vb = Ub.y();
    double wb = Ub.z();

    switch (coord) {
        case Coord::X: {
            // Eqn (15.124)
            double ac = 1 - (nx * nx);
            double b = ub * (1 - (nx * nx));
            b += (vc - vb) * ny * nx;
            b += -(wc - wb) * nz * nx;
            return {ac, b};
        }

        case Coord::Y: {
            // Eqn (15.125)
            double ac = 1 - (ny * ny);
            double b = (uc - ub) * nx * ny;
            b += vb * (1 - (ny * ny));
            b += (wc - wb) * nz * ny;
            return {ac, b};
        }
        case Coord::Z: {
            // Eqn (15.126)
            double ac = 1 - (nz * nz);
            double b = (uc - ub) * nx * nz;
            b += (vc - vb) * ny * nz;
            b += wb * (1 - (nz * nz));
            return {ac, b};
        }
    }
}

void NoSlip<Momentum>::apply(Momentum& eqn, const mesh::BoundaryPatch& patch) {
    field::VelocityComponent field = eqn.field();
    const auto& mesh = eqn.field().mesh();

    using F = field::VelocityComponent;
    using Convection = scheme::convection::IConvection<field::UniformScalar, F>;
    auto conv_scheme = castScheme<F, Convection>(eqn.convectionScheme());
    const auto& U = conv_scheme->U();

    // TODO: this is a bit of a hack, what if kappa is not uniform?
    using Diffusion = scheme::diffusion::IAppliedDiffusion<field::UniformScalar, F>;
    auto diff_scheme = castScheme<F, Diffusion>(eqn.diffusionScheme());
    const auto& mu = diff_scheme->kappa();

    LinearSystem sys(mesh.cellCount());

    for (std::size_t face_id : patch.facesIds()) {
        const auto& face = mesh.face(face_id);
        const auto& owner = mesh.cell(face.owner());
        const auto owner_id = owner.id();

        // we need to calculate the normal distance from the centroid of the boundary element to
        // the face (wall patch)
        Vector3d n = face.normal();
        Vector3d d_CB = face.center() - owner.center();
        double d_normal = d_CB.dot(n);

        Vector3d Uc = U.valueAtCell(owner);
        Vector3d Ub = U.valueAtFace(face);

        double g = mu.valueAtFace(face) * face.area() / d_normal;
        auto [ac, b] = contribution(field.coord().value(), Uc, Ub, n);
        sys.insert(owner_id, owner_id, g * ac);
        sys.rhs(owner_id) += g * b;
    }

    sys.collect();

    eqn.matrix() += sys.matrix();
    eqn.rhs() += sys.rhs();
}

} // namespace prism::eqn::boundary
