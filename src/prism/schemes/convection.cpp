#include "convection.h"

#include "prism/operations/operations.h"

namespace prism::scheme::convection {
IConvection::IConvection(field::Scalar rho, field::Velocity U, field::Scalar phi)
    : _rho(std::move(rho)),
      _U(std::move(U)),
      _phi(std::move(phi)),
      FVScheme(phi.mesh().nCells()) {
    // add default boundary handlers for IConvection based types
    using Scheme = std::remove_reference_t<decltype(*this)>;
    _bc_manager.template add_handler<scheme::boundary::Empty<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Fixed<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Outlet<Scheme>>();
    _bc_manager.template add_handler<scheme::boundary::Symmetry<Scheme>>();
}

void IConvection::apply() {
    apply_boundary();

    for (const auto& iface : _phi.mesh().interiorFaces()) {
        apply_interior(iface);
    }

    collect();
}

void IConvection::apply_interior(const mesh::Face& face) {
    const auto& mesh = _phi.mesh();
    const mesh::Cell& owner = mesh.cell(face.owner());
    const mesh::Cell& neighbor = mesh.cell(face.neighbor().value());

    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    const Vector3d& S_f = mesh::outward_area_vector(face, owner);
    const Vector3d U_f = _U.valueAtFace(face);
    const double rho_f = _rho.valueAtFace(face);
    const double m_dot_f = ops::faceFlowRate(rho_f, U_f, S_f);

    auto [a_C, a_N, b] = interpolate(m_dot_f, owner, neighbor, face);
    auto [x_C, x_N, s] = interpolate(-m_dot_f, neighbor, owner, face); // NOLINT

    insert(owner_id, owner_id, a_C);
    insert(owner_id, neighbor_id, a_N);

    insert(neighbor_id, neighbor_id, x_C);
    insert(neighbor_id, owner_id, x_N);

    rhs(owner_id) += b;
    rhs(neighbor_id) += s;
}

void IConvection::apply_boundary() {
    boundary::detail::apply_boundary("IConvection", *this);
}

} // namespace prism::scheme::convection