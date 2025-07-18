#pragma once

#include <cmath>

#include "boundary.h"
#include "convection_boundary.h"
#include "convection_hr.h"
#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/mesh/cell.h"
#include "scheme.h"

namespace prism::scheme::convection {

namespace detail {
// coefficients for the discretized convection equation
struct CoeffsTriplet {
    double ownerCoeff {0.0};    // cell
    double neighborCoeff {0.0}; // neighbor
    double rhs {0.0};           // source
};
} // namespace detail

// Basic base class for all convection schemes, without templating clutter.
class IConvection {};

// Finite volume scheme for the discretization of the convection term
template <field::IVectorBased ConvectiveField, typename Field>
class IAppliedConvection
    : public IConvection,
      public IFullScheme<Field>,
      public prism::boundary::BHManagerProvider<
          boundary::ISchemeBoundaryHandler<IAppliedConvection<ConvectiveField, Field>>> {
  public:
    IAppliedConvection(ConvectiveField U, Field phi);

    auto needsCorrection() const noexcept -> bool override { return true; }
    auto inline field() -> Field override { return _phi; }
    auto inline U() -> const ConvectiveField& { return _U; }

    using ConvectiveFieldType = ConvectiveField;
    using FieldType = Field;

  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;
    ConvectiveField _U;
    Field _phi;

    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;
};

// Concept for diffusion schemes that are based on IAppliedConvection.
template <typename T>
concept IAppliedConvectionBased =
    std::derived_from<T,
                      IAppliedConvection<typename T::ConvectiveFieldType, typename T::FieldType>>;

// Upwind scheme
template <field::IVectorBased ConvectiveField, typename Field>
class Upwind : public IAppliedConvection<ConvectiveField, Field> {
  public:
    Upwind(ConvectiveField U, Field phi) : IAppliedConvection<ConvectiveField, Field>(U, phi) {}

  protected:
    struct WeightsNVF {
        double l {0.0};
        double k {0.0};
    };

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
    virtual auto weightsNVF(double phi_tilde) const noexcept -> WeightsNVF;
};

// Central difference scheme
template <field::IVectorBased ConvectiveField, typename Field>
class CentralDifference : public Upwind<ConvectiveField, Field> {
  public:
    CentralDifference(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

// Second order upwind (linear upwind) scheme
template <field::IVectorBased ConvectiveField, typename Field>
class LinearUpwind : public Upwind<ConvectiveField, Field> {
  public:
    LinearUpwind(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

// QUICK scheme
template <field::IVectorBased ConvectiveField, typename Field>
class QUICK : public Upwind<ConvectiveField, Field> {
  public:
    QUICK(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

// FROMM scheme
template <field::IVectorBased ConvectiveField, typename Field>
class FROMM : public Upwind<ConvectiveField, Field> {
  public:
    FROMM(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

// MINMOD scheme
template <field::IVectorBased ConvectiveField, typename Field>
class MINMOD : public Upwind<ConvectiveField, Field> {
  public:
    MINMOD(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

// MUSCL scheme
template <field::IVectorBased ConvectiveField, typename Field>
class MUSCL : public Upwind<ConvectiveField, Field> {
  public:
    MUSCL(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

// SMART scheme
template <field::IVectorBased ConvectiveField, typename Field>
class SMART : public Upwind<ConvectiveField, Field> {
  public:
    SMART(ConvectiveField U, Field phi) : Upwind<ConvectiveField, Field>(U, phi) {}

  protected:
    auto weightsNVF(double phi_tilde) const noexcept
        -> Upwind<ConvectiveField, Field>::WeightsNVF override;
};

template <field::IVectorBased ConvectiveField, typename Field>
IAppliedConvection<ConvectiveField, Field>::IAppliedConvection(ConvectiveField U, Field phi)
    /// TODO: check why _phi.mesh()->cellCount() is not working, as phi should be in a moved
    /// state. Also, we need to avoid std::move and just make the constructor take a const
    /// reference.
    : _U(std::move(U)), _phi(std::move(phi)), IFullScheme<Field>(phi.mesh()->cellCount()) {
    // add default boundary handlers for IConvection based types
    using Scheme = std::remove_reference_t<decltype(*this)>;
    this->boundaryHandlersManager().template addHandler<boundary::Fixed<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Outlet<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::Symmetry<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::ZeroGradient<Scheme>>();
    this->boundaryHandlersManager().template addHandler<boundary::NoSlip<Scheme>>();
}

template <field::IVectorBased ConvectiveField, typename Field>
void IAppliedConvection<ConvectiveField, Field>::applyInterior(const mesh::Face& face) {
    const auto& mesh = _phi.mesh();
    const mesh::Cell& owner = mesh->cell(face.owner());
    const mesh::Cell& neighbor = mesh->cell(face.neighbor().value());
    const std::size_t owner_id = owner.id();
    const std::size_t neighbor_id = neighbor.id();

    // since face is owned by `owner`, the flux will be in the intended way, as the face normal
    // points from cell center to face center, so no need to flup the sign of the flux.
    const double mdot_f = _U.fluxAtFace(face);

    auto [a_C, a_N, b] = interpolate(mdot_f, owner, neighbor, face);
    auto [x_C, x_N, s] = interpolate(-mdot_f, neighbor, owner, face); // NOLINT

    this->insert(owner_id, owner_id, a_C);
    this->insert(owner_id, neighbor_id, a_N);

    this->insert(neighbor_id, neighbor_id, x_C);
    this->insert(neighbor_id, owner_id, x_N);

    this->rhs(owner_id) += b;
    this->rhs(neighbor_id) += s;
}

template <field::IVectorBased ConvectiveField, typename Field>
void IAppliedConvection<ConvectiveField, Field>::applyBoundary() {
    prism::boundary::detail::applyBoundary("prism::scheme::convection::IAppliedConvection",
                                           *this);
}

template <field::IVectorBased ConvectiveField, typename Field>
auto Upwind<ConvectiveField, Field>::interpolate(double m_dot,
                                                 const mesh::Cell& cell,     // NOLINT
                                                 const mesh::Cell& neighbor, // NOLINT
                                                 const mesh::Face& face)     // NOLINT
    -> detail::CoeffsTriplet {
    // when `cell` is the upwind cell
    double phi_tilde = phiTilde(this->field(), cell, neighbor);
    double phi_upwind = phiAtDummyUpwind(this->field(), cell, neighbor);
    const auto [l_plus, k_plus] = weightsNVF(phi_tilde);
    double a_C = std::max(m_dot, 0.0) * l_plus;
    double a_N = std::max(m_dot, 0.0) * k_plus;
    double b = -std::max(m_dot, 0.0) * (1 - l_plus - k_plus) * phi_upwind;

    // when `neighbor` is the upwind cell
    phi_tilde = phiTilde(this->field(), neighbor, cell);
    phi_upwind = phiAtDummyUpwind(this->field(), neighbor, cell);
    const auto [l_minus, k_minus] = weightsNVF(phi_tilde);
    a_C += -std::max(-m_dot, 0.0) * k_minus;
    a_N += -std::max(-m_dot, 0.0) * l_minus;
    b += std::max(-m_dot, 0.0) * (1 - l_minus - k_minus) * phi_upwind;
    return {a_C, a_N, b};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto Upwind<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept // NOLINT
    -> WeightsNVF {
    return {1.0, 0.0};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto CentralDifference<ConvectiveField, Field>::weightsNVF(
    double phi_tilde) const noexcept // NOLINT
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    return {0.5, 0.5};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto LinearUpwind<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept // NOLINT
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    return {1.5, 0.0};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto QUICK<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept // NOLINT
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    return {0.75, 3.0 / 8.0};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto FROMM<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept // NOLINT
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    return {1.0, 0.25};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto MINMOD<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    if (phi_tilde > 0.0 && phi_tilde < 0.5) {
        return {1.5, 0.0};
    }
    if (phi_tilde >= 0.5 && phi_tilde < 1.0) {
        return {0.5, 0.5};
    }
    return {1.0, 0};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto MUSCL<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    if (phi_tilde > 0.0 && phi_tilde < 0.25) {
        return {2, 0.0};
    }
    if (phi_tilde >= 0.25 && phi_tilde < 0.75) {
        return {1, 0.25};
    }
    if (phi_tilde >= 0.75 && phi_tilde < 1.0) {
        return {0, 1.0};
    }
    return {1.0, 0};
}

template <field::IVectorBased ConvectiveField, typename Field>
auto SMART<ConvectiveField, Field>::weightsNVF(double phi_tilde) const noexcept
    -> Upwind<ConvectiveField, Field>::WeightsNVF {
    if (phi_tilde > 0.0 && phi_tilde < 1 / 6) {
        return {4, 0.0};
    }
    if (phi_tilde >= 1 / 6 && phi_tilde < 5 / 6) {
        return {0.75, 3 / 8};
    }
    if (phi_tilde >= 5 / 6 && phi_tilde < 1.0) {
        return {0, 1.0};
    }
    return {1.0, 0};
}

} // namespace prism::scheme::convection
