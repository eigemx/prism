#pragma once

#include <cmath>

#include "convection_boundary.h"
#include "convection_hr.h"
#include "prism/boundary.h"
#include "prism/field/ifield.h"
#include "prism/field/scalar.h"
#include "prism/mesh/cell.h"
#include "prism/scheme/boundary.h"
#include "prism/scheme/scheme.h"

namespace prism::scheme::convection {

namespace detail {
// coefficients for the discretized convection equation
struct CoeffsTriplet {
    f64 ownerCoeff {0.0};    // cell
    f64 neighborCoeff {0.0}; // neighbor
    f64 rhs {0.0};           // source
};
} // namespace detail

// Basic base class for all convection schemes, without templating clutter.
class IConvection {};

// Finite volume scheme for the discretization of the convection term
class IAppliedConvection : public IConvection,
                           public IFullScheme,
                           public prism::boundary::BHManagerProvider<
                               boundary::ISchemeBoundaryHandler<IAppliedConvection>> {
  public:
    IAppliedConvection(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi);

    auto needsCorrection() const noexcept -> bool override { return true; }
    auto U() -> const SharedPtr<field::IVector>& { return _U; }

  private:
    void applyInterior(const mesh::Face& face) override;
    void applyBoundary() override;

    virtual auto interpolate(f64 m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face) -> detail::CoeffsTriplet = 0;

    SharedPtr<field::IVector> _U;
};

// Concept for diffusion schemes that are based on IAppliedConvection.
template <typename T>
concept IAppliedConvectionBased = std::derived_from<T, IAppliedConvection>;

// Upwind scheme
class Upwind : public IAppliedConvection {
  public:
    Upwind(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : IAppliedConvection(U, phi) {}

  protected:
    struct WeightsNVF {
        f64 l {0.0};
        f64 k {0.0};
    };

    auto interpolate(f64 m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face) -> detail::CoeffsTriplet override;
    virtual auto weightsNVF(f64 phi_tilde) const noexcept -> WeightsNVF;
};

// Central difference scheme
class CentralDifference : public Upwind {
  public:
    CentralDifference(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};

// Second order upwind (linear upwind) scheme
class LinearUpwind : public Upwind {
  public:
    LinearUpwind(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};

// QUICK scheme
class QUICK : public Upwind {
  public:
    QUICK(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};

// FROMM scheme
class FROMM : public Upwind {
  public:
    FROMM(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};

// MINMOD scheme
class MINMOD : public Upwind {
  public:
    MINMOD(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};

// MUSCL scheme
class MUSCL : public Upwind {
  public:
    MUSCL(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};

// SMART scheme
class SMART : public Upwind {
  public:
    SMART(const SharedPtr<field::IVector>& U, const SharedPtr<field::Scalar>& phi)
        : Upwind(U, phi) {}

  protected:
    auto weightsNVF(f64 phi_tilde) const noexcept -> Upwind::WeightsNVF override;
};
} // namespace prism::scheme::convection
