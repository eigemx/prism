#pragma once

#include "prism/boundary.h"
#include "prism/mesh/boundary.h"

/**
 * @brief Equation-level boundary handlers applied to the transport equation.
 *
 * These handlers modify the discretized equation's matrix/rhs directly, after the
 * finite-volume schemes have been applied by Transport::updateCoeffs(). Only the
 * symmetry handler is currently registered; the equation-level no-slip and outlet
 * handlers were removed:
 *
 * - No-slip: the no-slip wall is owned by the diffusion scheme (diffusion::NoSlip
 *   <Corrected>, a Dirichlet treatment), which is the OpenFOAM-consistent approach
 *   and is validated against OpenFOAM (<1% on the pipe test). The equation-level
 *   Moukalled 15.124 treatment was a redundant alternative: it coincides with the
 *   Dirichlet on the tangential component but gives ac = 0 on the wall-normal
 *   component (uFVM's formulation), diverging from OpenFOAM, so it was not adopted.
 * - Outlet: the outlet momentum condition is already handled at the scheme level
 *   (convection Outlet mass-flux insertion + diffusion Outlet zero-flux + field
 *   zero-gradient), matching both OpenFOAM and uFVM.
 *
 * The symmetry handler remains because it is the one piece the schemes cannot
 * express: the diffusion symmetry handler is a no-op, and the symmetry reflection
 * couples the velocity components, which a scalar diffusion scheme cannot assemble.
 */

namespace prism::eqn {

// forward declaration
class Transport;

class Momentum;

} // namespace prism::eqn

namespace prism::eqn::boundary {

template <typename Equation>
class IEquationBoundaryHandler : public prism::boundary::IBoundaryHandler {
  public:
    virtual auto name() const noexcept -> std::string = 0;
    virtual void apply(Equation& eqn, const mesh::BoundaryPatch& patch) = 0;
};

template <typename Equation>
class Symmetry : public IEquationBoundaryHandler<Equation> {
  public:
    auto name() const noexcept -> std::string override { return "symmetry"; }
    void apply(Equation& eqn, const mesh::BoundaryPatch& patch) override;
};

/** @brief Equation-level symmetry boundary handler for the momentum equation.
 *
 * Imposes the symmetry-plane momentum constraint of Moukalled et al. (2015),
 * Eqn (15.154), on the momentum matrix: at a symmetry face it adds the
 * wall-normal diagonal contribution g*ac (and the cross-component terms g*b to
 * the rhs), which the diffusion scheme cannot express because it couples the
 * velocity components. It is applied by Transport::updateCoeffs() after the
 * schemes; the diffusion scheme's symmetry handler is a no-op, so the two do
 * not double-count. */
template <>
class Symmetry<Transport> : public IEquationBoundaryHandler<Transport> {
  public:
    auto name() const noexcept -> std::string override { return "symmetry"; }
    void apply(Transport& eqn, const mesh::BoundaryPatch& patch) override;
};

} // namespace prism::eqn::boundary
