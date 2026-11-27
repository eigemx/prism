#pragma once

#include "prism/field/density.h"
#include "temporal.h"

namespace prism::scheme::temporal {

namespace detail {

/** @brief Assembles the incompressible backward-Euler diagonal and RHS into `system`.
 *
 * This lives in detail so AdamMoulton can reuse the Backward-Euler fallback without
 * constructing a temporary scheme and copying its matrix. */
void assembleBackwardEuler(LinearSystem& system,
                           const SharedPtr<field::Scalar>& field,
                           f64 dt,
                           const VectorXd& phi_prev);

} // namespace detail

class BackwardEuler : public ITemporal {
  public:
    BackwardEuler(const SharedPtr<field::Scalar>& phi);

    BackwardEuler(const SharedPtr<field::Scalar>& phi, double dt);

    BackwardEuler(const SharedPtr<field::Density>& rho, const SharedPtr<field::Scalar>& phi);

    BackwardEuler(const SharedPtr<field::Density>& rho,
                  const SharedPtr<field::Scalar>& phi,
                  double dt);

    void apply() override;

  private:
    void applyIncompressible();
    void applyCompressible();

    SharedPtr<field::Density> _rho;
};
} // namespace prism::scheme::temporal
