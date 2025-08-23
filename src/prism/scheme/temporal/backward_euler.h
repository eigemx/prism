#pragma once

#include "prism/field/density.h"
#include "temporal.h"

namespace prism::scheme::temporal {

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
