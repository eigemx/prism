#pragma once


#include "prism/field/density.h"
#include "temporal.h"

namespace prism::scheme::temporal {
class AdamMoulton : public ITemporal {
  public:
    AdamMoulton(const SharedPtr<field::Scalar>& phi);

    AdamMoulton(const SharedPtr<field::Scalar>& phi, double dt);

    AdamMoulton(const SharedPtr<field::Density>& rho, const SharedPtr<field::Scalar>& phi);

    AdamMoulton(const SharedPtr<field::Density>& rho,
                const SharedPtr<field::Scalar>& phi,
                double dt);

    void apply() override;

  private:
    void applyIncompressible();
    void applyCompressible();
    void applyBackwardEuler();

    SharedPtr<field::Density> _rho;
};

} // namespace prism::scheme::temporal
