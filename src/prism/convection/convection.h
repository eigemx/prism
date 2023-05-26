#pragma once

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::convection {
class ConvectionSchemeBase {};

class Upwind : public FVScheme, public ConvectionSchemeBase {
  public:
    Upwind(double rho, VectorField& U, ScalarField& phi)
        : _rho(rho),
          _U(U),
          _phi(phi),
          _mesh(U.mesh()),
          _gradient_scheme(std::make_unique<gradient::GreenGauss>(U)) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    double _rho;
    VectorField& _U;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    std::unique_ptr<gradient::GradientSchemeBase> _gradient_scheme;
    bool _main_coeffs_calculated {false};
};
} // namespace prism::convection