#pragma once

#include "../field.h"
#include "../fvscheme.h"

namespace prism::convection {

class ConvectionSchemeBase {};

class CentralDifference : public FVScheme, public ConvectionSchemeBase {
  public:
    CentralDifference(VectorField& U, ScalarField& phi) : _U(U), _phi(phi) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    VectorField& _U;
    ScalarField& _phi;
};

class Upwind : public FVScheme, public ConvectionSchemeBase {
  public:
    Upwind(VectorField& U, ScalarField& phi) : _U(U), _phi(phi) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    VectorField& _U;
    ScalarField& _phi;
};

class SecondOrderUpwind : public FVScheme, public ConvectionSchemeBase {
  public:
    SecondOrderUpwind(VectorField& U, ScalarField& phi) : _U(U), _phi(phi) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    VectorField& _U;
    ScalarField& _phi;
};


} // namespace prism::convection