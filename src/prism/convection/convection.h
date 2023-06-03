#pragma once

#include "../field.h"
#include "../fvscheme.h"
#include "../mesh/pmesh.h"

namespace prism::convection {

class ConvectionSchemeBase {};

class CentralDifference : public FVScheme, public ConvectionSchemeBase {
  public:
    CentralDifference(double rho, VectorField& U, ScalarField& phi)
        : _rho(rho), _U(U), _phi(phi), _mesh(phi.mesh()) {}

    void inline finalize() override {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void apply_boundary_outlet(const mesh::Cell& cell, const mesh::Face& face);

    double _rho;
    VectorField& _U;
    ScalarField& _phi;
    mesh::PMesh _mesh;
};

class Upwind : public FVScheme, public ConvectionSchemeBase {
  public:
    Upwind(double rho, VectorField& U, ScalarField& phi)
        : _rho(rho), _U(U), _phi(phi), _mesh(phi.mesh()) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    double _rho;
    VectorField& _U;
    ScalarField& _phi;
    mesh::PMesh _mesh;
};

class SecondOrderUpwind : public FVScheme, public ConvectionSchemeBase {
  public:
    SecondOrderUpwind(double rho, VectorField& U, ScalarField& phi)
        : _rho(rho), _U(U), _phi(phi), _mesh(phi.mesh()) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    double _rho;
    VectorField& _U;
    ScalarField& _phi;
    mesh::PMesh _mesh;
};


} // namespace prism::convection