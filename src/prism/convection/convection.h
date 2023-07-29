#pragma once

#include <cmath>
#include <memory>

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"

namespace prism::convection {

// coefficients for the discretized convection equation for a face
struct CoeffsTriplet {
    double a_C {}; // cell
    double a_N {}; // neighbor
    double b {};   // source
};

// Finite volume scheme for the discretization of the convection term
class ConvectionBase : public FVScheme {
  public:
    // Gradient scheme is not given, using Green-Gauss as a default explicit gradient scheme
    ConvectionBase(ScalarField& rho, VectorField& U, ScalarField& phi)
        : _rho(rho),
          _U(U),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::make_shared<gradient::GreenGauss>(phi)),
          FVScheme(phi.mesh().n_cells()) {}

    // Gradient scheme is given
    ConvectionBase(ScalarField& rho,
                   VectorField& U,
                   ScalarField& phi,
                   std::shared_ptr<gradient::GradientSchemeBase> grad_scheme)
        : _rho(rho),
          _U(U),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::move(grad_scheme)),
          FVScheme(phi.mesh().n_cells()) {}

    auto inline field() -> ScalarField& override { return _phi; }

  private:
    virtual auto interpolate(double m_dot,
                             const mesh::Cell& cell,
                             const mesh::Cell& neighbor,
                             const mesh::Face& face,
                             const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet = 0;

    void apply_interior(const mesh::Cell& cellell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cellell, const mesh::Face& face) override;
    void apply_boundary_fixed(const mesh::Cell& cellell, const mesh::Face& face);
    void apply_boundary_outlet(const mesh::Cell& cellell, const mesh::Face& face);
    auto boundary_face_velocity(const mesh::Face& face) const -> Vector3d;

    ScalarField& _rho;
    VectorField& _U;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    std::shared_ptr<gradient::GradientSchemeBase> _gradient_scheme;
};

// Central difference scheme
class CentralDifference : public ConvectionBase {
  public:
    CentralDifference(ScalarField& rho, VectorField& U, ScalarField& phi)
        : ConvectionBase(rho, U, phi) {}

    CentralDifference(ScalarField& rho,
                      VectorField& U,
                      ScalarField& phi,
                      std::shared_ptr<gradient::GradientSchemeBase> grad_scheme)
        : ConvectionBase(rho, U, phi, std::move(grad_scheme)) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// Upwind scheme
class Upwind : public ConvectionBase {
  public:
    Upwind(ScalarField& rho, VectorField& U, ScalarField& phi) : ConvectionBase(rho, U, phi) {}

    Upwind(ScalarField& rho,
           VectorField& U,
           ScalarField& phi,
           std::shared_ptr<gradient::GradientSchemeBase> grad_scheme)
        : ConvectionBase(rho, U, phi, std::move(grad_scheme)) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// Second order upwind scheme
class SecondOrderUpwind : public ConvectionBase {
  public:
    SecondOrderUpwind(ScalarField& rho, VectorField& U, ScalarField& phi)
        : ConvectionBase(rho, U, phi) {}

    SecondOrderUpwind(ScalarField& rho,
                      VectorField& U,
                      ScalarField& phi,
                      std::shared_ptr<gradient::GradientSchemeBase> grad_scheme)
        : ConvectionBase(rho, U, phi, std::move(grad_scheme)) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};

// QUICK scheme
class QUICK : public ConvectionBase {
  public:
    QUICK(ScalarField& rho, VectorField& U, ScalarField& phi) : ConvectionBase(rho, U, phi) {}

    QUICK(ScalarField& rho,
          VectorField& U,
          ScalarField& phi,
          std::shared_ptr<gradient::GradientSchemeBase> grad_scheme)
        : ConvectionBase(rho, U, phi, std::move(grad_scheme)) {}

    auto interpolate(double m_dot,
                     const mesh::Cell& cell,
                     const mesh::Cell& neighbor,
                     const mesh::Face& face,
                     const std::shared_ptr<gradient::GradientSchemeBase>& grad_scheme)
        -> CoeffsTriplet override;
};


auto inline face_mass_flow_rate(double rho, const Vector3d& U, const Vector3d& S) -> double {
    return rho * U.dot(S);
}

} // namespace prism::convection
