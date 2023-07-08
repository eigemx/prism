#pragma once

#include <memory>

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

// TODO: Implement a non-corrected version
class Diffusion : public FVScheme {
  public:
    // No gradient scheme is given, use Green-Gauss as a default explicit gradient scheme
    Diffusion(double kappa, ScalarField& phi)
        : _kappa(kappa),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::make_shared<gradient::GreenGauss>(phi)),
          FVScheme(phi.mesh().n_cells()) {}


    // Explicit gradient scheme is provided by the user
    Diffusion(double kappa,
              ScalarField& phi,
              std::shared_ptr<gradient::GradientSchemeBase> gradient_scheme)
        : _kappa(kappa),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::move(gradient_scheme)),
          FVScheme(phi.mesh().n_cells()) {}

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void correct_non_orhto_boundary_fixed(const mesh::Cell& cell,
                                          const mesh::Face& face,
                                          const Vector3d& T_f);

    void apply_boundary_gradient(const mesh::Cell& cell, const mesh::Face& face);

    double _kappa;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    std::shared_ptr<gradient::GradientSchemeBase> _gradient_scheme;
};

} // namespace prism::diffusion