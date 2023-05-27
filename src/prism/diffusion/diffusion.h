#pragma once

#include <memory>

#include "../field.h"
#include "../fvscheme.h"
#include "../gradient/gradient.h"
#include "../mesh/pmesh.h"
#include "../types.h"

namespace prism::diffusion {

class DiffusionSchemeBase {};

class Linear : public FVScheme, public DiffusionSchemeBase {
  public:
    // Default explicit gradient scheme is Green-Gauss
    Linear(double kappa, ScalarField& phi)
        : _kappa(kappa),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::make_unique<gradient::GreenGauss>(phi)) {}

    // Explicit gradient scheme is provided by the user
    template <typename G,
              typename = std::enable_if_t<std::is_base_of_v<gradient::GradientSchemeBase, G>>>
    Linear(double kappa, ScalarField& phi, const G& gradient_scheme)
        : _kappa(kappa),
          _phi(phi),
          _mesh(phi.mesh()),
          _gradient_scheme(std::make_unique<G>(gradient_scheme)) {}

    inline void finalize() override { _main_coeffs_calculated = true; }

  private:
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override;
    void correct_non_orhto_interior(const mesh::Cell& cell,
                                    const mesh::Cell& nei_cell,
                                    const mesh::Face& face,
                                    const Vector3d& T_f);

    void apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override;

    void apply_boundary_fixed(const mesh::Cell& cell, const mesh::Face& face);
    void correct_non_orhto_boundary_fixed(const mesh::Cell& cell,
                                          const mesh::Face& face,
                                          const Vector3d& T_f);

    void apply_boundary_gradient(const mesh::Cell& cell, const mesh::Face& face);

    double _kappa;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    std::unique_ptr<gradient::GradientSchemeBase> _gradient_scheme;
    bool _main_coeffs_calculated {false};
};

} // namespace prism::diffusion