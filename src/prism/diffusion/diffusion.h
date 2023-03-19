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
    Linear(double kappa, ScalarField& phi);

    template <typename GradientScheme>
    Linear(double kappa, ScalarField& phi, const GradientScheme& gradient_scheme);

    void finalize() override;
    inline auto mesh() const -> const mesh::PMesh& override { return _mesh; }

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

    double _kappa;
    ScalarField& _phi;
    const mesh::PMesh& _mesh;
    std::unique_ptr<gradient::GradientSchemeBase> _gradient_scheme;
};

} // namespace prism::diffusion