#pragma once

#include "fvscheme.h"

namespace prism::transient {
class AbstractTransientScheme {};

class BackwardEuler : public FVScheme, public AbstractTransientScheme {
  public:
    BackwardEuler(ScalarField& rho, ScalarField& phi, double dt);

    void apply() override;
    auto inline field() -> std::optional<ScalarField> override { return _phi; }
    auto inline time_step() const -> double { return _dt; }
    void set_time_step(double dt);

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    ScalarField _rho;
    ScalarField _phi;
    ScalarField _rho_prev;
    ScalarField _phi_prev;
    VectorXd _volume_field;
    double _dt {1e-8};
};

} // namespace prism::transient
