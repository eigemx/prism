#pragma once

#include "fvscheme.h"

namespace prism::transient {
class AbstractTransientScheme {};

class BackwardEuler : public FVScheme, public AbstractTransientScheme {
  public:
    BackwardEuler(field::Scalar& rho, field::Scalar& phi, double dt);

    void apply() override;
    auto inline field() -> std::optional<field::Scalar> override { return _phi; }
    auto inline time_step() const -> double { return _dt; }
    void set_time_step(double dt);

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Face& face) override {}

    field::Scalar _rho;
    field::Scalar _phi;
    field::Scalar _rho_prev;
    field::Scalar _phi_prev;
    VectorXd _volume_field;
    double _dt {1e-8};
};

} // namespace prism::transient
