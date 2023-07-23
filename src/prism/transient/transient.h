#pragma once

#include "../field.h"
#include "../fvscheme.h"

namespace prism::transient {
class TransientBase {
  public:
    TransientBase(double dt) : _time_step(dt) {}
    virtual auto time_step() const -> double = 0;
    virtual auto time_step() -> double& { return _time_step; }

  private:
    double _time_step;
};

class ImplicitFirstOrder : public TransientBase, public FVScheme {
  public:
    ImplicitFirstOrder() = delete;
    ImplicitFirstOrder(const ScalarField& rho, ScalarField& phi, double dt)
        : FVScheme(phi.mesh().n_cells()), TransientBase(dt), _rho(rho), _phi(phi) {}

    auto inline field() -> ScalarField& override { return _phi; }

  private:
    // TODO: same as commented on source::ConstantScalar::apply_interior()
    // this can be vectorized in one function call
    // NOLINTNEXTLINE (`face` is unused, but required by interface)
    void apply_interior(const mesh::Cell& cell, const mesh::Face& face) override {
        auto cell_id = cell.id();
        auto phi = _phi[cell_id];
        auto rho = _rho[cell_id];
        auto volume = cell.volume();
        auto dt = time_step();

        coeff_matrix().coeffRef(cell_id, cell_id) = rho * volume / dt;
        rhs_vector()(cell_id) = rho * volume * phi / dt;
    }

    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    const ScalarField& _rho;
    ScalarField& _phi;
};
} // namespace prism::transient