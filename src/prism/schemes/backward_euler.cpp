#include "transient.h"

namespace prism::scheme::transient {
BackwardEuler::BackwardEuler(field::Scalar& rho, field::Scalar& phi, double dt)
    : _rho(rho),
      _phi(phi),
      _dt(dt),
      _phi_prev(phi),
      _rho_prev(_rho),
      FVScheme(phi.mesh().nCells()) {
    if (dt <= 0.0) {
        throw std::runtime_error(
            "BackwardEuler::BackwardEuler() constructor was given a zero or negative time step.");
    }
    _volume_field.resize(phi.mesh().nCells());

    for (const auto& cell : phi.mesh().cells()) {
        _volume_field[cell.id()] = cell.volume();
    }
}

void BackwardEuler::set_time_step(double dt) {
    if (dt <= 0.0) {
        throw std::runtime_error(
            "BackwardEuler::set_time_step() was given a zero or negative time step.");
    }
    _dt = dt;
}


void BackwardEuler::apply() {
    // TODO: Test this.
    // TODO: implement finalize() to set phi_prev = phi and same for rho.
    matrix().setIdentity();
    matrix().diagonal() = (_rho.values().array() * _volume_field.array()) / _dt;
    rhs() =
        (_rho_prev.values().array() * _volume_field.array() * _phi_prev.values().array()) / _dt;
}

} // namespace prism::scheme::transient