#include "pressure.h"

namespace prism::field {
void PressureBHManagerSetter::set(IScalarBHManager& manager) {
    manager.addHandler<field::boundary::Fixed<Pressure>>();
    manager.addHandler<field::boundary::Empty<Pressure>>();
    manager.addHandler<field::boundary::Symmetry<Pressure>>();
    manager.addHandler<field::boundary::Outlet<Pressure>>();
    manager.addHandler<field::boundary::FixedGradient<Pressure>>();
    // this->boundaryHandlersManager().addHandler<field::boundary::NoSlip<Scalar>>();
}

} // namespace prism::field