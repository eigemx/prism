#include "convection_hr.h"
#include "prism/field/scalar.h"

namespace prism::scheme::convection {
auto phiAtDummyUpwind(SharedPtr<field::Scalar>& field,
                      const mesh::Cell& cell,
                      const mesh::Cell& downwind) -> double {
    // The following is based on equation (12.66) from Moukalled et. al (2016)
    auto phi_downwind = field->valueAtCell(downwind);
    auto d_CD = downwind.center() - cell.center();
    return phi_downwind - (2 * field->gradAtCell(cell).dot(d_CD));
}

auto phiTilde(SharedPtr<field::Scalar>& field,
              const mesh::Cell& cell,
              const mesh::Cell& downwind) -> double {
    // The following is based on equation (12.1) from Moukalled et. al (2016)
    const auto phi_upwind = phiAtDummyUpwind(field, cell, downwind);
    auto a = field->valueAtCell(cell) - phi_upwind;
    auto b = field->valueAtCell(downwind) - phi_upwind;
    return a / b;
}
} // namespace prism::scheme::convection
