#pragma once

#include "prism/field/ifield.h"
#include "prism/mesh/cell.h"

namespace prism::scheme::convection {

template <field::IScalarBased Field>
auto phiAtDummyUpwind(const Field& field, const mesh::Cell& cell, const mesh::Cell& downwind)
    -> double;

template <field::IScalarBased Field>
auto phiTilde(const Field& field, const mesh::Cell& cell, const mesh::Cell& downwind) -> double;


template <field::IScalarBased Field>
auto phiAtDummyUpwind(const Field& field, const mesh::Cell& cell, const mesh::Cell& downwind)
    -> double {
    // The following is based on equation (12.66) from Moukalled et. al (2016)
    auto phi_downwind = field.valueAtCell(downwind);
    auto d_CD = downwind.center() - cell.center();
    return phi_downwind - (2 * field.gradAtCell(cell).dot(d_CD));
}

template <field::IScalarBased Field>
auto phiTilde(const Field& field, const mesh::Cell& cell, const mesh::Cell& downwind) -> double {
    // The following is based on equation (12.1) from Moukalled et. al (2016)
    const auto phi_upwind = phiAtDummyUpwind(field, cell, downwind);
    auto a = field.valueAtCell(cell) - phi_upwind;
    auto b = field.valueAtCell(downwind) - phi_upwind;
    return a / b;
}
} // namespace prism::scheme::convection
