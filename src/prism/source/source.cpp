#include "source.h"

namespace prism::source {

template <>
void ConstantScalar<SourceSign::Positive>::apply() {
    rhs() = _phi.data().array() * _volume_field.array();
}

template <>
void ConstantScalar<SourceSign::Negative>::apply() {
    rhs() = -_phi.data().array() * _volume_field.array();
}

template <>
void ImplicitPhi<SourceSign::Positive>::apply() {
    matrix().setIdentity();
    matrix() *= -_coeff;
}

template <>
void ImplicitPhi<SourceSign::Negative>::apply() {
    matrix().setIdentity();
    matrix() *= _coeff;
}

} // namespace prism::source