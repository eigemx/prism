#pragma once

// TODO: rename this to constant.h and move ImplicitPhi to somewher else

#include "fvscheme.h"
#include "prism/field.h"

namespace prism::source {
// source term is assumed to be always on the right hand side of the conserved equation.
// we can control the sign of the source term by using the SourceSign template parameter.
enum class SourceSign { Positive, Negative };

// Discretized constant source/sink term (like gravity), takes a scalar field
// and adds it to the right hand side of the system of equation
template <SourceSign Sign = SourceSign::Positive>
class ConstantScalar : public FVScheme {
  public:
    ConstantScalar(ScalarField& phi);
    ConstantScalar(ScalarField&& phi) = delete;

    auto requires_correction() const -> bool override { return false; }

    void apply() override;

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    ScalarField _phi;
    VectorXd _volume_field;
};

// TODO: Test this!
template <SourceSign Sign = SourceSign::Positive>
class ImplicitPhi : public FVScheme {
  public:
    ImplicitPhi(ScalarField& phi, double coeff = 1.0)
        : _phi(phi), _coeff(coeff), FVScheme(phi.mesh().n_cells()) {}
    ImplicitPhi(ScalarField&& phi) = delete;

    void apply() override;
    auto inline field() -> std::optional<ScalarField> override { return _phi; }

  private:
    void inline apply_interior(const mesh::Face& face) override {}
    void inline apply_boundary(const mesh::Cell& cell, const mesh::Face& face) override {}

    auto inline requires_correction() const -> bool override { return false; }

    ScalarField _phi;
    double _coeff {1.0};
};

template <SourceSign Sign>
ConstantScalar<Sign>::ConstantScalar(ScalarField& phi)
    : _phi(phi), FVScheme(phi.mesh().n_cells()) {
    _volume_field.resize(phi.mesh().n_cells());

    for (const auto& cell : phi.mesh().cells()) {
        _volume_field[cell.id()] = cell.volume();
    }
}

template <>
void inline ConstantScalar<SourceSign::Positive>::apply() {
    rhs() = _phi.data().array() * _volume_field.array();
}

template <>
void inline ConstantScalar<SourceSign::Negative>::apply() {
    rhs() = -_phi.data().array() * _volume_field.array();
}

template <>
void inline ImplicitPhi<SourceSign::Positive>::apply() {
    matrix().setIdentity();
    matrix() *= -_coeff;
}

template <>
void inline ImplicitPhi<SourceSign::Negative>::apply() {
    matrix().setIdentity();
    matrix() *= _coeff;
}

} // namespace prism::source