#pragma once

#include "prism/scheme/scheme.h"

namespace prism::scheme::source {
// source term is assumed to be always on the right hand side of the conserved equation.
// we can control the sign of the source term by using the Sign template parameter.

class ISource {};

class IImplicitSource : public ISource {};

template <typename T>
concept IImplicitSourceBased = std::derived_from<T, IImplicitSource>;

// We inherit from FVScheme with field::Scalar as template specialization because for explicit
// sources the type of the field won't matter, because we are not contributing to the matrix of
// coefficients for the linear system of the conserved equation.
class IExplicitSource : public ISource, public IPartialScheme {
  public:
    IExplicitSource(std::size_t n_cells) : IPartialScheme(n_cells) {}
};

template <typename T>
concept IExplicitSourceBased = std::derived_from<T, IExplicitSource>;


} // namespace prism::scheme::source
