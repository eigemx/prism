#pragma once

#include "prism/field.h"
#include "prism/fvscheme.h"
#include "prism/mesh/cell.h"
#include "prism/mesh/face.h"

namespace prism::source {

// Calculates the divergence of a vector field U: ∇.U
auto div(const VectorField& U) -> ScalarField;

// Calculates the gradient of a scalar field in the x-direction
template <typename GradientScheme>
auto grad_x(const ScalarField& phi) -> ScalarField;

// Calculates the gradient of a scalar field in the y-direction
template <typename GradientScheme>
auto grad_y(const ScalarField& phi) -> ScalarField;

// Calculates the gradient of a scalar field in the z-direction
template <typename GradientScheme>
auto grad_z(const ScalarField& phi) -> ScalarField;

// Calculates the gradient field of a scalar ϕ
template <typename GradientScheme>
auto grad(const ScalarField& phi) -> VectorField;

} // namespace prism::source
