#pragma once

#include "prism/field.h"

namespace prism::ops {

// Calculates the divergence of a vector field U: ∇.U
auto div(const VectorField& U, bool return_face_data = true) -> ScalarField;

// Calculates the laplacian of a scalar field ϕ: ∇.∇ϕ
auto laplacian(const ScalarField& phi, bool return_face_data = true) -> ScalarField;

// Calculates the magnitude of a vector field U
auto mag(const VectorField& U, bool return_face_data = true) -> ScalarField;

// Calculates the curl of a vector field U
auto curl(const VectorField& U, bool return_face_data = true) -> VectorField;

} // namespace prism::ops
