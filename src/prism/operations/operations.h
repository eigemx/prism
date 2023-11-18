#pragma once

#include "prism/field.h"

namespace prism::ops {

// Calculates the divergence of a vector field U: ∇.U
auto div(const VectorField& U, bool return_face_data = true) -> ScalarField;

} // namespace prism::ops
