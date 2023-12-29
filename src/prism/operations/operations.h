#pragma once

#include "prism/field.h"

namespace prism::ops {

// Calculates the divergence of a vector field U: ∇.U
auto div(const field::Vector& U, bool return_face_data = true) -> field::Scalar;

// Calculates the laplacian of a scalar field ϕ: ∇.∇ϕ
auto laplacian(const field::Scalar& phi, bool return_face_data = true) -> field::Scalar;

// Calculates the magnitude of a vector field U
auto mag(const field::Vector& U, bool return_face_data = true) -> field::Scalar;

// Calculates the curl of a vector field U
auto curl(const field::Vector& U, bool return_face_data = true) -> field::Vector;

} // namespace prism::ops
