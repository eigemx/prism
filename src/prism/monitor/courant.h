#pragma once

#include "prism/field/ifield.h"
#include "prism/field/scalar.h"
#include "prism/types.h"

namespace prism::monitor {

/** @brief Computes a cell Courant number field from a face flux field.
 *
 * Each cell gets Co = 0.5 * dt * sum_faces(|mdot_face|) / V_cell, where the sum runs over
 * all faces of the cell. The 0.5 factor accounts for each interior face contributing its
 * flux to both neighbouring cells.
 *
 * @param mdot The face flux field.
 * @param dt The time step size.
 * @return A scalar field holding the Courant number of each cell. */
auto courantNumber(const field::IVector& mdot, f64 dt) -> field::Scalar;

} // namespace prism::monitor
