#pragma once


#include "cell.h"
#include "face.h"
#include "prism/types.h"

namespace prism::mesh {
// TODO: move definitions to .cpp file
/**
 * @brief Calculates the area vector of a face, pointing out of the cell.
 * The area vector is calculated as the area of the face times the unit normal
 * vector of the face, pointing out of the cell.
 *
 * @param face Face to calculate the area vector of.
 * @param cell Cell that the face belongs to.
 * @return Vector3d Area vector of the face.
 */
auto outwardAreaVector(const Face& face, const Cell& cell) -> Vector3d;

/**
 * @brief Calculates the geometric weighting factor between two cells `c` and `n`
 *
 * @param c First cell sharing face `f`.
 * @param n Second cell sharing face `f`.
 * @param f Face that is shared by the two cells.
 * @return double Geometric weighting factor between the two cells.
 */
auto geometricWeight(const Cell& c, const Cell& n, const Face& f) -> double;
} // namespace prism::mesh