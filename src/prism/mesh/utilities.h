#pragma once

#include <cassert>

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
auto inline outwardAreaVector(const Face& face, const Cell& cell) -> Vector3d {
    bool is_neighbor = !face.isOwnedBy(cell.id());
    return face.areaVector() * std::pow(-1., static_cast<int>(is_neighbor));
}

/**
 * @brief Calculates the geometric weighting factor between two cells `c` and `n`
 *
 * @param c First cell sharing face `f`.
 * @param n Second cell sharing face `f`.
 * @param f Face that is shared by the two cells.
 * @return double Geometric weighting factor between the two cells.
 */
auto inline geometricWeight(const Cell& c, const Cell& n, const Face& f) -> double {
    auto gc = (n.center() - f.center()).norm() / (n.center() - c.center()).norm();
    assert(gc > 0 && "geometricWeight() returned a negative weight factor");
    assert(gc <= 1.0 && "geometricWeight() returned a weight factor higher than 1.0");
    return gc;
}

} // namespace prism::mesh