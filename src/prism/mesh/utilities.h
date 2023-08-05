#pragma once

#include "cell.h"
#include "face.h"
#include "pmesh.h"

namespace prism::mesh {
/**
* @brief Calculates the area vector of a face, pointing out of the cell.
* The area vector is calculated as the area of the face times the unit normal
* vector of the face, pointing out of the cell.
* 
* @param face Face to calculate the area vector of.
* @param cell Cell that the face belongs to.
* @return Vector3d Area vector of the face.
*/
auto inline outward_area_vector(const Face& face, const Cell& cell) -> Vector3d {
    bool is_owned = face.is_owned_by(cell.id());
    auto neighbor_cell_id = is_owned ? face.neighbor().value() : face.owner();

    return face.area_vector() * std::pow(-1., static_cast<int>(is_owned));
}

/**
 * @brief Calculates the geometric weighting factor between two cells `c` and `n`
 * 
 * @param c Cell to calculate the weighting factor from.
 * @param n Cell to calculate the weighting factor to.
 * @param f Face that is shared by the two cells.
 * @return double Geometric weighting factor between the two cells.
 */
auto inline cells_weighting_factor(const Cell& c, const Cell& n, const Face& f) -> double {
    return (n.center() - f.center()).norm() / (n.center() - c.center()).norm();
}

} // namespace prism::mesh