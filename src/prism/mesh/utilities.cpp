#include "utilities.h"

#include <cassert>
#include <cmath>

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
auto outwardAreaVector(const Face& face, const Cell& cell) -> Vector3d {
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
auto geometricWeight(const Cell& c, const Cell& n, const Face& face) -> double {
    // This is based on equation (6.30) in the book
    const Vector3d d_Cf = face.center() - c.center();
    const Vector3d d_fN = n.center() - face.center();
    const Vector3d Sf = outwardAreaVector(face, c);
    const Vector3d ef = Sf / Sf.norm();
    auto gc = d_Cf.dot(ef) / (d_Cf.dot(ef) + d_fN.dot(ef));
    // auto gc = (n.center() - face.center()).norm() / (n.center() - c.center()).norm();
    assert(gc > 0 && "geometricWeight() returned a negative weight factor");
    assert(gc <= 1.0 && "geometricWeight() returned a weight factor higher than 1.0");
    return gc;
}

} // namespace prism::mesh