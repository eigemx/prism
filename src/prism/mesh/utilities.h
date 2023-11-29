#pragma once

#include <cassert>
#include <iostream>

#include "cell.h"
#include "face.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

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
    bool is_neighbor = !face.is_owned_by(cell.id());
    return face.area_vector() * std::pow(-1., static_cast<int>(is_neighbor));
}

/**
 * @brief Calculates the geometric weighting factor between two cells `c` and `n`
 * 
 * @param c First cell sharing face `f`.
 * @param n Second cell sharing face `f`.
 * @param f Face that is shared by the two cells.
 * @return double Geometric weighting factor between the two cells.
 */
auto inline geo_weight(const Cell& c, const Cell& n, const Face& f) -> double {
    auto gc = (n.center() - f.center()).norm() / (n.center() - c.center()).norm();
    assert(gc > 0 && "geo_weight() returned a negative weight factor");
    assert(gc <= 1.0 && "geo_weight() returned a weight factor higher than 1.0");
    return gc;
}

auto inline cells_volume_vec(const mesh::PMesh& mesh) -> VectorXd {
    VectorXd vec;
    vec.resize(mesh.n_cells());

    for (const auto& cell : mesh.cells()) {
        vec[cell.id()] = cell.volume();
    }

    return vec;
}

} // namespace prism::mesh