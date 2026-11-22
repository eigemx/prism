#include "least_squares.h"

#include <unordered_map>

#include "prism/constants.h"

namespace prism::gradient {
namespace {
using CacheEntry = std::pair<std::weak_ptr<mesh::PMesh>, std::shared_ptr<const LSQGeometry>>;
std::unordered_map<const mesh::PMesh*, CacheEntry> g_geometry_cache;

auto buildGeometry(const SharedPtr<mesh::PMesh>& mesh) -> std::shared_ptr<const LSQGeometry> {
    // This function is based on section 9.3 'Least-Square Gradient'
    auto geometry = std::make_shared<LSQGeometry>();
    geometry->pinv_matrices.resize(mesh->cellCount());
    geometry->cell_face_cache.resize(mesh->cellCount());

    for (const auto& cell : mesh->cells()) {
        // A 3x3 matrix of the left hand side of equation (9.27)
        // for the k-th cell, we calculate this distance matrix D
        // and push the pseudo-inverse [(D * D^T)^{-1} * D^T] to _pinv_matrix vector
        Matrix3d d_matrix = Matrix3d::Zero();

        for (auto face_id : cell.facesIds()) {
            const auto& face = mesh->face(face_id);

            // This will hold the distance vector from neighbor cell center to k-th cell
            // center, or in case we have a boundary face, r_CF will be the distance vector
            // from boundary face center to the k-th cell center.
            // check equation (9.22)
            Vector3d r_CF = {.0, .0, .0};

            if (face.isInterior()) {
                // interior face
                const auto neighbor = mesh->otherSharingCell(cell, face);
                r_CF = neighbor.center() - cell.center();

                // weight factor defined in equation (9.28)
                const double wk = 1 / (r_CF.norm() + EPSILON);
                Vector3d w_dir = r_CF * wk;
                d_matrix += w_dir * r_CF.transpose();
                geometry->cell_face_cache[cell.id()].push_back({true, neighbor.id(), w_dir});

            } else {
                // boundary face
                r_CF = face.center() - cell.center();
                const f64 wk = 1 / (r_CF.norm() + EPSILON);
                Vector3d w_dir = r_CF * wk;
                d_matrix += w_dir * r_CF.transpose();

                const auto& patch = mesh->faceBoundaryPatch(face);
                if (!patch.isEmpty()) {
                    geometry->cell_face_cache[cell.id()].push_back({false, face_id, w_dir});
                }
            }
        }

        const auto& d_matrix_t = d_matrix.transpose();
        geometry->pinv_matrices[cell.id()] = (d_matrix_t * d_matrix).inverse() * d_matrix_t;
    }
    return geometry;
}

auto sharedGeometry(const SharedPtr<mesh::PMesh>& mesh) -> std::shared_ptr<const LSQGeometry> {
    const auto* key = mesh.get();
    const auto it = g_geometry_cache.find(key);

    if (it != g_geometry_cache.end() && !it->second.first.expired()) {
        return it->second.second;
    }

    auto geometry = buildGeometry(mesh);
    g_geometry_cache[key] = {mesh, geometry};

    return geometry;
}
} // namespace

LeastSquares::LeastSquares(const SharedPtr<mesh::PMesh>& mesh) {
    _cell_gradients.resize(mesh->cellCount(), Vector3d::Zero());
    _geometry = sharedGeometry(mesh);
}

auto LeastSquares::gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d {
    // Lazy cache: recompute over the whole mesh only when the field's cell values changed.
    if (field.eventNo() != _computed_event) {
        computeAllCellGradients(field);
    }
    return _cell_gradients[cell.id()];
}

void LeastSquares::computeAllCellGradients(const field::IScalar& field) {
    for (const auto& cell : field.mesh()->cells()) {
        computeGradAtCell(cell, field);
    }
    _computed_event = field.eventNo();
}

auto LeastSquares::computeGradAtCell(const mesh::Cell& cell, const field::IScalar& field)
    -> Vector3d {
    // right hand side of equation (9.27)
    Vector3d b {0.0, 0.0, 0.0};
    auto phi_cell = field.valueAtCell(cell);

    for (const auto& gf : _geometry->cell_face_cache[cell.id()]) {
        f64 other_phi {};
        if (gf.is_interior) {
            other_phi = field.valueAtCell(gf.other_id);
        } else {
            other_phi = field.valueAtFace(gf.other_id);
        }
        // w_dir is the weighted direction vector (dx*wk, dy*wk, dz*wk)
        // from equation (9.24) and (9.28)
        b += gf.w_dir * (other_phi - phi_cell);
    }

    Vector3d grad = _geometry->pinv_matrices[cell.id()] * b;
    _cell_gradients[cell.id()] = grad;

    return grad;
}
} // namespace prism::gradient
