#pragma once
#include <limits>

#include "gradient.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::gradient {

struct FaceContrib {
    bool is_interior;
    size_t other_id; // neighbor cell ID (interior) or face ID (boundary)
    Vector3d w_dir;  // weighted direction (dx*wk, dy*wk, dz*wk)
};

// Geometry-only data of the least-squares scheme. It depends only on the mesh, so it is
// built once and shared by every field's gradient scheme instance (the ops::grad/div
// helpers construct temporary fields, which would otherwise rebuild it every call).
struct LSQGeometry {
    // pseudoinverse matrices for each cell, used to compute the least-squares gradient.
    std::vector<Matrix3d> pinv_matrices;
    std::vector<std::vector<FaceContrib>> cell_face_cache;
};

class LeastSquares : public IGradient {
  public:
    explicit LeastSquares(const SharedPtr<mesh::PMesh>& mesh);

    auto gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d override;
    void computeAllCellGradients(const field::IScalar& field) override;

  protected:
    auto computeGradAtCell(const mesh::Cell& cell, const field::IScalar& field)
        -> Vector3d override;

  private:
    std::vector<Vector3d> _cell_gradients;
    SharedPtr<const LSQGeometry> _geometry;
    size_t _computed_event {std::numeric_limits<size_t>::max()};
};

} // namespace prism::gradient
