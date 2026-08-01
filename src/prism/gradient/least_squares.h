#pragma once
#include <limits>

#include "gradient.h"
#include "prism/mesh/pmesh.h"

namespace prism::gradient {

class LeastSquares : public IGradient {
  public:
    explicit LeastSquares(const SharedPtr<mesh::PMesh>& mesh);

    auto gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d override;
    void computeAllCellGradients(const field::IScalar& field) override;

  protected:
    auto computeGradAtCell(const mesh::Cell& cell, const field::IScalar& field)
        -> Vector3d override;

  private:
    struct FaceContrib {
        bool is_interior;
        size_t other_id; // neighbor cell ID (interior) or face ID (boundary)
        Vector3d w_dir;  // weighted direction (dx*wk, dy*wk, dz*wk)
    };

    void setPseudoInvMatrices(const SharedPtr<mesh::PMesh>& mesh);

    std::vector<Vector3d> _cell_gradients;
    std::vector<Matrix3d> _pinv_matrices;
    std::vector<std::vector<FaceContrib>> _cell_face_cache;
    size_t _computed_event {std::numeric_limits<size_t>::max()};
};

} // namespace prism::gradient
