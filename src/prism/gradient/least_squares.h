#pragma once
#include "gradient.h"
#include "prism/mesh/pmesh.h"

namespace prism::gradient {

class LeastSquares : public IGradient {
  public:
    explicit LeastSquares(const SharedPtr<mesh::PMesh>& mesh);

    auto gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell,
                          const field::IScalar& field) -> Vector3d override;

  private:
    void setPseudoInvMatrices(const SharedPtr<mesh::PMesh>& mesh);

    std::vector<Vector3d> _cell_gradients;
    std::vector<Matrix3d> _pinv_matrices; // pseudo-inverse matrices
};

} // namespace prism::gradient
