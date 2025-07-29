#pragma once
#include "gradient.h"

namespace prism::gradient {

class LeastSquares : public IGradient {
  public:
    LeastSquares(field::IScalar* field);

    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) -> Vector3d override;

  private:
    void setPseudoInvMatrices();

    std::vector<Vector3d> _cell_gradients;
    std::vector<Matrix3d> _pinv_matrices; // pseudo-inverse matrices
};

} // namespace prism::gradient
