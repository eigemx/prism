#pragma once

#include <vector>

#include "gradient.h"

namespace prism::gradient {

class GreenGauss : public IGradient {
  public:
    GreenGauss(field::IScalar* field);

    auto gradAtCell(const mesh::Cell& cell) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell) -> Vector3d override;

  private:
    auto correctSkewness(const mesh::Face& face,
                         const mesh::Cell& cell,
                         const mesh::Cell& nei) const -> double;


    auto boundaryFaceIntegral(const mesh::Face& f) -> Vector3d;

    std::vector<Vector3d> _cell_gradients;
};

} // namespace prism::gradient
