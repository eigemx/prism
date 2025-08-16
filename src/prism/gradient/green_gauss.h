#pragma once

#include <vector>

#include "gradient.h"

namespace prism::gradient {

class GreenGauss : public IGradient {
  public:
    explicit GreenGauss(const SharedPtr<mesh::PMesh>& mesh);

    auto gradAtCell(const mesh::Cell& cell,  field::IScalar& field) -> Vector3d override;
    auto gradAtCellStored(const mesh::Cell& cell,
                          const field::IScalar& field) -> Vector3d override;

  private:
    auto correctSkewness(const mesh::Face& face,
                         const mesh::Cell& cell,
                         const mesh::Cell& nei) const -> double;

    std::vector<Vector3d> _cell_gradients;
};

} // namespace prism::gradient
