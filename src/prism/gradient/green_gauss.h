#pragma once

#include <limits>
#include <vector>

#include "gradient.h"

namespace prism::gradient {

class GreenGauss : public IGradient {
  public:
    explicit GreenGauss(const SharedPtr<mesh::PMesh>& mesh);

    auto gradAtCell(const mesh::Cell& cell, const field::IScalar& field) -> Vector3d override;
    void computeAllCellGradients(const field::IScalar& field) override;

  protected:
    auto computeGradAtCell(const mesh::Cell& cell, const field::IScalar& field)
        -> Vector3d override;

  private:
    auto correctSkewness(const mesh::Face& face,
                         const mesh::Cell& cell,
                         const mesh::Cell& nei) const -> double;

    std::vector<Vector3d> _cell_gradients;
    size_t _computed_event {std::numeric_limits<size_t>::max()};
};

} // namespace prism::gradient
