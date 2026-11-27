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

    /** @brief Skewness correction added to the linearly interpolated value of an interior face.
     *
     * The plain Green-Gauss scheme applies no correction. Derived schemes override this. */
    virtual auto skewnessCorrection(const mesh::Face& face,
                                    const mesh::Cell& cell,
                                    const mesh::Cell& nei) const -> f64;

    std::vector<Vector3d> _cell_gradients;
    size_t _computed_event {std::numeric_limits<size_t>::max()};
};

/** @brief Green-Gauss gradient scheme with skewness correction.
 *
 * The correction uses the cell gradients from the previous gradient sweep. */
class SkewCorrectedGreenGauss : public GreenGauss {
  public:
    using GreenGauss::GreenGauss;

  protected:
    auto skewnessCorrection(const mesh::Face& face,
                            const mesh::Cell& cell,
                            const mesh::Cell& nei) const -> f64 override;
};

} // namespace prism::gradient
