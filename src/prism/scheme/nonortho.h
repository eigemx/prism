#pragma once

#include "Eigen/Core"
#include "prism/constants.h"
#include "prism/types.h"

namespace prism::nonortho {

struct NonOrthoPair {
    Vector3d Ef, Tf;
};

class INonOrthoCorrector {
  public:
    INonOrthoCorrector() = default;
    virtual ~INonOrthoCorrector() = default;
    INonOrthoCorrector(const INonOrthoCorrector&) = default;
    INonOrthoCorrector(INonOrthoCorrector&&) noexcept = default;
    auto operator=(const INonOrthoCorrector&) -> INonOrthoCorrector& = default;
    auto operator=(INonOrthoCorrector&&) noexcept -> INonOrthoCorrector& = default;

    virtual auto decompose(const Vector3d& Sf, const Vector3d& e) const -> NonOrthoPair = 0;
};

class OverRelaxedCorrector : public INonOrthoCorrector {
  public:
    auto decompose(const Vector3d& Sf, const Vector3d& e) const -> NonOrthoPair override;
};

class MinimumCorrector : public INonOrthoCorrector {
  public:
    auto decompose(const Vector3d& Sf, const Vector3d& e) const -> NonOrthoPair override;
};

class OrthogonalCorrector : public INonOrthoCorrector {
  public:
    auto decompose(const Vector3d& Sf, const Vector3d& e) const -> NonOrthoPair override;
};

} // namespace prism::nonortho