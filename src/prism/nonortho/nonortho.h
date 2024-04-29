#pragma once

#include "prism/constants.h"
#include "prism/field.h"
#include "prism/gradient/gradient.h"
#include "prism/mesh/pmesh.h"
#include "prism/types.h"

namespace prism::nonortho {

struct NonOrthoTriplet {
    Vector3d Sf {}, Ef {}, Tf {};
};

class AbstractCorrector {
  public:
    AbstractCorrector() = default;
    virtual ~AbstractCorrector() = default;
    AbstractCorrector(const AbstractCorrector&) = default;
    AbstractCorrector(AbstractCorrector&&) noexcept = default;
    auto operator=(const AbstractCorrector&) -> AbstractCorrector& = default;
    auto operator=(AbstractCorrector&&) noexcept -> AbstractCorrector& = default;

    virtual auto interior_triplet(const mesh::Cell& owner,
                                  const mesh::Cell& neighbor,
                                  const mesh::Face& face) const -> NonOrthoTriplet;

    virtual auto boundary_triplet(const mesh::Cell& owner, const mesh::Face& face) const
        -> NonOrthoTriplet;

  private:
    virtual auto decompose(const Vector3d& Sf, const Vector3d& e) const -> Vector3d = 0;
};


class OverRelaxedCorrector : public AbstractCorrector {
  private:
    auto decompose(const Vector3d& Sf, const Vector3d& e) const -> Vector3d override;
};

class MinimumCorrector : public AbstractCorrector {
  private:
    auto decompose(const Vector3d& Sf, const Vector3d& e) const -> Vector3d override;
};

class OrthogonalCorrector : public AbstractCorrector {
  private:
    auto decompose(const Vector3d& Sf, const Vector3d& e) const -> Vector3d override;
};

auto inline AbstractCorrector::interior_triplet(const mesh::Cell& owner,
                                                const mesh::Cell& neighbor,
                                                const mesh::Face& face) const -> NonOrthoTriplet {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the two cells
    auto d_CF = neighbor.center() - owner.center();
    auto d_CF_norm = d_CF.norm();

    // unit vector in d_CF direction
    auto e = d_CF / d_CF_norm;

    auto Ef = decompose(Sf, e);

    return {Sf, Ef, Sf - Ef};
}

auto inline AbstractCorrector::boundary_triplet(const mesh::Cell& owner,
                                                const mesh::Face& face) const -> NonOrthoTriplet {
    const auto& Sf = face.area_vector();

    // vector joining the centers of the cell and the face
    auto d_Cf = face.center() - owner.center();
    auto d_Cf_norm = d_Cf.norm();
    auto e = d_Cf / d_Cf_norm;

    auto Ef = decompose(Sf, e);

    return {Sf, Ef, Sf - Ef};
}

auto inline OverRelaxedCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const
    -> Vector3d {
    return ((Sf.dot(Sf) / (e.dot(Sf) + EPSILON))) * e;
}

auto inline MinimumCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const -> Vector3d {
    return (e.dot(Sf)) * e;
}

auto inline OrthogonalCorrector::decompose(const Vector3d& Sf, const Vector3d& e) const
    -> Vector3d {
    return Sf.norm() * e;
}


} // namespace prism::nonortho